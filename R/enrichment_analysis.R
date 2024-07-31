utils::globalVariables(c("ENTREZID"))
#' Perform over-representation analysis
#'
#' `do_ora()` performs over-representation analysis (ORA) using the clusterProfiler package.
#' It also produces useful plots to visualize the results.
#'
#' @param gene_list A character vector containing the gene names.
#' @param database The database to perform the ORA. It can be either "KEGG", "GO", or "Reactome".
#' @param pval_lim The p-value threshold to consider a term as significant.
#' @param ncateg The number of categories to show in the plots.
#' @param fontsize The font size for the plots.
#'
#' @return A list containing the results of the ORA.
#' @export
#'
#' @details The ontology option used when database = "GO" is "BP" (Biological Process).
#' @examples
#' # Perform Differential Expression Analysis
#' de_res <- do_limma(example_data, example_metadata, "AML", wide = FALSE)
#'
#' # Extract the up-regulated proteins for AML
#' sig_up_proteins_aml <- de_res$de_results$AML |>
#'   dplyr::filter(sig == "significant up") |>
#'   dplyr::pull(Assay)
#'
#' # Perform ORA with GO database
#' do_ora(sig_up_proteins_aml, database = "GO", ncateg = 5)
do_ora <- function(gene_list,
                   database = c("KEGG", "GO", "Reactome"),
                   pval_lim = 0.05,
                   ncateg = 10,
                   fontsize = 10) {
  database <- match.arg(database)

  # From gene name to ENTREZID
  gene_conversion <- clusterProfiler::bitr(gene_list,
                                           fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Hs.eg.db::org.Hs.eg.db)

  gene_list <- gene_conversion |> dplyr::pull(ENTREZID) |> unique()

  if (database == "KEGG") {
    # Perform KEGG enrichment analysis
    enrichment <- clusterProfiler::enrichKEGG(gene = gene_list, organism = "hsa")
  } else if (database == "GO") {
    # Perform GO enrichment analysis
    enrichment <- clusterProfiler::enrichGO(gene = gene_list,
                                            OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                            ont = "BP")
  } else if (database == "Reactome") {
    # Perform Reactome enrichment analysis
    enrichment <- ReactomePA::enrichPathway(gene = gene_list,
                                            organism = "human",
                                            pvalueCutoff = pval_lim)
  }

  if (!any(enrichment@result$p.adjust < pval_lim)) {
    message("No significant terms found.")
    return(NULL)
  }

  # Visualize the results
  dotplot <- clusterProfiler::dotplot(enrichment,
                                      showCategory = ncateg,
                                      font.size = fontsize)

  barplot <- barplot(enrichment,
                     drop = TRUE,
                     showCategory = ncateg,
                     font.size = fontsize)

  goplot <- clusterProfiler::goplot(enrichment,
                                    showCategory = ncateg)

  enrichment <- clusterProfiler::setReadable(enrichment, OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  cnetplot <- clusterProfiler::cnetplot(enrichment,
                                        showCategory = ncateg,
                                        categorySize = "pvalue",
                                        color.params = list(foldChange = gene_list),
                                        cex.params = list(category_label = (fontsize + 2)/12,
                                                          gene_label = (fontsize)/12))

  return(list("enrichment" = enrichment,
              "dotplot" = dotplot,
              "barplot" = barplot,
              "goplot" = goplot,
              "cnetplot" = cnetplot))
}


#' Perform gene set enrichment analysis
#'
#' This function performs gene set enrichment analysis (GSEA) using the clusterProfiler package.
#' It also produces useful plots to visualize the results.
#'
#' @param de_results A tibble containing the results of a differential expression analysis.
#' @param database The database to perform the GSEA. It can be either "KEGG", "GO", or "Reactome".
#' @param background A character vector containing the background genes.
#' @param pval_lim The p-value threshold to consider a term as significant.
#' @param ncateg The number of categories to show in the plots.
#' @param fontsize The font size for the plots.
#'
#' @return A list containing the results of the GSEA.
#' @export
#'
#' @details The ontology option used when database = "GO" is "ALL".
#' When Reactome is used, background functionality is not available.
#' @examples
#' # Run Differential Expression Analysis and extract results
#' de_res <- do_limma(example_data, example_metadata, "AML", wide = FALSE)
#' de_results <- de_res$de_results
#'
#' # Run GSEA with Reactome database
#' do_gsea(de_results,
#'         database = "GO",
#'         pval_lim = 0.9,  # Remember that the data is artificial
#'         ncateg = 7,
#'         fontsize = 7)
do_gsea <- function(de_results,
                    database = c("KEGG", "GO", "Reactome"),
                    background = NULL,
                    pval_lim = 0.05,
                    ncateg = 10,
                    fontsize = 10) {
  database <- match.arg(database)

  # Prepare sorted_protein_list
  protein_list <- stats::setNames(de_results$logFC,
                                  de_results$Assay)
  sorted_protein_list <- sort(protein_list, decreasing = TRUE)

  # From gene name to ENTREZID
  protein_conversion <- clusterProfiler::bitr(names(sorted_protein_list),
                                              fromType = "SYMBOL",
                                              toType = "ENTREZID",
                                              OrgDb = org.Hs.eg.db::org.Hs.eg.db)

  protein_list <- stats::setNames(sorted_protein_list, protein_conversion$ENTREZID)

  if (database == "KEGG") {
    # Perform GSEA for KEGG
    if (!is.null(background)) {
      enrichment <- clusterProfiler::gseKEGG(geneList = protein_list,
                                             organism = "hsa",
                                             pvalueCutoff = pval_lim,
                                             pAdjustMethod = "BH",
                                             universe = background,
                                             minGSSize = 10,
                                             maxGSSize = 500)
    } else {
      enrichment <- clusterProfiler::gseKEGG(geneList = protein_list,
                                             organism = "hsa",
                                             pvalueCutoff = pval_lim,
                                             pAdjustMethod = "BH",
                                             minGSSize = 10,
                                             maxGSSize = 500)
    }
  } else if (database == "GO") {
    # Perform GSEA for GO
    if (!is.null(background)) {
      enrichment <- clusterProfiler::gseGO(geneList = protein_list,
                                         OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                         ont = "BP",
                                         pvalueCutoff = pval_lim,
                                         pAdjustMethod = "BH",
                                         universe = background,
                                         minGSSize = 10,
                                         maxGSSize = 500)
    } else {
      enrichment <- clusterProfiler::gseGO(geneList = protein_list,
                                         OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                         ont = "BP",
                                         pvalueCutoff = pval_lim,
                                         pAdjustMethod = "BH",
                                         minGSSize = 10,
                                         maxGSSize = 500)
    }
  } else if (database == "Reactome") {
    # Perform GSEA for Reactome
    enrichment <- ReactomePA::gsePathway(protein_list,
                                         organism = "human",
                                         pvalueCutoff = pval_lim,
                                         pAdjustMethod = "BH",
                                         verbose = FALSE)
  }

  if (!any(enrichment@result$p.adjust < pval_lim)) {
    message("No significant terms found.")
    return(NULL)
  }

  # Visualize the results
  dotplot <- clusterProfiler::dotplot(enrichment,
                                      showCategory = ncateg,
                                      font.size = fontsize,
                                      split=".sign") +
    ggplot2::facet_grid(.~.sign)

  ridgeplot <- clusterProfiler::ridgeplot(enrichment, showCategory = ncateg) +
    ggplot2::labs(x = "enrichment distribution") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = fontsize),
                   axis.text.y = ggplot2::element_text(size = fontsize),
                   text = ggplot2::element_text(size = fontsize))

  gseaplot <- clusterProfiler::gseaplot(enrichment,
                                        by = "all",
                                        title = enrichment$Description[1],
                                        geneSetID = 1)

  enrichment <- clusterProfiler::setReadable(enrichment, OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  cnetplot <- clusterProfiler::cnetplot(enrichment,
                                        showCategory = ncateg,
                                        categorySize = "pvalue",
                                        color.params = list(foldChange = protein_list),
                                        cex.params = list(category_label = (fontsize + 2)/12,
                                                          gene_label = (fontsize)/12))

  return(list("enrichment" = enrichment,
              "dotplot" = dotplot,
              "cnetplot" = cnetplot,
              "ridgeplot" = ridgeplot,
              "gseaplot" = gseaplot))
}
