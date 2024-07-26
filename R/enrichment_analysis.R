utils::globalVariables(c("ENTREZID"))
#' Perform over-representation analysis
#'
#' `do_ora()` performs over-representation analysis (ORA) using the clusterProfiler package.
#' It also produces useful plots to visualize the results.
#'
#' @param gene_list A character vector containing the gene names.
#' @param database The database to perform the ORA. It can be either "KEGG" or "GO".
#'
#' @return A list containing the results of the ORA.
#' @export
#'
#' @examples
#' # Perform Differential Expression Analysis
#' de_res <- do_limma(example_data, example_metadata, wide = FALSE)
#'
#' # Extract the up-regulated proteins for AML
#' sig_up_proteins_aml <- de_res$de_results$AML |>
#'   dplyr::filter(sig == "significant up") |>
#'   dplyr::pull(Assay)
#'
#' # Perform ORA with GO database
#' do_ora(sig_up_proteins_aml, database = "GO")
do_ora <- function(gene_list, database = c("KEGG", "GO")) {
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
    # Perform GO enrichment analysis (Biological Process)
    enrichment <- clusterProfiler::enrichGO(gene = gene_list,
                                            OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                            ont = "BP")
  }

  if (!any(enrichment@result$p.adjust < 0.05)) {
    message("No significant terms found.")
    return(NULL)
  }

  # Visualize the results
  dotplot <- clusterProfiler::dotplot(enrichment)
  barplot <- barplot(enrichment, drop = TRUE,
                     showCategory = 10,
                     font.size = 8)
  goplot <- clusterProfiler::goplot(enrichment,
                                    showCategory = 10)
  cnetplot <- clusterProfiler::cnetplot(enrichment,
                                        categorySize = "pvalue",
                                        color.params = list(foldChange = gene_list))

  return(list("enrichment" = enrichment,
              "dotplot" = dotplot,
              "barplot" = barplot,
              "goplot" = goplot,
              "cnetplot" = cnetplot))
}


#' Perform gene set enrichment analysis
#'
#' This function performs gene set enrichment analysis (GSEA) using the clusterProfiler package.
#'
#' @param de_results (list). A list containing the results of the differential expression analysis.
#' @param database (character). The database to perform the GSEA. It can be either "KEGG" or "GO".
#'
#' @return enrichment_results (list). A list containing the results of the GSEA.
#' @export
#'
#' @examples
#' #enrichment_results <- do_gsea(de_results, database = "KEGG")
do_gsea <- function(de_results, database = c("KEGG", "GO")) {
  database <- match.arg(database)
  diseases <- de_results$de_results |> names()
  enrichment_results <- list()

  for (disease in diseases) {
    protein_list <- stats::setNames(de_results$de_results[[disease]]$logFC,
                                    de_results$de_results[[disease]]$Assay)
    sorted_proteins <- sort(protein_list, decreasing = TRUE)

    # From gene name to ENTREZID
    protein_conversion <- clusterProfiler::bitr(names(sorted_proteins),
                                                fromType = "SYMBOL",
                                                toType = "ENTREZID",
                                                OrgDb = org.Hs.eg.db::org.Hs.eg.db)

    protein_list <- stats::setNames(sorted_proteins, protein_conversion$ENTREZID)

    if (database == "KEGG") {
      # Perform GSEA for KEGG
      enrichment <- clusterProfiler::gseKEGG(geneList = protein_list,
                                             organism = "hsa",
                                             pvalueCutoff = 0.05,
                                             pAdjustMethod = "BH",
                                             minGSSize = 10,
                                             maxGSSize = 500)
    } else if (database == "GO") {
      # Perform GSEA for GO (Biological Process)
      enrichment <- clusterProfiler::gseGO(geneList = protein_list,
                                           OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                           ont = "BP",
                                           pvalueCutoff = 0.05,
                                           pAdjustMethod = "BH",
                                           minGSSize = 10,
                                           maxGSSize = 500)
    }

    enrichment_results[[disease]] <- enrichment
  }

  return(enrichment_results)
}
