utils::globalVariables(c("ENTREZID"))
#' Perform over-representation analysis (ORA) using clusterProfiler
#'
#' @param gene_list (character). A vector containing the gene names.
#' @param significance (character). The significance of the genes. It can be either "up", "down", or "all".
#' @param database (character). The database to perform the ORA. It can be either "KEGG" or "GO".
#'
#' @return enrichment (list). A list containing the results of the ORA.
#' @export
#'
#' @examples
#' #enrichment <- do_ora(c("TP53", "BRCA1", "BRCA2"), database = "KEGG")
do_ora <- function(gene_list, significance = c("up", "down", "all"), database = c("KEGG", "GO")) {
  database <- match.arg(database)
  significance <- match.arg(significance)

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

  return(enrichment)
}


#' Plot the results of the ORA
#'
#' @param enrichment (list). The results of the ORA.
#'
#' @return dotplot (plot). A dotplot showing the results of the ORA.
#' @export
#'
#' @examples
#' enrichment <- do_ora(c("TP53", "BRCA1", "BRCA2"), database = "KEGG")
#' plots <- plot_enrichment(enrichment)
plot_enrichment <- function(enrichment) {
  dotplot <- clusterProfiler::dotplot(enrichment)
}


#' Perform gene set enrichment analysis (GSEA) using clusterProfiler
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
