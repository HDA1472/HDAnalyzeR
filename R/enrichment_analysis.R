utils::globalVariables(c("ENTREZID"))
#' Perform over-representation analysis (ORA) using clusterProfiler
#'
#' @param gene_list (vector). A list of gene names.
#' @param database (character). The database to perform the ORA. It can be either "KEGG" or "GO".
#'
#' @return enrichment (list). A list containing the results of the ORA.
#' @export
#'
#' @examples
#' enrichment <- do_ora(c("TP53", "BRCA1", "BRCA2"), database = "KEGG")
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

# kegg_result <- do_ora(gene_list, database = "KEGG")
# clusterProfiler::dotplot(kegg_result)
