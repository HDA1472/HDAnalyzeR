utils::globalVariables(c("pmid", "year", "journal", "firstname", "lastname", "title", "First_author"))
#' Automated PubMed literature search
#'
#' This function searches for articles on proteins and diseases in PubMed.
#' A list of proteins and diseases is provided as input. The function retrieves the articles for each protein-disease pair.
#' The input should be in the correct format: a list with diseases as names and protein vectors associated with each disease as elements.
#' The disease and gene names should be correct in order for the query to be successful.
#'
#' @param prot_dis_list (list). A list of proteins and diseases. The names of the list are the diseases and the elements are vectors of proteins.
#' @param max_articles (numeric). The maximum number of articles to retrieve for each protein-disease pair. Default is 10.
#'
#' @return articles_database (list). A list of dataframes. Each dataframe contains the articles found for a protein-disease pair.
#' @export
#'
#' @examples
#' prot_dis_list <- list("acute myeloid leukemia" = c("FLT3", "EPO"),
#'                       "chronic lymphocytic leukemia" = c("PARP1", "FCER2"))
#' articles_database <- literature_search(prot_dis_list, max_articles = 2)
literature_search <- function(prot_dis_list, max_articles = 10) {

  diseases <- names(prot_dis_list)
  proteins <- prot_dis_list

  articles_database <- list()
  for (disease in diseases) {
    articles_database[[disease]] <- list()
    for (protein in proteins[[disease]]) {
      # Collect articles from PubMed
      ids <- easyPubMed::get_pubmed_ids(paste0(protein, '[All Fields] AND "', disease, '[All Fields]"'))
      message(paste0("Searching for articles on ", protein, " and ", disease))
      abstracts_xml <- easyPubMed::fetch_pubmed_data(pubmed_id_list = ids, retmax = max_articles)

      # Extract data and store them in dataframe
      if(!is.null(abstracts_xml)) {
        abstracts_list <- easyPubMed::articles_to_list(abstracts_xml)

        abstract_table <- abstracts_list |>
          purrr::map_df(function(abstract) {
            easyPubMed::article_to_df(pubmedArticle = abstract, autofill = FALSE) |>
              utils::head(1) |>
              dplyr::select(pmid, year, journal, firstname, lastname, title) |>
              dplyr::mutate(First_author = paste0(lastname, ", ", firstname)) |>
              dplyr::relocate(First_author, .before = year) |>
              dplyr::select(-firstname, -lastname)
          })
      }
      articles_database[[disease]][[protein]] <- abstract_table
    }
  }
  return(articles_database)
}

