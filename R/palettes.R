#' HPA color palettes
#'
#' This function returns a list of color palettes used by the Human Protein Atlas (HPA) project.
#'
#' @return A list of color palettes.
#' @keywords internal
get_hpa_palettes <- function() {
  palettes <- list(
    sex = c(
      "Male" = "blue",
      "Female" = "red"
    ),
    diff_exp = c(
      "not significant" = "grey",
      "significant down" = "blue",
      "significant up" = "red"
    ),
    cancers12 = c(
      "AML" = "#A6CEE3",
      "CLL" = "#2271B5",
      "LYMPH" = "#08585A",
      "MYEL" = "#66C2A5",
      "CRC" = "#B89B74",
      "LUNGC" = "#ADC74F",
      "GLIOM" = "#FFD321",
      "BRC" = "#E8A29A",
      "CVX" = "#9E0142",
      "ENDC" = "#B195AE",
      "OVC" = "#603479",
      "PRC" = "#E7662B"
    ),
    cancers15 = c(
      "AML" = "#A6CEE3",
      "CLL" = "#2271B5",
      "LYMPH" = "#08585A",
      "MYEL" = "#66C2A5",
      "CRC" = "#B89B74",
      "LUNGC" = "#ADC74F",
      "GLIOM" = "#FFD321",
      "BRC" = "#E8A29A",
      "CVX" = "#9E0142",
      "ENDC" = "#B195AE",
      "OVC" = "#603479",
      "PRC" = "#E7662B",
      "MENI" = "#FFFF80",
      "SI-NET" = "#504538",
      "PIT-NET" = "#FFFF00"
    )
  )

  return(palettes)
}


#' HPA color scales
#'
#' This function creates a ggplot2 scale for color aesthetics using the color palettes from the Human Protein Atlas (HPA) project.
#'
#' @param palette (character). The name of the palette to use. Available palettes are: "sex", "diff_exp", "cancers12", "cancers15".
#'
#' @return A ggplot2 scale for color aesthetics.
#' @export
#'
#' @examples
#' data <- data.frame(
#'   var1 = 1:10,
#'   var2 = seq(2, 20, by = 2),
#'   Sex = rep(c("Male", "Female"), each = 5)
#' )
#'
#' plot <- ggplot2::ggplot(data, ggplot2::aes(x = var1, y = var2, color = Sex)) +
#'   ggplot2::geom_point() +
#'   scale_color_hpa("sex")
scale_color_hpa <- function(palette) {
  hpa_palettes <- get_hpa_palettes()

  if (!palette %in% names(hpa_palettes)) {
    stop("Palette not found. Available palettes are: ", paste(names(hpa_palettes), collapse = ", "))
  }

  ggplot2::scale_color_manual(values = hpa_palettes[[palette]])
}


#' HPA fill scales
#'
#' @param palette (character). The name of the palette to use. Available palettes are: "sex", "diff_exp", "cancers12", "cancers15".
#'
#' @return A ggplot2 scale for fill aesthetics.
#' @export
#'
#' @examples
#' data <- data.frame(
#'   Sex = c("Male", "Female"),
#'   Count = c(60, 40)
#' )
#'
#' plot <- ggplot2::ggplot(data, ggplot2::aes(x = Sex, y = Count, fill = Sex)) +
#'   ggplot2::geom_bar(stat = "identity", position = "dodge") +
#'   scale_fill_hpa("sex")
scale_fill_hpa <- function(palette) {
  hpa_palettes <- get_hpa_palettes()

  if (!palette %in% names(hpa_palettes)) {
    stop("Palette not found. Available palettes are: ", paste(names(hpa_palettes), collapse = ", "))
  }

  ggplot2::scale_fill_manual(values = hpa_palettes[[palette]])
}
