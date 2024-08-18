#' @title Visualization of Venn
#'
#' @description
#' After processing with the `venn()` function, you can use this function
#' for visualization.
#'
#' @param data Plotting data.
#' @param color_scheme (character) Color scheme.
#'
#' @param show_percentage (logical) Whether to display the percentage of each
#' set.
#' @param digits (numeric) Number of decimal places to retain for percentages.
#'
#' @param size_set_name (numeric) Title size of each dataset.
#' @param size_text (numeric) Font size of the text within the dataset.
#'
#' @param filename (character) File name for saving.
#' @param file_width (numeric) Width of the image.
#' @param file_height (numeric) Height of the image.
#'
#' @return ggvenn plot
#' @export
#'
#' @examples
#' \dontrun{
#' venn1 <- venn(otu, metadata, 1, "group")
#' venn_plot(venn1, color_scheme = c('#aaf200','#0082ff',"#d23aa4","#c777ff", "#79ff79"))
#' }
#'
#' @importFrom ggvenn ggvenn
#' @importFrom ggplot2 ggsave theme
#'
venn_plot <- function(data, color_scheme = NULL, show_percentage = T, digits = 2,
                      size_set_name = 12, size_text = 7, filename = NULL,
                      file_width = 12, file_height = 9
) {
  otu <- as.data.frame(data > 0)
  otu <- otu[rowSums(otu) > 0, ]


  ##
  if(is.null(color_scheme)) {
    color_scheme = c('#fff200','#0082ff',"#ff2d34","#7777ff", "#79ff79")
  }


  ##
  p1 <- ggvenn::ggvenn(
    otu,

    #
    fill_color = color_scheme,
    fill_alpha = 0.45,

    #
    stroke_color = "white",
    stroke_alpha = 1,
    stroke_size = 0.1,
    stroke_linetype = 1,

    #
    set_name_color = "black",
    set_name_size = size_set_name,

    #
    text_color = "black",
    text_size = size_text,

    #
    show_percentage = show_percentage,
    digits = digits,
  )


  ##
  p2 <- p1 +
    ggplot2::theme(plot.margin = ggplot2::margin(
      t = 35, r = 5, b = 15, l = 5, unit = "pt"))

  p2

  if(is.null(filename)) { filename = "Venn" }
  ggplot2::ggsave(filename = paste0(filename, ".png"), plot = p2,
                  width = file_width, height = file_height)
  ggplot2::ggsave(filename = paste0(filename, ".jpg"), plot = p2,
                  width = file_width, height = file_height)
  ggplot2::ggsave(filename = paste0(filename, ".pdf"), plot = p2,
                  width = file_width, height = file_height)


  ##
  cat("\033[32mtaxa_bar: success!\033[0m\n")
  cat("\033[0;32m", "The file has been saved to \n",
      getwd(), "\033[0m\n", sep = "")

  return(p2)
}
