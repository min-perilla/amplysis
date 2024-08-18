#' @title Chord Diagram Visualization
#' @description
#' After processing with the `chord()` function, you can use this function
#' for visualization.
#'
#' @param data Plotting data.
#' @param color_scheme (character) Color scheme.
#'
#' @param size_axis (numeric) Font size on the sector axis.
#' @param size_label (numeric) Font size of sector labels.
#' @param label_height (numeric) Track height of label text.
#'
#' @param filename (character) File name for saving.
#' @param file_width (numeric) Width of the image.
#' @param file_height (numeric) Height of the image.
#'
#' @return Chord Diagram
#' @export
#'
#' @examples
#' \dontrun{
#' chord_plot(
#' data = chord_data,
#' color_scheme = c("#27e6ff", "#42ff0e", "#33BEB7", "#F66320", "#FBA127",
#'                  "#A463D7", "#DB3937", "#ffaec8", "#828282"))
#' }
#'
#' @importFrom graphics par plot.new legend
#@importFrom circlize circos.axis circos.clear circos.text circos.track
#circos.par get.cell.meta.data
chord_plot = function(data, color_scheme = NULL, size_axis = 0.8,
         size_label = 1, label_height = 0.35, filename = "chord",
         file_width = 16, file_height = 10)
{

  ###
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("The 'circlize' package is required but not installed. Please install
         it using 'install.packages(\"circlize\")'.")
  }


  ##
  col_name = unique(data[["from"]])
  if(!is.null(color_scheme)) {
    color = NULL
    color[col_name] = color_scheme[1: length(col_name)]
    grid_col = color
  } else {
    grid_col = NULL
  }

  filename = paste0(filename, ".pdf")
  pdf(file = filename, width = file_width, height = file_height, pointsize = 16)
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(5, 2))
  graphics::par(mar = c(1, 1, 1, 1))
  circlize::circos.par(track.height = 0.1)

  p1 = circlize::chordDiagram(
    x = data,
    grid.col = grid_col,
    annotationTrack = c("grid"),
    annotationTrackHeight = 0.05,
    preAllocateTracks = list(track.height = max(graphics::strwidth(unlist(dimnames(data)))))
  )


  ##
  circlize::circos.track(
    track.index = 2,
    panel.fun = function(x, y) {
      circlize::circos.axis(
        h = "top",
        labels.cex = size_axis
      )
    },
    bg.border = NA
  )

  #
  circlize::circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      xlim = circlize::get.cell.meta.data("xlim")
      ylim = circlize::get.cell.meta.data("ylim")
      sector.name = circlize::get.cell.meta.data("sector.index")

      circlize::circos.text(
        mean(xlim),
        ylim[1] + label_height,
        sector.name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        col = "black",
        cex = size_label
      )}, bg.border = NA
  )

  circlize::circos.clear()

  graphics::par(mar = c(1, 1, 1, 1))
  graphics::plot.new()


  ##
  graphics::legend(
    "left",
    pch = 20,
    legend = col_name,
    col = color,
    bty = "n",
    cex = 1.3,
    pt.cex = 3,
    border = "black",
    y.intersp = 1,
    adj = 0
  )


  ##
  grDevices::dev.off()
  graphics::plot.new()


  ##
  return(p1)
}
