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
# @importFrom circlize circos.axis circos.clear circos.text circos.track
# circos.par get.cell.meta.data
chord_plot = function(data, color_scheme = NULL, size_axis = 0.8,
         size_label = 1, label_height = 0.35, filename = "chord",
         file_width = 16, file_height = 10)
{

  {
    # Extract unique values from the "from" column
    col_name = unique(data[["from"]])

    # Custom colors
    if(!is.null(color_scheme)) {
      # Apply custom colors
      color = NULL                       # Initialize
      color[col_name] = color_scheme[1: length(col_name)]  # Assign colors
      grid_col = color
    } else {
      grid_col = NULL
    }
  }

  # ?circlize::chordDiagram

  # Set output filename
  filename = paste0(filename, ".pdf")
  pdf(file = filename, width = file_width, height = file_height, pointsize = 16)

  ##
  # Chord Diagram

  circlize::circos.par(track.height = 0.1)

  p1 = circlize::chordDiagram(
    x = data,             # Plotting data
    grid.col = grid_col,  # Track colors

    # Add annotation tracks, options: "grid" (labels, annotations, or other data),
    #                                "name" (sector name display)
    #                                "axis" (sector axis display)
    annotationTrack = c("grid"),

    # Set track height
    annotationTrackHeight = 0.05,  # Adjust track width
    preAllocateTracks = list(track.height = max(graphics::strwidth(unlist(dimnames(data)))))
  )

  ##
  # Add circos.axis to adjust axis label size
  circlize::circos.track(
    track.index = 2,
    panel.fun = function(x, y) {
      circlize::circos.axis(
        h = "top",               # Axis position
        labels.cex = size_axis   # Axis label size (adjust as needed)
      )
    },
    bg.border = NA
  )

  # Create or modify the first track
  circlize::circos.track(
    track.index = 1,  # Track index 1

    # Define a custom drawing function with x and y coordinates
    panel.fun = function(x, y) {
      xlim = circlize::get.cell.meta.data("xlim")  # Get x-axis limits for the current sector
      ylim = circlize::get.cell.meta.data("ylim")  # Get y-axis limits for the current sector
      sector.name = circlize::get.cell.meta.data("sector.index")  # Get the sector name (index)

      # Draw text within the current sector
      circlize::circos.text(
        mean(xlim),              # X position: midpoint of x-axis range
        ylim[1] + label_height,  # Y position: lower boundary of y-axis range
        sector.name,             # Text content (sector name)
        facing = "clockwise",    # Text orientation: clockwise
        niceFacing = TRUE,       # Auto-adjust text direction for better readability
        adj = c(0, 0.5),         # Left-aligned, vertically centered
        col = "black",           # Text color: black
        cex = size_label         # Text size
      )}, bg.border = NA         # Set track background border to NA (no border)
  )

  circlize::circos.clear()

  # Close the graphics device to ensure output is saved
  if (dev.cur() > 1) dev.off()


  ##
  # Return chord diagram
  return(p1)
}

