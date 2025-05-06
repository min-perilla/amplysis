#' @title Heatmap Visualization
#'
#' @description
#' After performing heatmap analysis using the `heatmap()` function,
#' you will obtain a dataframe. This function allows you to visualize the
#' dataframe, thereby generating a high-quality heatmap.
#' If you want to customize the x-axis order of the heatmap, please input a
#' vector into the parameter cluster_cols. For example,
#' cluster_cols = c(4, 5, 6, 1, 2, 3), indicates sorting according to column
#' numbers 4 5 6 1 2 3.
#'
#' @param data (data.frame) Plotting Data
#' @param scale (character) scale is used to set normalization. "row" represents
#' row-wise normalization, "column" represents column-wise normalization, and
#' "none" represents no normalization.
#' @param cellwidth (numeric) Translation: Represents the width of a single
#' cell, default is "NA".
#' @param cellheight (numeric) Translation: Represents the height of a single
#' cell, default is "NA".
#'
#' @param color (character) The heatmap cell colors are generated automatically
#' with a gradient.
#' For example: c("#2196f3", "#a8d1f2", "#f4faff", "#ec9fa2", "#ec1c24")
#'
#' @param gaps_row (numeric vector) Used only when row clustering is not performed,
#' indicating the break positions in the row direction of the heatmap.
#' @param gaps_col (numeric vector) Used only when column clustering is not performed,
#' indicating the break positions in the column direction of the heatmap.
#'
#' @param custom_order (numeric or character vector) Custom order of column titles.
#' Only takes effect when column clustering is not enabled. When not `NULL`,
#' column clustering will be automatically disabled.
#' @param cluster_cols (logical) Whether to enable column clustering.
#' Custom ordering is not possible when clustering is enabled.
#'
#' @param clustering_distance_cols (character) Optional parameters for
#' `clustering_distance_rows` and `clustering_distance_cols`:
#'
#' 1 - "euclidean"   : Euclidean distance, the most commonly used distance metric.
#'
#' 2 - "correlation" : Correlation-based distance, computed using Pearson correlation
#' coefficients.
#'
#' 3 - "maximum"     : Chebyshev distance (maximum distance), calculated as the
#' maximum coordinate difference between points.
#'
#' 4 - "manhattan"   : Manhattan distance (city block distance), calculated as
#' the sum of absolute coordinate differences.
#'
#' 5 - "canberra"    : Canberra distance, based on the ratio of absolute coordinate
#' differences to their sum, suitable for sparse data.
#'
#' 6 - "binary"      : Binary distance, measuring differences in the binary
#' representation of data.
#'
#' 7 - "minkowski"   : Minkowski distance, a generalization of Euclidean distance.
#' @param cutree_cols (numeric) Number of clusters to divide the columns into
#' based on hierarchical clustering.
#'
#' @param cluster_rows (logical) Whether to enable row clustering.
#' @param clustering_distance_rows (character) Optional parameters for
#' `clustering_distance_rows` and `clustering_distance_cols`:
#'
#' 1 - "euclidean"   : Euclidean distance, the most commonly used distance metric.
#'
#' 2 - "correlation" : Correlation-based distance, computed using Pearson correlation
#' coefficients.
#'
#' 3 - "maximum"     : Chebyshev distance (maximum distance), calculated as the
#' maximum coordinate difference between points.
#'
#' 4 - "manhattan"   : Manhattan distance (city block distance), calculated as
#' the sum of absolute coordinate differences.
#'
#' 5 - "canberra"    : Canberra distance, based on the ratio of absolute coordinate
#' differences to their sum, suitable for sparse data.
#'
#' 6 - "binary"      : Binary distance, measuring differences in the binary
#' representation of data.
#'
#' 7 - "minkowski"   : Minkowski distance, a generalization of Euclidean distance.
#' @param cutree_rows (numeric) Number of clusters to divide the rows into
#' based on hierarchical clustering.
#'
#' @param clustering_method (character) Clustering method, options: 'ward.D',
#' 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'.
#'
#' @param treeheight_row (numeric) Row clustering tree height
#' @param treeheight_col (numeric) Col clustering tree height
#'
#' @param annotation_col (data.frame) Column annotation. When set to NULL or NA, the default column annotation generated based on group2 is used.
#' Custom input data can also be provided, for example:
#'
#' annotation_col <- data.frame(
#'
#'   Groups = c("A", "A", "A", "B", "B", "B", "D", "D", "D", "R", "R", "R"),
#'
#'   row.names = c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5",
#'   "Sample6", "Sample7", "Sample8", "Sample9", "Sample10", "Sample11", "Sample12")
#'
#' )
#'
#' @param annotation_row (data.frame) Row  annotation. Please provide a two-column data frame:
#' the first column should contain the row names of the heatmap, and the second
#' column should contain the annotation information.
#' @param annotation_colors (list or vector) Colors for annotation tracks. Can be provided as a list or a vector. Example:
#'
#' Vector format (column annotation): annotation_colors = c("A" = "purple", "B" = "orange", "R" = "pink", "D" = "yellow")
#'
#' List format (including both column and row annotations):
#'
#' annotation_colors = list(
#'
#'   Groups = c("A" = "purple", "B" = "orange", "R" = "pink", "D" = "yellow")  # Column annotation colors
#'
#'   Type = c("A" = "blue", "B" = "green", "C" = "red")  # Row annotation colors
#'
#' )
#'
#' @param title (character) The main title of the heatmap.
#' @param angle_col (numeric) Rotation angle of column labels. Options:
#' 0, 45, 90, 270, 315.
#'
#' @param fontsize (numeric) Base font size for the heatmap.
#' @param fontsize_row (numeric) Font size for row names.
#' @param fontsize_col (numeric) Font size for column names.
#' @param row_fontface_italic (logical) Whether to italicize row labels. Set to
#' TRUE for italics.
#'
#' @param legend_breaks (numeric vector) Defines the legend breakpoints, e.g.,
#' c(-1.2, 0, 1.2).
#' @param legend_labels (character vector) Labels for the legend breakpoints,
#' corresponding to legend_breaks, e.g., c("Low", "Medium", "High").
#'
#' @param display_numbers (logical) Whether to display numbers inside each cell.
#' @param number_format (numeric or character) Format for numbers inside cells.
#' Can be a numeric value (indicating the number of decimal places) or a format
#' string, such as "%.2f".
#' @param number_color (character) Color of the numbers inside cells.
#' @param fontsize_number (numeric) Font size for numbers inside cells.
#'
#' @param filename (character) File name for saving
#' @param file_width (numeric) Image width
#' @param file_height (numeric) Image height
#'
#' @return pheatmap
#' @export
#'
#' @examples
#' \dontrun{
#' heatmap1 = heatmap(otu = otu, tax = tax, metadata = metadata, id_col = 1,
#' group1 = "group", group2 = NULL, tax_cla = "genus", replicate_method = "mean",
#' row_n = 50)
#'
#' heatmap_plot(
#'   data = heatmap1, scale = "row", cellwidth = NA, cellheight = NA,
#'
#'   color = c("#2196f3", "#a8d1f2", "#f4faff", "#ec9fa2", "#ec1c24"),
#'
#'   gaps_row = NULL, gaps_col = NULL,
#'
#'   custom_order = NULL, cluster_cols = F, clustering_distance_cols = "euclidean",
#'   cutree_cols = NA,
#'
#'   cluster_rows = T, clustering_distance_rows = 1, cutree_rows = NA,
#'   clustering_method = "ward.D",
#'   treeheight_row = 100, treeheight_col = 50,
#'
#'   annotation_col = NULL,
#'   annotation_row = NA,
#'   annotation_colors = NULL,
#'
#'   title = "Heatmap",
#'   angle_col = 0,
#'
#'   fontsize = 18, fontsize_row = 16, fontsize_col = 18, row_fontface_italic = T,
#'   legend_breaks = NA, legend_labels = NA,
#'   display_numbers = F, number_format = 2, number_color = "grey30",
#'   fontsize_number = 0.8 * fontsize,
#'   filename = "heatmap", file_width = 16, file_height = 12)
#' }
#'
#' @importFrom stats dist hclust
#' @importFrom grDevices dev.cur
#'
# Please use tools:: showNonASCIIfile(file.R) to check for the presence of non ASCII characters.
# tools::showNonASCIIfile(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "heatmap_plot.R"))
heatmap_plot <- function(data, scale = "row", cellwidth = NA, cellheight = NA,
         color = c("#2196f3", "#a8d1f2", "#f4faff", "#ec9fa2", "#ec1c24"),
         gaps_row = NULL, gaps_col = NULL, custom_order = NULL, cluster_cols = F,
         clustering_distance_cols = "euclidean", cutree_cols = NA,
         cluster_rows = T, clustering_distance_rows = "euclidean",
         cutree_rows = NA, clustering_method = "ward.D", treeheight_row = 50,
         treeheight_col = 50, annotation_col = NA, annotation_row = NA,
         annotation_colors = NA, title = "Heatmap", angle_col = 0,
         fontsize = 18, fontsize_row = 16, fontsize_col = 18,
         row_fontface_italic = T, legend_breaks = NA, legend_labels = NA,
         display_numbers = F, number_format = "%.2f", number_color = "grey30",
         fontsize_number = 0.8 * fontsize, filename = "heatmap", file_width = 12,
         file_height = 12)
{
  # Close the graphics device to prevent unwanted output
  if (dev.cur() > 1) dev.off()

  # Plotting data
  mat = data[["data"]]

  ## Column annotation data
  flag_annotation_col = FALSE  # Flag to indicate whether custom column annotation is enabled

  # annotation_col = read_data("annotation_col.csv")

  # Check if annotation_col is valid (not NULL, empty string, whitespace-only string, or contains NA)
  if (is.null(annotation_col) ||
      (is.character(annotation_col) && trimws(annotation_col) == "") ||  # Check if it's an empty or whitespace-only string
      any(is.na(annotation_col)) ||
      (is.data.frame(annotation_col) && nrow(annotation_col) == 0)) {
    # Enable default annotation
    annotation_col = data[["annotation_col"]]
    cat("Using default column annotation", "\n", sep = "")
    if (is.null(annotation_col)) {
      annotation_col = NA
    }


    ## Custom column annotation
  } else {
    flag_annotation_col = TRUE  # Flag to indicate custom column annotation is enabled
    cat("Using custom column annotation", "\n", sep = "")

    # annotation_col = read_data("annotation_col.csv")

    # Process the annotation data
    annotation_col = as.data.frame(annotation_col)
    rownames(annotation_col) = annotation_col[, 1]       # Set the first column as row names
    annotation_col = annotation_col[, -1, drop = FALSE]  # Remove the first column
  }


  ##
  # Annotation track color processing (convert to a list format)

  # If custom column annotation is not enabled
  if (isFALSE(flag_annotation_col)) {
    col_name = colnames(data[["annotation_col"]])[1]  # Get the first column name

  } else {  # If custom column annotation is enabled
    col_name = colnames(annotation_col)[1]  # Get the first column name
  }


  # Check if annotation_colors is NULL or NA
  if (is.null(annotation_colors) || all(is.na(annotation_colors))) {
    annotation_colors = NULL
  } else {
    # Determine if annotation_colors is a vector or a list
    if (!is.list(annotation_colors)) {
      # If it's a vector, convert it to a list and set the name
      annotation_colors = setNames(list(annotation_colors), col_name)
    } else {
      # If it's already a list, update the name directly
      names(annotation_colors)[1] = col_name
    }
  }


  ## Color settings
  color = grDevices::colorRampPalette(color)(100) # Gradient colors
  # ----------------------------------------------------------------------------

  # custom_order = c("A", "B", "R", "D")
  # custom_order = c(1, 2, 4, 3)
  # custom_order = NULL
  # If custom_order is not NULL, apply manual sorting
  if (!is.null(custom_order)) {
    # Disable column clustering
    cluster_cols = FALSE
    cat("\033[31mColumn clustering: Disabled (Custom column ordering enabled)\033[0m\n")

    # If custom_order is a character vector (i.e., column names), convert it to column indices
    if (is.character(custom_order)) {
      if (!all(custom_order %in% colnames(mat))) {
        stop("Error: The specified custom_order contains non-existent column names. Please check your input!")
      }
      # Convert to column indices
      cat("Custom order: ", custom_order, "\n", sep = " ")
      custom_order <- match(custom_order, colnames(mat))
      cat("Custom order: ", custom_order, "(column number)\n", sep = " ")

      ## If already column indices, output corresponding column names for user reference
    } else {
      # Output column indices
      cat("Custom order: ", paste(custom_order, collapse = " "), "\n", sep = "")
      # Convert to corresponding column names
      custom_order_names <- colnames(mat)[custom_order]
      # Output column names
      cat("Custom order: ", paste(custom_order_names, collapse = " "), " (column name)\n", sep = "")
    }


    ##
    # Standardize by row, round to 2 decimal places, and transpose
    matrix2 <- round(t(apply(mat, MARGIN = 1, base::scale)), 2)

    # Restore column names for matrix2
    colnames(matrix2) <- colnames(mat)

    # Convert matrix2 to a data frame
    exprTable <- as.data.frame(t(matrix2))

    # Compute Euclidean distance
    row_dist <- stats::dist(exprTable, method = "euclidean")

    # Perform hierarchical clustering
    hclust_1 <- stats::hclust(row_dist)

    # Set custom column order
    hclust_1[[3]] <- custom_order

    # Assign to cluster_cols
    cluster_cols2 <- hclust_1


    ## Custom ordering not enabled
  } else {

    ## Disable column clustering
    if(isFALSE(cluster_cols)) {
      cluster_cols2 = cluster_cols
      cat("\033[31m", "Column clustering: Disabled (Manually disabled by user)", "\033[0m", "\n", sep = "")
    }

    cat("\033[31m", "Custom ordering: Disabled", "\033[0m", "\n", sep = "")
    cat("To enable custom column ordering for the heatmap, use the parameter \"custom_order\"\n",
        "e.g., custom_order = c(\"A\", \"B\", \"R\", \"D\")\n", sep = "")
  }


  ## Output the current column clustering status
  # Enable column clustering
  if (isTRUE(cluster_cols)) {
    cluster_cols2 = cluster_cols
    cat("\033[32m", "Column Clustering: Enabled", "\033[0m", "\n", sep = "")

    ## Check clustering method format
    valid_methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")

    if (!(clustering_method %in% valid_methods)) {
      stop("Invalid clustering method. Please choose from \n",
           "'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'.")
    } else {
      cat("Clustering method for columns: ", clustering_method, "\n", sep = "")
    }

    cat("To customize heatmap column ordering, use the parameter \"custom_order\"\n",
        "For example: custom_order = c(\"A\", \"B\", \"R\", \"D\")\n", sep = "")
  }
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # Distance metrics for row and column clustering: represented as a vector
  distance_methods <- c(
    "euclidean",   # 1: Euclidean distance
    "correlation", # 2: Pearson correlation distance
    "maximum",     # 3: Chebyshev distance
    "manhattan",   # 4: Manhattan distance
    "canberra",    # 5: Canberra distance
    "binary",      # 6: Binary distance
    "minkowski"    # 7: Minkowski distance
  )

  # Function to update the distance metric
  set_clustering_distance <- function(distance_param) {
    # If the input is a string, match it directly
    if (is.character(distance_param)) {
      if (distance_param %in% distance_methods) {
        return(distance_param)
      } else {
        stop("Error: Invalid distance metric provided. Please enter a valid string: 'euclidean', 'correlation', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'.")
      }
    }

    # If the input is numeric, convert it to the corresponding string
    if (is.numeric(distance_param)) {
      # Use modulo operation to ensure the number is within the valid range
      distance_method <- distance_methods[(distance_param - 1) %% length(distance_methods) + 1]
      return(distance_method)
    }

    stop("Error: Invalid input type. Please enter a string or a number.")
  }

  # Update distance metrics using numeric or string mapping
  clustering_distance_cols <- set_clustering_distance(clustering_distance_cols)
  clustering_distance_rows <- set_clustering_distance(clustering_distance_rows)

  # Output the selected distance metrics (including index and name)
  cat("Distance metric for column clustering: ", which(distance_methods == clustering_distance_cols), " - ", clustering_distance_cols, "\n", sep = "")
  cat("Distance metric for row clustering: ", which(distance_methods == clustering_distance_rows), " - ", clustering_distance_rows, "\n", sep = "")

  # Store the selected distance metrics
  selected_methods <- c(clustering_distance_cols, clustering_distance_rows)

  # Output the remaining available distance metrics with brief descriptions
  cat("\nOther available distance metrics:\n")
  for (i in 1:length(distance_methods)) {
    # Skip if the method has already been selected
    if (distance_methods[i] %in% selected_methods) {
      next
    }

    # Print the unselected distance metrics
    # Use sprintf for aligned output to ensure proper formatting
    cat(sprintf("%-2d- %-12s: %s\n", i, distance_methods[i],
                switch(i,
                       "1" = "Euclidean distance",
                       "2" = "Pearson correlation distance",
                       "3" = "Chebyshev distance, computes the maximum coordinate difference",
                       "4" = "Manhattan distance, sum of absolute coordinate differences",
                       "5" = "Canberra distance, suitable for measuring sparse data",
                       "6" = "Binary distance, computes differences in binary form",
                       "7" = "Minkowski distance, a generalization of Euclidean distance"
                )), sep = "")
  }
  cat("\n")
  # ----------------------------------------------------------------------------



  ###
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    message("The 'pheatmap' package is required but not installed.")
    message("Installing 'pheatmap' package...")
    install.packages("pheatmap")
  }


  ##
  # Preprocessing: Numbers inside the heatmap
  # Dynamically generate format based on user input
  # number_format = "%.3f"

  # Check input type and process
  if (is.numeric(number_format)) {
    number_format <- paste0("%.", number_format, "f")
  } else if (is.character(number_format) && grepl("^%\\.\\d+f$", number_format)) {
    number_format <- number_format
  } else {
    stop("Error: Invalid input format. Please enter an integer (number of decimal places) or a format string (e.g., '%.2f').")
  }

  # Preprocessing: Filename
  filename = paste0(filename, ".png")

  # Convert all row names to italics
  if(isTRUE(row_fontface_italic)) {
    labels_row <- parse(text = paste0("italic('", rownames(mat), "')"))
  } else {
    # Standard font
    labels_row = NULL
  }

  # Convert column names to bold
  labels_col <- parse(text = paste0("bold('", colnames(mat), "')"))

  ## Visualization
  # Close the graphics device to prevent no output
  if (dev.cur() > 1) dev.off()

  # Draw heatmap
  p <- pheatmap::pheatmap(
    mat = mat,
    scale = scale,

    # Cell styling
    border_color = "white",   # Set grid border color to white
    cellwidth = cellwidth,    # Width/height of a single cell, default is "NA"
    cellheight = cellheight,  # Width/height of a single cell, default is "NA"
    color = color,  # Color scheme
    # gaps_row = NULL,        # Used only when row clustering is disabled, specifies row gaps in the heatmap
    # gaps_col = c(1, 2, 3, 4),  # Used only when column clustering is disabled, specifies column gaps in the heatmap
    gaps_row = gaps_row,      # Used only when row clustering is disabled, specifies row gaps in the heatmap
    gaps_col = gaps_col,      # Used only when column clustering is disabled, specifies column gaps in the heatmap

    # Clustering settings
    # Columns
    cluster_cols = cluster_cols2,  # Column clustering
    clustering_distance_cols = clustering_distance_cols,  # Distance metric for column clustering
    cutree_cols = cutree_cols,        # Number of clusters to cut columns into

    # Rows
    cluster_rows = cluster_rows,   # Row clustering, enabled by default
    clustering_distance_rows = clustering_distance_rows,  # Distance metric for row clustering
    cutree_rows = cutree_rows,        # Number of clusters to cut rows into

    # Clustering method
    clustering_method = clustering_method,     # Clustering method options: 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'

    # Dendrogram settings
    treeheight_row = treeheight_row,  # Adjust row dendrogram height
    treeheight_col = treeheight_col,  # Adjust column dendrogram height

    ## Annotations
    # Column annotation (can be developed to use the second grouping in metadata)
    annotation_col = annotation_col,

    # Row annotation (reserved interface, can use a separate CSV file for annotations)
    annotation_row = annotation_row,

    # Annotation colors
    annotation_colors = annotation_colors,

    # Title
    main = title,                              # Heatmap title

    # Font size
    fontsize = fontsize,                       # Font size in the heatmap
    fontsize_row = fontsize_row,               # Row name font size
    fontsize_col = fontsize_col,               # Column name font size

    # Axis labels
    labels_row = labels_row,                   # Use row labels instead of row names
    labels_col = labels_col,                   # Use column labels instead of column names
    angle_col = angle_col,                     # Column name rotation angle

    # Legend
    # legend = F,                                # Remove legend
    # legend_breaks = c(-1.2, 0, 1.2),           # Set legend range
    # legend_labels = c("Low", "Medium", "High"),  # Labels for legend breakpoints
    legend_breaks = legend_breaks,             # Set legend range
    legend_labels = legend_labels,             # Labels for legend breakpoints

    # Numbers inside the heatmap
    display_numbers = display_numbers,  # Display numbers inside cells
    number_format = number_format,      # Number format inside cells
    number_color = number_color,        # Number color inside cells
    fontsize_number = fontsize_number,  # Font size for numbers inside cells

    # silent = TRUE,

    # Save file
    filename = filename, width = file_width, height = file_height
  )

  cat("\033[32mheatmap: success!\033[0m\n")
  cat("\033[0;32m", "The file \"", filename, "\" has been saved to \n",
      getwd(), "/", filename, "\033[0m\n", sep = "")

  # Close the graphics device to prevent no output
  if (dev.cur() > 1) dev.off()
  print(p)

  return(p)
}
