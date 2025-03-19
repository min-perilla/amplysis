#' @title Heatmap Analysis
#'
#' @description
#' To perform heatmap analysis, you will need three files:
#' 1. OTU table (containing abundance information);
#' 2. Taxonomy table (containing species annotation information);
#' 3. Metadata table (containing grouping information).
#' This function undergoes a series of processing steps,
#' such as clustering based on specified taxonomic levels using column numbers
#' or names from the tax table, to generate a dataframe directly usable for
#' pheatmap plotting. For quick plotting, you can use the built-in
#' `heatmap_plot()` function to generate high-quality plots efficiently.
#'
#' @param otu otu table
#' @param tax tax table
#' @param metadata metadata table containing grouping information.
#' @param id_col (integer) The column number of the OTU ID column in the
#' OTU table, by default, is 1.
#' @param tax_cla (character) Taxonomic level. Only column names from the tax
#' table can be entered, for example tax_cla = 'genus'.
#' @param group1 (character) Group 1, please enter the column name or column
#' number of the grouping information in the metadata table.
#' @param group2 (character) Group 2 for facetting plots, please enter the
#' column name or column number of the grouping information in the metadata
#' table.
#' @param parallel_method (character) Parallel sample processing method,
#' defaulting to mean. Options: mean (average), sum (summation),
#' median (median).
#' @param row_n (integer) Preserve the top N taxa (including the Nth) based on
#' abundance.
#'
#' @return A list: Heatmap plotting data and column annotation file
#' @export
#'
#' @examples
#' \dontrun{
#' heatmap(otu = otu, tax = tax, metadata = metadata, id_col = 1,
#' group1 = "group", group2 = NULL, tax_cla = "genus", parallel_method = "mean",
#' row_n = 50)}
#'
#' @importFrom dplyr across arrange desc group_by group_by_at left_join ungroup
#' where select slice sym summarise_if
#'
#' @importFrom tidyr all_of gather spread
#'
# Please use tools:: showNonASCIIfile(file.R) to check for the presence of non ASCII characters.
# tools::showNonASCIIfile(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "heatmap.R"))
heatmap <- function(otu, tax, metadata, id_col = 1, tax_cla = "genus",
                    group1 = "group", group2 = NULL, parallel_method = "mean", row_n = 35)
{
  # Process data
  ## Check arguments group1 and group2
  # Check if group1 and group2 are present in the column names of the metadata table
  if (!all(group1 %in% colnames(metadata)) || !all(group2 %in% colnames(metadata))) {
    stop(paste("Some values in", ifelse(!all(group1 %in% colnames(metadata)), "group1", "group2"), "are not present in metadata column names."))
  } else {
    cat("\033[32mgroup1: `", group1, "`\n",
        "group2: `", group2, "`\033[0m\n",
        sep = "")
  }

  ## Format check
  # Check if the metadata dataframe contains "sample" and "parallel"
  if ("sample" %in% base::tolower(colnames(metadata)) &&
      "parallel" %in% base::tolower(colnames(metadata))) {
    cat("metadata --> DONE\n")
  } else {
    stop("Please ensure that the metadata table contains `sample` and `parallel` columns!",
         "\nsample  : Sample ID (unique)",
         "\nparallel: Parallel sample identifier")
  }

  ## Process metadata table
  # Extract columns from metadata table: "sample", "parallel", argument group1, and argument group2
  metadata2 <- metadata[, c("sample", "parallel", group1, group2)]

  # Discard rows with NA values
  na_rows <- apply(metadata2, 1, function(row) any(is.na(row)))
  # Output and discard rows with NA values
  if (any(na_rows)) {
    cat("The following rows contain NA values and have been discarded:\n")
    cat(which(na_rows), "\n")
    metadata2 <- metadata2[!na_rows, ]
  } else {
    cat("No `NA` values found in the grouping information")
  }


  ## Process otu table
  # For otu, discard corresponding columns based on the grouping information
  sample_values <- c(names(otu)[1], metadata2[["sample"]])  # Get values from "sample" column of metadata2
  keep_columns <- colnames(otu) %in% sample_values  # Create a logical vector indicating which columns to keep
  otu2 <- otu[, keep_columns]    # Keep columns in otu where the column name is TRUE


  ##
  # Check parallel sample processing method
  if(parallel_method == "mean") {
    cat("In the `metadata` table, samples with the same `parallel` value are considered parallel samples\n")
    cat("Parallel sample processing method: mean\n")

  } else if(parallel_method == "sum") {
    cat("In the `metadata` table, samples with the same `parallel` value are considered parallel samples\n")
    cat("Parallel sample processing method: sum\n")

  } else if(parallel_method == "median") {
    cat("In the `metadata` table, samples with the same `parallel` value are considered parallel samples\n")
    cat("Parallel sample processing method: median\n")

  } else if(parallel_method == "none") {
    cat("No parallel sample processing\n")

  } else {
    stop("Please enter a valid parameter for the parallel_method argument:\n",
         "Process based on the `parallel` column in the `metadata` table. Samples with the same `parallel` value are considered parallel samples.\n",
         "`mean`  : Take the mean\n",
         "`sum`   : Sum the values\n",
         "`median`: Take the median\n",
         "`none`  : Do not process parallel samples\n")
  }


  ##
  # If the first column of the OTU table is OTU ID, it needs to be converted to row names,
  # so that the OTU table becomes a pure numeric matrix.
  if(id_col > 0) {

    # Convert to data.frame
    otu2 <- as.data.frame(otu2)
    tax <- as.data.frame(tax)


    ##
    # Extract OTU sample columns based on the values in the "sample" column of metadata2
    otu_colnames <- (metadata2$sample)
    matching_columns <- which(colnames(otu2) %in% otu_colnames)  # Match column names in the OTU table
    otu2 <- otu2[, c(id_col, matching_columns)]


    ##
    # Merge tax table and OTU table based on the `OTU ID` column (left join)
    otu2 <- merge(tax, otu2, by = id_col, all.x = T, sort = F)

    # Rename row names
    rownames(otu2) <- otu2[, id_col]

    # Remove OTU ID column
    otu2 <- otu2[, -id_col, drop = FALSE]

    cat("otu2 ---> DONE\n")


  } else {
    # No OTU_ID column

    # Convert to data.frame
    otu <- as.data.frame(otu)
    tax <- as.data.frame(tax)


    ##
    # Extract OTU sample columns based on the values in the "sample" column of metadata2
    otu_colnames <- (metadata2$sample)
    matching_columns <- which(colnames(otu) %in% otu_colnames)  # Match column names in the OTU table
    otu2 <- otu[, c(matching_columns)]

    # Merge tax table and OTU table based on row names (left join)
    otu2 <- merge(tax, otu2, by = "row.names", all.x = T, sort = F)

    # Convert the first column to row names and remove the first column
    rownames(otu2) <- otu2[, 1]

    # Remove OTU ID column
    otu2 <- otu2[, -id_col, drop = FALSE]

    cat("otu2 ---> DONE\n")
  }
  ###


  ##
  # Merge based on classification level (tax_cla)
  otu3 <- otu2 %>%
    dplyr::group_by(dplyr::select(otu2, tidyr::all_of(tax_cla))) %>%  # Group by classification level in tax
    dplyr::summarise_if(is.numeric, sum) %>%       # Sum the data for the same classification level
    dplyr::arrange(dplyr::desc(rowSums(dplyr::across(dplyr::where(is.numeric))))) %>%  # Sort by row sum in descending order
    dplyr::ungroup() %>%
    slice(1:row_n)  # Select the top n species by abundance, default is 50
  cat("otu3 ---> DONE\n")
  ###



  ## Convert to long format and left join with metadata2 table
  otu4 <- otu3 %>%
    # Convert the data frame from wide format to long format,
    # converting all columns except classification columns into two columns: sample and abun
    tidyr::gather(key = "sample", value = "abun", -1) %>%
    dplyr::left_join(metadata2, by = c("sample" = "sample"))  # Left join with metadata2 based on sample column
  cat("otu4 ---> DONE\n")
  ##


  ##
  # Process parallel samples
  if (parallel_method != "none") {
    otu5 <- otu4 %>%
      # Group by parallel and classification columns
      dplyr::group_by_at(dplyr::vars(all_of(tax_cla), "parallel")) %>%
      dplyr::select(all_of(tax_cla), "parallel", "abun") %>%
      dplyr::summarise_if(is.numeric, match.fun(parallel_method)) %>%
      dplyr::ungroup()

    # Left join again with metadata2 table
    metadata3 <- metadata2[, -which(names(metadata2) == "sample")]
    metadata3 <- unique(metadata3) # Remove duplicates
    otu5 <- merge(otu5, metadata3, by = "parallel", all.x = T, all.y = F, sort = F)

    cat("parallel_method --> DONE\n")
    cat("otu5 ---> DONE\n")
  } else {
    # Do not process parallel samples
    otu5 <- otu4
    cat("\033[31mparallelMethod --> NONE\033[0m\n")
    cat("otu5 ---> DONE\n")
  }


  # Merge based on groups
  otu6 <- otu5 %>%
    dplyr::group_by(dplyr::select(otu5, tidyr::all_of(tax_cla)),
                    dplyr::select(otu5, tidyr::all_of(group1))) %>%  # Group by classification level in tax
    dplyr::summarise_if(is.numeric, sum) %>%       # Sum the data for the same classification level
    dplyr::arrange(dplyr::desc(rowSums(dplyr::across(dplyr::where(is.numeric))))) %>%  # Sort by row sum in descending order
    dplyr::ungroup()
  cat("otu6 ---> DONE\n")


  # Convert data from long back to wide format
  otu7 <- otu6 %>%
    # dplyr::select(-dplyr::all_of("parallel")) %>%
    spread(key = {{ group1 }}, value = !!dplyr::sym("abun"))

  ##
  otu7 <- as.data.frame(otu7)         # Convert to data frame
  rownames(otu7) <- otu7[, 1]         # Set row names to the first column's value
  otu7 <- otu7[, -1]                  # Remove the original first column

  # Round all numeric columns to 3 decimal places
  otu7[] <- lapply(otu7, function(x) if(is.numeric(x)) round(x, 3) else x)
  cat("otu7 ---> DONE\n")
  ##


  ## Generate an annotation file based on the group2 parameter
  annotation_col = NULL  # Initialize annotation

  if (!is.null(group2) && nzchar(trimws(group2))) {
    # Extract the columns corresponding to group1 and group2 from metadata3
    annotation_col <- metadata3[, c(group1, group2), drop = FALSE]

    # Remove duplicate values in the column specified by group1
    annotation_col <- annotation_col[!duplicated(annotation_col[[group1]]), ]

    # Convert the first column into row names
    annotation_col = as.data.frame(annotation_col)
    rownames(annotation_col) = annotation_col[, 1]

    # Remove the first column
    annotation_col = annotation_col[, -1, drop = FALSE]

    # Rename the column to "Groups"
    colnames(annotation_col) = "Groups"
  }

  ## Print a message prompting the user to use the heatmap_plot() function for visualization
  cat("\033[32m--- Please use the `heatmap_plot()` function for visualization. ---\n\033[0m")

  result = NULL

  result = list(
    data = otu7,                     # Plotting data
    annotation_col = annotation_col  # Column annotation file
  )
  return(result)
}
