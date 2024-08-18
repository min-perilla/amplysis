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
#' `heat_map_plot()` function to generate high-quality plots efficiently.
#'
#' @param otu otu table
#' @param tax tax table
#' @param metadata metadata table containing grouping information.
#' @param id_col (integer) The column number of the OTU ID column in the
#' OTU table, by default, is 1.
#' @param tax_cla (character) Taxonomic level. Only column names from the tax
#' table can be entered, for example tax_cla = 'genus'.
#' @param group (character) Group 1, please enter the column name or column
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
#' @return Heatmap analysis dataframe
#' @export
#'
#' @examples
#' \dontrun{heat_map(otu = otu, tax = tax, metadata = metadata, id_col = 1,
#' group = "group", group2 = NULL, tax_cla = "genus", parallel_method = "mean", row_n = 50)}
#'
#' @importFrom dplyr across arrange desc group_by group_by_at left_join ungroup
#' where select slice sym summarise_if
#'
#' @importFrom tidyr all_of gather spread
#'
# Please use tools:: showNonASCIIfile(file.R) to check for the presence of non ASCII characters.
# tools::showNonASCIIfile(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "heat_map.R"))
heat_map <- function(otu, tax, metadata, id_col = 1, tax_cla = "genus",
         group = "group", group2 = NULL, parallel_method = "mean", row_n = 35)
{
  if (!all(group %in% colnames(metadata)) ||
      !all(group2 %in% colnames(metadata))) {
    stop(paste("Some values in", ifelse(
      !all(group %in% colnames(metadata)), "group1", "group2"),
      "are not present in metadata column names."))
  } else {
    cat("\033[32mgroup1: `", group, "`\n",
        "group2: `", group2, "`\033[0m\n",
        sep = "")
  }

  if ("sample" %in% base::tolower(colnames(metadata)) &&
      "parallel" %in% base::tolower(colnames(metadata))) {
    cat("metadata --> DONE\n")
  } else {
    stop("Please ensure that the metadata table contains the `sample` column and the `parallel` column!",
         "\nsample: Sample ID (unique)",
         "\nparallel: Parallel sample identifier")
  }

  metadata2 <- metadata[, c("sample", "parallel", group, group2)]

  na_rows <- apply(metadata2, 1, function(row) any(is.na(row)))
  if (any(na_rows)) {
    cat("The following row numbers contain NA values and have been discarded:\n")
    cat(which(na_rows), "\n")
    metadata2 <- metadata2[!na_rows, ]
  } else {
    cat("No 'NA' values found in the grouping information.")
  }


  ##
  if(parallel_method == "mean") {
    cat("In the `metadata` table, samples with the same `parallel` value are considered parallel samples.")
    cat("Parallel sample treatment method: mean")

  } else if(parallel_method == "sum") {
    cat("In the `metadata` table, samples with the same `parallel` value are considered parallel samples.")
    cat("Parallel sample treatment method: sum")

  } else if(parallel_method == "median") {
    cat("In the `metadata` table, samples with the same `parallel` value are considered parallel samples.")
    cat("Parallel sample treatment method: median")

  } else if(parallel_method == "none") {
    cat("No processing for parallel samples.")

  } else {
    stop("Please enter the correct parameter for the parallel_method:\n",
         "Process according to the `parallel` column in the `metadata` table, where samples with the same `parallel` value are considered parallel samples.\n",
         "`mean`: Take the mean\n",
         "`sum`: Summation\n",
         "`median`: Take the median\n",
         "`none`: No processing for parallel samples\n")
  }


  ##
  if(id_col > 0) {

    otu <- as.data.frame(otu)
    tax <- as.data.frame(tax)


    ##
    otu_colnames <- (metadata2$sample)
    matching_columns <- which(colnames(otu) %in% otu_colnames)
    otu2 <- otu[, c(id_col, matching_columns)]


    ##
    otu2 <- merge(tax, otu2, by = id_col, all.x = T, sort = F)
    rownames(otu2) <- otu2[, id_col]
    otu2 <- otu2[, -id_col, drop = FALSE]
    cat("otu2 ---> DONE\n")


  } else {

    otu <- as.data.frame(otu)
    tax <- as.data.frame(tax)


    ##
    otu_colnames <- (metadata2$sample)
    matching_columns <- which(colnames(otu) %in% otu_colnames)
    otu2 <- otu[, c(matching_columns)]
    otu2 <- merge(tax, otu2, by = "row.names", all.x = T, sort = F)
    rownames(otu2) <- otu2[, 1]
    otu2 <- otu2[, -id_col, drop = FALSE]
    cat("otu2 ---> DONE\n")
  }
  ###


  ##
  otu3 <- otu2 %>%
    dplyr::group_by(dplyr::select(otu2, tidyr::all_of(tax_cla))) %>%
    dplyr::summarise_if(is.numeric, sum) %>%
    dplyr::arrange(dplyr::desc(rowSums(dplyr::across(dplyr::where(is.numeric))))) %>%
    dplyr::ungroup() %>%
    dplyr::slice(1:row_n)
  cat("otu3 ---> DONE\n")
  ###


  ##
  otu4 <- otu3 %>%
    tidyr::gather(key = "sample", value = "abun", -1) %>%
    dplyr::left_join(metadata2, by = c("sample" = "sample"))
  cat("otu4 ---> DONE\n")
  ##


  ##
  if (parallel_method != "none") {
    otu5 <- otu4 %>%
      dplyr::group_by_at(dplyr::vars(all_of(tax_cla), "parallel")) %>%
      dplyr::select(all_of(tax_cla), "parallel", "abun") %>%
      dplyr::summarise_if(is.numeric, match.fun(parallel_method)) %>%
      dplyr::ungroup()

    metadata3 <- metadata2[, -which(names(metadata2) == "sample")]
    metadata3 <- unique(metadata3)
    otu5 <- merge(otu5, metadata3, by = "parallel", all.x = T, all.y = F, sort = F)

    cat("parallel_method --> DONE\n")
    cat("otu5 ---> DONE\n")
  } else {
    otu5 <- otu5
    cat("\033[31mparallelMethod --> NONE\033[0m\n")
    cat("otu5 ---> DONE\n")
  }


  otu6 <- otu5 %>%
    dplyr::group_by(dplyr::select(otu5, tidyr::all_of(tax_cla)),
                    dplyr::select(otu5, tidyr::all_of(group))) %>%
    dplyr::summarise_if(is.numeric, sum) %>%
    dplyr::arrange(dplyr::desc(rowSums(dplyr::across(dplyr::where(is.numeric))))) %>%
    dplyr::ungroup()
  cat("otu6 ---> DONE\n")



  otu7 <- otu6 %>%
    spread(key = {{ group }}, value = !!dplyr::sym("abun"))

  ##
  otu7 <- as.data.frame(otu7)
  rownames(otu7) <- otu7[, 1]
  otu7 <- otu7[, -1]
  cat("otu7 ---> DONE\n")
  ##

  cat("\033[32m--- Please use the `heat_map_plot()` function for visualization. ---\n\033[0m")
  return(otu7)
}
