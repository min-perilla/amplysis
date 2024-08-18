#' @title Upset Diagram Analysis
#'
#' @description
#' When the number of samples is greater than or equal to 5, the set diagram can
#' more effectively display the intersection information of OTUs among the
#' samples. Although the `Upset()` function allows the `metadata` parameter to
#' be NULL, we still recommend using the `group` column in the metadata file to
#' control the grouping information.
#'
#'
#' @param otu otu table
#' @param metadata metadata table. The metadata parameter can be `NULL`. When it
#' is NULL, the function will perform the analysis based on the OTU table.
#' When it is not `NULL`, the function will preprocess the data according to the
#' group information in the metadata.
#' @param id_col (integer) The column number of the OTU ID column in the
#' OTU table, by default, is 1.
#' @param group (character) Group 1, please enter the column name of
#' the grouping information in the metadata table.
#' @param parallel_method (character) Sample processing methods for the same group:
#' average, sum, median, none.
#'
#' @return otu table
#' @export
#'
#' @examples
#' \dontrun{
#' Upset(otu, metadata, 1, "group", "mean")
#' }
#'
#' @importFrom dplyr group_by group_cols summarize_at
#' @importFrom ggplot2 vars
#'
Upset <- function(otu, metadata = NULL, id_col = 1,
         group = "group", parallel_method = "mean")
{
  ##
  # Obtain the column number of the OTU_ID.
  if(id_col > 0) {
    cat("The column number for 'OTU ID Column' is: ", id_col, "\n", sep = "")
    otu <- as.data.frame(otu)            # Convert to data.frame
    rownames(otu) <- otu[, id_col]       # Rename row names
    otu <- otu[, -id_col, drop = FALSE]  # Remove the OTU ID column

  } else {}  # No OTU_ID column


  ##
  if(!is.null(metadata)) {

    if (!all(group %in% colnames(metadata))) {
      stop(paste("Some values in", ifelse(
        !all(group %in% colnames(metadata)), "group"),
        "are not present in metadata column names."))
    } else {
      cat("\033[32mgroup: `", group, "`\n\033[30m", sep = "")
    }


    ##
    if ("sample" %in% base::tolower(colnames(metadata)) &&
        "parallel" %in% base::tolower(colnames(metadata))) {
      cat("\033[32mmetadata --> DONE\n\033[30m")
    } else {
      stop("Please ensure that the metadata table contains the `sample` column
           and the `parallel` column!",
           "\nsample: Sample ID (unique)",
           "\nparallel: Parallel sample identifier")
    }


    ##
    metadata2 <- metadata[, c("sample", "parallel", group)]
    na_rows <- apply(metadata2, 1, function(row) any(is.na(row)))
    if (any(na_rows)) {
      cat("metadata: The following row numbers contain NA values and have been
          discarded:\n")
      cat(which(na_rows), "\n")
      metadata2 <- metadata2[!na_rows, ]

    } else {
      cat("metadata: No 'NA' values found in the grouping information.")
    }

    #
    sample_values <- c(names(otu)[1], metadata2[["sample"]])
    keep_columns <- colnames(otu) %in% sample_values
    otu2 <- otu[, keep_columns]


    ##
    allowedMethods <- base::tolower(c("mean", "sum", "median", "none"))
    parallel_method <- base::tolower(parallel_method)
    if(!parallel_method %in% allowedMethods) {
      stop("Please enter the correct argument for the parameter 'parallel_method':\n",
           "Process according to the 'parallel' column in the `metadata` table,
           samples with the same 'parallel' value are considered parallel
           samples\n",
           "`mean`  : Calculate the average\n",
           "`sum`   : Calculate the sum\n",
           "`median`: Calculate the median\n",
           "`none`  : Do not process parallel samples\n")
    } else {
      cat("\033[32mParallel parallel_method: `", parallel_method, "`\n\033[30m", sep = "")
    }


    ###
    otu_t <- as.data.frame(t(otu2))
    otu_t_g <- merge(otu_t, metadata2, by.x = "row.names", by.y = "sample",
                     all.x = TRUE, sort = F)
    row.names(otu_t_g) <- otu_t_g[, 1]
    otu_t_g <- otu_t_g[, -1]
    metadata3 <- otu_t_g[-c(1:ncol(otu_t))]


    ##
    otu_t2 <- otu_t %>%
      dplyr::group_by(metadata3[[group]]) %>%
      dplyr::summarize_at(ggplot2::vars(-dplyr::group_cols()), parallel_method)


    ##
    otu_t2 <- as.data.frame(otu_t2)
    row.names(otu_t2) <- otu_t2[, 1]
    otu_t2 <- otu_t2[, -1]
    otu3 <- as.data.frame(t(otu_t2))


    ###
  } else {
    #
    otu3 <- otu
  }


  ##
  return(otu3)
}
