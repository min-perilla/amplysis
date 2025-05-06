#' @title Principal Component Analysis
#' @description
#' In bioinformatics analysis, Principal Component Analysis (PCA) is a commonly
#' used statistical technique for exploring and visualizing patterns and
#' structures within high-dimensional datasets. The "pca()" function utilizes
#' the built-in R function "stats::prcomp()" to analyze and calculate the
#' explanatory power of each principal component.
#'
#' @param otu otu table
#' @param metadata metadata table
#' @param id_col (integer) The OTU_ID column is in which column,
#' defaulting to 0 means there is no OTU_ID column, and the data is already
#' numeric.
#' @param group (Required, character) Grouping information. please enter the column name of
#' the grouping information in the metadata table.
#' @param replicate_method (character) Sample processing methods for the same group:
#' mean, sum, median, none.
#'
#' @return a list (consisting of two columns of data: plot data and the
#' explanatory power of each principal component)
#' @export
#'
#' @examples
#' \dontrun{pca(otu = otu, metadata = metadata, id_col = 1, group = "group",
#'              replicate_method = "none")}
#'
#' @importFrom stats prcomp var
#'
# Please use tools:: showNonASCIIfile(file.R) to check for the presence of non ASCII characters.
# tools::showNonASCIIfile(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "pca.R"))
pca <- function(otu, metadata, id_col = 1, group = "group", replicate_method = "none")
{
  # Check if the column name specified by 'group' exists in metadata
  if (!all(group %in% colnames(metadata))) {
    stop(paste("Some values in", ifelse(
      !all(group %in% colnames(metadata)), "group"),
      "are not present in metadata column names."))
  } else {
    cat("\033[32mgroup: `", group, "`\n\033[30m", sep = "")
  }

  # Format validation
  # Check if metadata contains "sample" and "replicate"
  if ("sample" %in% base::tolower(colnames(metadata)) &&
      "replicate" %in% base::tolower(colnames(metadata))) {
    cat("metadata --> DONE\n")
  } else {
    stop("Please ensure that the metadata table contains the `sample` column and the `replicate` column!",
         "\nsample: Sample ID (unique)",
         "\nreplicate: replicate sample identifier")
  }

  # Process metadata table
  metadata2 <- metadata[, c("sample", "replicate", group)]

  # Remove rows containing NA values
  na_rows <- apply(metadata2, 1, function(row) any(is.na(row)))
  if (any(na_rows)) {
    cat("The following row numbers contain NA values and have been discarded:\n")
    cat(which(na_rows), "\n")
    metadata2 <- metadata2[!na_rows, ]
  } else {
    cat("No 'NA' values found in the grouping information.\n")
  }

  # Process OTU table
  sample_values <- c(names(otu)[1], metadata2[["sample"]])  # Retrieve "sample" column values from metadata2
  keep_columns <- colnames(otu) %in% sample_values  # Create a logical vector indicating which columns to keep
  otu2 <- otu[, keep_columns]  # Retain only the selected columns in OTU

  # replicate sample processing
  allowedMethods <- base::tolower(c("mean", "sum", "median", "none"))
  replicate_method <- base::tolower(replicate_method)

  # Check replicate sample processing method
  if(!replicate_method %in% allowedMethods) {
    stop("Please enter the correct parameter for replicate_method:\n",
         "Processing is based on the `replicate` column in `metadata`. Samples with the same `replicate` value are considered replicates\n",
         "`mean`  : Take the average\n",
         "`sum`   : Sum\n",
         "`median`: Take the median\n",
         "`none`  : Do not process replicates\n")
  } else {
    cat("\033[32mreplicate replicate_method: `", replicate_method, "`\n\033[30m", sep = "")
  }

  # Process replicates
  if (replicate_method != "none") {
    otu3 <- otu2 %>%
      tidyr::gather(key = "sample", value = "abun", -1) %>%
      dplyr::left_join(metadata2, by = c("sample" = "sample"))
    cat("\033[32motu3 ---> DONE\n\033[30m")

    otu4 <- otu3 %>%
      dplyr::group_by_at(dplyr::vars(names(otu3)[1], dplyr::all_of("replicate"))) %>%
      dplyr::select(names(otu3)[1], dplyr::all_of("replicate"), dplyr::all_of(group), dplyr::all_of("abun")) %>%
      dplyr::summarise_if(is.numeric, ~round(match.fun(replicate_method)(.), 1)) %>%
      dplyr::ungroup()
    cat("\033[32motu4 ---> DONE\n\033[30m")

    # Convert to wide format
    otu5 <- otu4 %>%
      tidyr::spread(key = names(otu4)[1], value = "abun")
    otu5 <- data.frame(otu5)
    colnames(otu5)[1] <- "#OTU ID"
    otu5 <- as.data.frame(t(otu5))
    row_names <- row.names(otu5)
    otu5 <- data.frame(sample = row_names, otu5)
    row.names(otu5) <- NULL
    colnames(otu5) <- otu5[1, ]
    otu5 <- otu5[-1, ]

    # Restore original order
    otu6 <- otu5 %>%
      dplyr::arrange(match(otu5[[1]], otu2[[1]]))
    cat("\033[32motu6 ---> DONE\n\033[30m")

    # Synchronize metadata
    metadata3 = metadata2 %>%
      dplyr::select(dplyr::all_of("replicate"), dplyr::all_of(group)) %>%
      dplyr::distinct(replicate, .keep_all = TRUE)
    colnames(metadata3)[1] <- "sample"
  } else {
    otu6 = otu2
    metadata3 = metadata2
    cat("\033[32motu6 ---> DONE2\n\033[30m")
    cat("\033[32mmetadata3 ---> DONE2\n\033[30m")
  }

  # Obtain the column number of the OTU_ID
  if(id_col > 0) {
    cat("The column number for 'OTU ID Column' is: ", id_col, "\n", sep = "")
    otu6 <- as.data.frame(otu6)
    rownames(otu6) <- otu6[, id_col]
    otu6 <- otu6[, -id_col, drop = FALSE]
  }

  # Convert to numeric
  otu6[] <- lapply(otu6, function(x) as.numeric(trimws(x)))

  # Transpose OTU table
  otu_t <- t(otu6)
  otu_t <- otu_t[, apply(otu_t, 2, function(x) var(x) != 0)]

  # Perform PCA analysis
  PCA <- prcomp(otu_t, scale = TRUE)
  df_PCA_sum <- summary(PCA)
  PoA <- df_PCA_sum$importance[2, ] * 100

  PC12 <- as.data.frame(PCA$x[, 1:2])
  PC12 <- cbind(sample = rownames(PC12), PC12)
  PC12 <- merge(PC12, metadata3, by = "sample", sort = F)
  colnames(PC12)[colnames(PC12) == group] <- "group"

  # Return results
  result <- list(PCA = PC12,  # PCA data
                 PoA = PoA)   # Percentage of variance explained
  return(result)
}

