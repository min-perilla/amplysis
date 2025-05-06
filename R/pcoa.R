#' @title Principal Co-ordinates Analysis
#' @description
#' Principal Coordinates Analysis(PCoA) is a visualization method for studying
#' data similarities or dissimilarities, enabling the observation of differences
#' between individuals or groups. The "pcoa()" function calculates the
#' Bray-Curtis distance using "vegan::vegdist()", then performs PCoA using the
#' built-in R function "stats::cmdscale()", etaining the eigenvalues.
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
#' \dontrun{pcoa(otu = otu, metadata = metadata, id_col = 1, group = "group",
#'               replicate_method = "none")}
#'
#' @importFrom vegan vegdist
#' @importFrom stats cmdscale
#'
pcoa <- function(otu, metadata, id_col = 1, group = "group", replicate_method = "none")
{
  # Check if the specified column exists in metadata
  if (!all(group %in% colnames(metadata))) {
    stop(paste("Some values in", ifelse(
      !all(group %in% colnames(metadata)), "group"),
      "are not present in metadata column names."))
  } else {
    cat("\033[32mgroup: `", group, "`\n", sep = "")
  }

  ## Format validation
  # Check if metadata contains "sample" and "replicate"
  if ("sample" %in% base::tolower(colnames(metadata)) &&
      "replicate" %in% base::tolower(colnames(metadata))) {
    cat("metadata --> DONE\n")
  } else {
    stop("Please ensure that the metadata table contains the `sample` and `replicate` columns!",
         "\nsample: Sample ID (unique)",
         "\nreplicate: replicate sample identifier")
  }

  ## Process metadata
  # Extract relevant columns from metadata
  metadata2 <- metadata[, c("sample", "replicate", group)]

  # Remove rows with NA values
  na_rows <- apply(metadata2, 1, function(row) any(is.na(row)))
  if (any(na_rows)) {
    cat("The following row numbers contain NA values and have been discarded:\n")
    cat(which(na_rows), "\n")
    metadata2 <- metadata2[!na_rows, ]
  } else {
    cat("No 'NA' values found in the grouping information.\n")
  }

  ## Process OTU table
  # Retain only relevant columns based on sample IDs in metadata
  sample_values <- c(names(otu)[1], metadata2[["sample"]])
  keep_columns <- colnames(otu) %in% sample_values
  otu2 <- otu[, keep_columns]

  ## Handle technical replicates
  allowedMethods <- base::tolower(c("mean", "sum", "median", "none"))
  replicate_method <- base::tolower(replicate_method)

  if(!replicate_method %in% allowedMethods) {
    stop("Please enter a valid replicate_method:\n",
         "Based on the `replicate` column in the metadata table, samples with the same `replicate` value are considered technical replicates\n",
         "`mean`  : Calculate the mean\n",
         "`sum`   : Calculate the sum\n",
         "`median`: Calculate the median\n",
         "`none`  : No processing of technical replicates\n")
  } else {
    cat("\033[32mreplicate replicate_method: `", replicate_method, "`\n\033[30m", sep = "")
  }

  if (replicate_method != "none") {
    ## Convert to long format and merge with metadata
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

    ## Convert back to wide format
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

    ## Restore original order
    otu6 <- otu5 %>%
      dplyr::arrange(match(otu5[[1]], otu2[[1]]))
    cat("\033[32motu6 ---> DONE\n\033[30m")

    ## Sync metadata
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

  ## Handle OTU ID column
  if(id_col > 0) {
    cat("The column number for 'OTU ID Column' is: ", id_col, "\n", sep = "")
    otu6 <- as.data.frame(otu6)
    rownames(otu6) <- otu6[, id_col]
    otu6 <- otu6[, -id_col, drop = FALSE]
    cat("\033[32mid_col ---> DONE\n\033[30m")
  }

  # Convert all values to numeric
  otu6[] <- lapply(otu6, function(x) as.numeric(trimws(x)))

  ## Perform PCoA
  otu_t <- t(otu6)
  otu_t <- otu_t[, apply(otu_t, 2, function(x) var(x) != 0)]
  otu.distance <- vegan::vegdist(otu_t)
  pc <- cmdscale(otu.distance, eig = TRUE)

  PoA <- round(pc$eig / sum(pc$eig) * 100, digits = 2)
  PC12 <- as.data.frame(pc$points[, 1:2])
  colnames(PC12) <- c("PC1", "PC2")
  PC12 <- cbind(sample = rownames(PC12), PC12)
  PC12 <- merge(PC12, metadata3, by = "sample", sort = F)
  colnames(PC12)[colnames(PC12) == group] <- "group"

  result <- list(PCoA = PC12, PoA = PoA)
  return(result)
}

