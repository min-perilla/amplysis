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
#'
#' @return a list (consisting of two columns of data: plot data and the
#' explanatory power of each principal component)
#' @export
#'
#' @examples
#' \dontrun{pca(otu, 1, metadata)}
#'
#' @importFrom stats prcomp
#'
# Please use tools:: showNonASCIIfile(file.R) to check for the presence of non ASCII characters.
# tools::showNonASCIIfile(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "pca.R"))
pca <- function(otu, metadata, id_col = 1, group = "group")
{
  #
  if (!all(group %in% colnames(metadata))) {
    stop(paste("Some values in", ifelse(
      !all(group %in% colnames(metadata)), "group"),
      "are not present in metadata column names."))
  } else {
    cat("\033[32mgroup: `", group, "`\n", sep = "")
  }


  ##
  #
  if ("sample" %in% base::tolower(colnames(metadata)) &&
      "parallel" %in% base::tolower(colnames(metadata))) {
    cat("metadata --> DONE\n")
  } else {
    stop("Please ensure that the metadata table contains the `sample` column and
         the `parallel` column!",
         "\nsample: Sample ID (unique)",
         "\nparallel: Parallel sample identifier")
  }


  ##
  #
  metadata2 <- metadata[, c("sample", "parallel", group)]

  #
  na_rows <- apply(metadata2, 1, function(row) any(is.na(row)))
  #
  if (any(na_rows)) {
    cat("The following row numbers contain NA values and have been discarded:\n")
    cat(which(na_rows), "\n")
    metadata2 <- metadata2[!na_rows, ]
  } else {
    cat("No 'NA' values found in the grouping information.")
  }


  ##
  # Obtain the column number of the OTU_ID.
  if(id_col > 0) {
    cat("The column number for 'OTU ID Column' is: ", id_col, "\n", sep = "")
    otu <- as.data.frame(otu)            # Convert to data.frame
    rownames(otu) <- otu[, id_col]       # Rename row names
    otu <- otu[, -id_col, drop = FALSE]  # Remove the OTU ID column

  } else {}  # No OTU_ID column

  otu_t <- t(otu)
  PCA <- prcomp(otu_t, scal = TRUE)

  df_PCA_sum <- summary(PCA)
  PoA <- df_PCA_sum$importance[2, ] * 100

  PC12 <- as.data.frame(PCA$x[, 1:2])
  PC12 <- cbind(sample = rownames(PC12), PC12)
  PC12 <- merge(PC12, metadata2, by = "sample", sort = F)

  #
  colnames(PC12)[colnames(PC12) == group] <- "group"

  #
  result <- list(PCA = PC12,
                 PoA = PoA)
  return(result)
}
