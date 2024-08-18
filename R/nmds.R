#' @title Non-metric Multi-dimensional Scaling
#' @description
#' In bioinformatics analysis, Non-Metric Multidimensional Scaling (NMDS) is
#' used to reduce the dimensionality of high-dimensional data and represent the
#' similarity or dissimilarity between samples in a lower-dimensional space.
#' The "nmds()" function utilizes "vegan::vegdist" to calculate the Bray-Curtis
#' distance and employs "vegan::metaMDS" for conducting NMDS ordination analysis.
#'
#' @param otu otu table
#' @param metadata metadata table
#' @param id_col (integer) The OTU_ID column is in which column, defaulting to 0
#' means there is no OTU_ID column, and the data is already numeric.
#' @param group (Required, character) Grouping information. please enter the
#' column name of the grouping information in the metadata table.
#'
#' @return A list (containing plot data and the stress values).
#' @export
#'
#' @examples
#' \dontrun{nmds(otu, metadata, 1, "group2")}
#'
#' @importFrom vegan vegdist metaMDS stressplot
#'
nmds <- function(otu, metadata, id_col = 1, group = "group")
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


  ##
  otu_t <- t(otu)  # 转置 OTU 表
  otu.distance <- vegan::vegdist(otu_t, method = "bray")
  df_NMDS <- vegan::metaMDS(otu.distance, k = 2, trymax = 20)
  Stress <- df_NMDS$stress
  cat("In general, Stress values < 0.2 are acceptable; < 0.1 indicates better performance.\n",
      "Stress = ", Stress,
      if(Stress < 0.2 && Stress > 0.1) {
        paste0("\033[34m", " (Acceptable)", "\033[39m")
      } else if(Stress < 0.1) {
        paste0("\033[32m", " (Good)", "\033[39m")
      } else {
        paste0("\033[31m", " (Poor)", "\033[39m")
      }, "\n", sep = "")
  vegan::stressplot(df_NMDS)


  #
  NMDS <- as.data.frame(df_NMDS$points)
  NMDS <- cbind(sample = rownames(NMDS), NMDS)
  NMDS <- merge(NMDS, metadata2, by = "sample", sort = F)

  #
  colnames(NMDS)[colnames(NMDS) == group] <- "group"

  #
  result <- list(NMDS = NMDS,
                 stress = Stress
  )
  return(result)
}
