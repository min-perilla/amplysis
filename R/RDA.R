#' @title Redundancy Analysis (RDA)
#' @description
#' Redundancy Analysis (RDA) is a statistical method used to explore the
#' relationship between response variables and explanatory variables in
#' multivariate datasets. This function references the `rda()` function in the
#' `vegan` package for analysis. Prior to conducting RDA analysis, it is
#' recommended to perform Detrended Correspondence Analysis (DCA) detection
#' using the `decorana` function from the `vegan` package. Selection is based
#' on the Axis Lengths value of DCA1 as follows: (1) If >4.0, choose CCA;
#' (2) If between 3.0-4.0, both RDA and CCA are acceptable;
#' (3) If <3.0, choose RDA.
#'
#' @param otu otu table
#' @param env Environment Factor Data Table
#' @param metadata metadata table containing grouping information.
#' @param id_col (integer) The column number of the OTU ID column in the
#' OTU table, by default, is 1.
#' @param group (character) Group 1, please enter the column name of
#' the grouping information in the metadata table.
#'
#' @return A table containing: (1) RDA analysis data, (2) RDA1, (3) RDA2,
#' (4) species scores, (5) environmental factor scores, (6) Analysis of
#' differences between environmental factors and community structure.
#' @export
#'
#' @examples
#' \dontrun{
#' RDA(otu, env, metadata, 1, "group")
#' }
#'
#' @importFrom vegan decorana envfit rda
#' @importFrom dplyr left_join
#'
# Please use tools:: showNonASCIIfile(file.R) to check for the presence of non ASCII characters.
# tools::showNonASCIIfile(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "RDA.R"))
RDA <- function(otu, env, metadata, id_col = 1, group = "group")
{
  ##
  if (!group %in% colnames(metadata)) {
    stop(paste("Some values in",
               ifelse(!all(group %in% colnames(metadata)), "group"),
               "are not present in metadata column names."))
  } else {
    cat("\033[32mgroup1: `", group, "`\n", "\033[0m", sep = "")
  }


  ##
  if ("sample" %in% base::tolower(colnames(metadata))) {
    cat("metadata --> DONE\n")
  } else {
    stop("Please ensure that the metadata table contains the `sample` column!",
         "\nsample: Sample ID (unique)")
  }


  ##
  metadata2 <- metadata[, c("sample", group)]
  na_rows <- apply(metadata2, 1, function(row) any(is.na(row)))
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
    env <- as.data.frame(env)            # Convert to data.frame

    rownames(otu) <- otu[, id_col]       # Rename row names
    rownames(env) <- env[, id_col]       # Rename row names

    otu <- otu[, -id_col, drop = FALSE]  # Remove the OTU ID column
    env <- env[, -id_col, drop = FALSE]  # Remove the OTU ID column

  } else {}  # No OTU_ID column


  ##
  # DCA Detection
  # Using the decorana function to check if the data is suitable for RDA analysis
  # Selection based on the Axis Lengths value of DCA1
  # If >4.0, choose CCA
  # If between 3.0-4.0, both RDA and CCA are acceptable
  # If <3.0, choose RDA
  DCA <- vegan::decorana(t(otu))
  print(DCA)
  cat("\033[31mWhen `Axis lengths` value > 4.0, choose CCA;
When `Axis lengths` value is between 3.0 and 4.0, either CCA or RDA can be chosen;
When `Axis lengths` value < 3.0, choose RDA;\n\033[0m")

  ##
  df_rda <- vegan::rda(t(otu), env, scale = T)

  #
  df_rda2 <- data.frame(df_rda$CCA$u[, 1:2], rownames(env))
  colnames(df_rda2) = c("RDA1", "RDA2", "sample")
  df_rda2 <- dplyr::left_join(df_rda2, metadata2, by = c("sample" = "sample"))


  ##
  df_rda_score <- data.frame(df_rda$CCA$v[,1:2])
  df_rda_env <- as.data.frame(df_rda$CCA$biplot[,1:2])


  ##
  RDA1 = round(df_rda$CCA$eig[1] / sum(df_rda$CCA$eig) * 100, 2)
  RDA2 = round(df_rda$CCA$eig[2] / sum(df_rda$CCA$eig) * 100, 2)


  ##
  df_envfit <- vegan::envfit(df_rda, env, permu = 999)
  cor_com <- data.frame(tax = colnames(env), r = df_envfit$vectors$r,
                        p = df_envfit$vectors$pvals)
  cor_com[1:5,3] = cor_com[,3] > 0.05


  # RDA Data
  # Axis Eigenvector 1
  # Axis Eigenvector 2
  # Extract Species Scores
  # Extract Environmental Factor Scores
  result <- NULL
  result <- list(
    data = df_rda2,
    RDA1 = RDA1,
    RDA2 = RDA2,
    score = df_rda_score,
    env = df_rda_env,
    cor_com = cor_com
  )

  ##
  cat("\033[32m--- Please use the `RDA_plot()` function for visualization. ---\n\033[0m")

  return(result)
}
