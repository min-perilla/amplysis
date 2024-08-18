#' @title Canonical Correspondence Analysis (CCA)
#' @description
#' Canonical Correspondence Analysis (CCA) is a multivariate statistical method
#' used to study the relationship between two sets of variables. It can reveal
#' the inherent connections between the two sets of variables. This function
#' references the `cca()` function in the `vegan` package for analysis. Prior to
#' conducting RDA analysis, it is recommended to perform Detrended
#' Correspondence Analysis (DCA) detection using the `decorana` function from
#' the `vegan` package. Selection is based on the Axis Lengths value of DCA1 as
#' follows: (1) If >4.0, choose CCA; (2) If between 3.0-4.0, both RDA and CCA
#' are acceptable; (3) If <3.0, choose RDA.
#'
#' @param otu otu table
#' @param env Environment Factor Data Table
#' @param metadata metadata table containing grouping information.
#' @param id_col (integer)The column number of the OTU ID column in the
#' OTU table, by default, is 1.
#' @param group (character) Group 1, please enter the column name of
#' the grouping information in the metadata table.
#'
#' @return A table containing: (1) CCA analysis data, (2) CCA1, (3) CCA2,
#' (4) environmental factor scores, (5) Analysis of differences between
#' environmental factors and community structure.
#' @export
#'
#' @examples
#' \dontrun{
#' CCA(otu, env, metadata, 1, "group")
#' }
#'
#' @importFrom vegan cca decorana envfit RsquareAdj
#' @importFrom dplyr left_join
#'
# Please use tools:: showNonASCIIfile(file.R) to check for the presence of non ASCII characters.
# tools::showNonASCIIfile(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "CCA.R"))
CCA <- function(otu, env, metadata, id_col = 1, group = "group")
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
  colnames(metadata2) = c("sample", "group")

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
  df_cca <- vegan::cca(t(otu) ~ ., env)
  scaling1 <- summary(df_cca, scaling = 1)


  ##
  R2 <- vegan::RsquareAdj(df_cca)
  R2_adj <- R2$adj.r.squared

  #
  R2_adj_exp <- R2_adj * df_cca$CCA$eig/sum(df_cca$CCA$eig)

  #
  CCA1 <- paste0("CCA1(",round(R2_adj_exp[1]*100, 1),"%)")
  CCA2 <- paste0("CCA2(",round(R2_adj_exp[2]*100, 1),"%)")


  ##
  #
  sites <- data.frame(scaling1$sites)[1:2]

  #
  sites <- data.frame(sites, rownames(sites))
  colnames(sites) = c("CCA1", "CCA2", "sample")
  sites <- dplyr::left_join(sites, metadata2, by = c("sample" = "sample"))

  #
  df_env <- data.frame(scaling1$biplot)[1:2]


  ##
  df_envfit <- vegan::envfit(df_cca, env, permu = 999)
  cor_com <- data.frame(tax = colnames(env), r = df_envfit$vectors$r,
                        p = df_envfit$vectors$pvals)
  cor_com[1:5,3] = cor_com[,3] > 0.05


  # CCA Data
  # Axis Eigenvector 1
  # Axis Eigenvector 2
  # Extract Environmental Factor Scores
  # Analysis of differences between environmental factors and community structure: Significance calculation
  result <- NULL
  result <- list(
    data = sites,
    CCA1 = CCA1,
    CCA2 = CCA2,
    env = df_env,
    cor_com = cor_com
  )

  ##
  cat("\033[32m--- Please use the `CCA_plot()` function for visualization. ---\n\033[0m")

  return(result)
}
