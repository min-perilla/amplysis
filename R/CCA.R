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
#' @param parallel_method (character) Sample processing methods for the same group:
#' mean, sum, median, none.
#'
#' @return A table containing: (1) CCA analysis data, (2) CCA1, (3) CCA2,
#' (4) environmental factor scores, (5) Analysis of differences between
#' environmental factors and community structure.
#' @export
#'
#' @examples
#' \dontrun{
#' CCA(otu = otu, env = env, metadata = metadata, id_col = 1,
#' group = "group", parallel_method = "none")
#' }
#'
#' @importFrom vegan cca decorana envfit RsquareAdj
#' @importFrom dplyr left_join
#'
# Please use tools:: showNonASCIIfile(file.R) to check for the presence of non ASCII characters.
# tools::showNonASCIIfile(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "CCA.R"))
CCA <- function(otu, env, metadata, id_col = 1, group = "group",
                parallel_method = "none")
{

  # Check if the column names represented by the argument 'group' exist in the metadata
  if (!all(group %in% colnames(metadata))) {
    stop(paste("Some values in", ifelse(
      !all(group %in% colnames(metadata)), "group"),
      "are not present in metadata column names."))
  } else {
    cat("\033[32mgroup: `", group, "`\n\033[30m", sep = "")
  }

  ## Format check
  # Check if the metadata dataframe contains "sample", "parallel"
  if ("sample" %in% base::tolower(colnames(metadata)) &&
      "parallel" %in% base::tolower(colnames(metadata))) {
    cat("metadata --> DONE\n")
  } else {
    stop("Please ensure that the metadata table contains the `sample` column and the `parallel` column!",
         "\nsample: Sample ID (unique)",
         "\nparallel: Parallel sample identifier")
  }

  ## Process the metadata table
  # Extract columns named "sample", "parallel", the value of argument 'group1', and the value of argument 'group2'
  metadata2 <- metadata[, c("sample", "parallel", group)]

  # Discard rows with NA values
  na_rows <- apply(metadata2, 1, function(row) any(is.na(row)))
  # Output and discard rows with NA values
  if (any(na_rows)) {
    cat("The following row numbers contain NA values and have been discarded:\n")
    cat(which(na_rows), "\n")
    metadata2 <- metadata2[!na_rows, ]
  } else {
    cat("No 'NA' values found in the grouping information.\n")
  }

  ## Process the otu table
  # Discard columns based on grouping information for the otu table
  sample_values <- c(names(otu)[id_col], metadata2[["sample"]])  # Get values from the "sample" column in metadata2
  keep_columns <- colnames(otu) %in% sample_values  # Create a logical vector indicating which columns to keep
  otu2 <- otu[, keep_columns]    # Keep columns in otu where the value is TRUE
  cat("\033[32motu2 ---> DONE\n\033[30m")

  ## Process the env table
  # Discard columns based on grouping information for the env table
  sample_values_env <- c(names(env)[id_col], metadata2[["sample"]])  # Get values from the "sample" column in metadata2
  keep_columns_env <- env[["sample"]] %in% sample_values_env  # Create a logical vector indicating which columns to keep
  env2 <- env[keep_columns_env, ]    # Keep rows in env where the value is TRUE
  cat("\033[32menv2 ---> DONE\n\033[30m")


  ## Parallel sample processing
  # Define allowed methods
  allowedMethods <- base::tolower(c("mean", "sum", "median", "none"))

  # Convert to lowercase
  parallel_method <- base::tolower(parallel_method)

  # Check the parallel sample processing method
  if(!parallel_method %in% allowedMethods) {
    stop("Please enter the correct parameter for argument parallel_method:\n",
         "Process based on the `parallel` column in the `metadata` table, where samples with the same `parallel` value are considered parallel samples\n",
         "`mean`  : Take the average\n",
         "`sum`   : Take the sum\n",
         "`median`: Take the median\n",
         "`none`  : Do not process parallel samples\n")
  } else {
    cat("\033[32mParallel parallel_method: `", parallel_method, "`\n\033[30m", sep = "")
  }

  ##
  # Process parallel samples
  if (parallel_method != "none") {
    ## Convert to long format and left join with metadata2 table
    otu3 <- otu2 %>%
      # Convert the dataframe from wide format to long format
      tidyr::gather(key = "sample", value = "abun", -1) %>%
      dplyr::left_join(metadata2, by = c("sample" = "sample"))  # Left join the otu3 dataframe with metadata2 by the "sample" column
    cat("\033[32motu3 ---> DONE\n\033[30m")


    otu4 <- otu3 %>%
      # Perform grouping
      dplyr::group_by_at(dplyr::vars(names(otu3)[1], dplyr::all_of("parallel"))) %>%
      dplyr::select(names(otu3)[1], dplyr::all_of("parallel"), dplyr::all_of(group), dplyr::all_of("abun")) %>%
      dplyr::summarise_if(is.numeric, ~round(match.fun(parallel_method)(.), 1)) %>%
      dplyr::ungroup()
    cat("\033[32motu4 ---> DONE\n\033[30m")

    ## Convert back to wide format
    # Convert the long format dataframe otu4 back to wide format
    otu5 <- otu4 %>%
      tidyr::spread(key = names(otu4)[1], value = "abun")

    ##
    otu5 <- data.frame(otu5)        # Convert to dataframe
    colnames(otu5)[1] <- "#OTU ID"  # Change the name of the first column to "#OTU_ID"

    # Transpose
    otu5 <- as.data.frame(t(otu5))

    # Convert row names to the first column
    row_names <- row.names(otu5)  # Get the row names
    otu5 <- data.frame(sample = row_names, otu5)  # Add row names as a new first column to the dataframe
    row.names(otu5) <- NULL  # Reset row names

    # Use the first row as column names
    colnames(otu5) <- otu5[1, ]
    # Remove the first row
    otu5 <- otu5[-1, ]


    ## Restore original sorting
    otu6 <- otu5 %>%
      dplyr::arrange(match(otu5[[1]], otu2[[1]]))

    cat("\033[32motu6 ---> DONE\n\033[30m")
    #---------------------------------------------------------------------------


    #---------------------------------------------------------------------------
    ## Process env

    # Transpose env
    env_t = t(env2)

    # Process the row names and column names of env_t
    colnames(env_t) = env_t[1, ]  # Convert the first row into column names
    env_t = env_t[-1, ]           # Remove the first row
    cat("\033[32menv_t ---> DONE\n\033[30m")

    env_t_rownames = as.data.frame(rownames(env_t))  # Extract row names into a separate column
    colnames(env_t_rownames) = "env"      # Rename the column name
    env_t2 = cbind(env_t_rownames, env_t)  # Combine
    rownames(env_t2) <- NULL  # Reset row names
    cat("\033[32menv_t2 ---> DONE\n\033[30m")

    ## Convert to long format and left join with metadata2 table
    env_t3 <- env_t2 %>%
      # Convert the data frame from wide format to long format
      tidyr::gather(key = "sample", value = "abun", -1) %>%
      dplyr::left_join(metadata2, by = c("sample" = "sample"))  # Left join the env_t3 and metadata2 data frames by the sample column
    cat("\033[32menv_t3 ---> DONE\n\033[30m")

    # Convert the "abun" column to numeric
    env_t3["abun"] <- as.numeric(env_t3[["abun"]])

    # str(env_t3)

    env_t4 <- env_t3 %>%
      # Group by
      dplyr::group_by_at(dplyr::vars(names(env_t3)[1], dplyr::all_of("parallel"))) %>%
      dplyr::select(names(env_t3)[1], dplyr::all_of("parallel"), dplyr::all_of(group), dplyr::all_of("abun")) %>%
      dplyr::summarise_if(is.numeric, ~round(match.fun(parallel_method)(.), 2)) %>%
      dplyr::ungroup()
    cat("\033[32menv_t4 ---> DONE\n\033[30m")

    ## Convert back to wide data
    # Convert the long format data frame env_t4 back to wide format
    env_t5 <- env_t4 %>%
      tidyr::spread(key = names(env_t4)[1], value = "abun")

    ##
    env_t5 <- data.frame(env_t5)     # Convert to a data frame
    colnames(env_t5)[1] <- "env"  # Change the first column name to "env"

    # Transpose
    env_t5 <- as.data.frame(t(env_t5))

    # Convert row names to the first column
    row_names <- row.names(env_t5)  # Get the row names
    env_t5 <- data.frame(sample = row_names, env_t5)  # Add the row names as the new first column in the data frame
    row.names(env_t5) <- NULL  # Reset row names

    # Use the first row as column names
    colnames(env_t5) <- env_t5[1, ]
    # Remove the first row
    env_t5 <- env_t5[-1, ]
    cat("\033[32menv_t5 ---> DONE\n\033[30m")

    ## Restore the original order
    env_t6 <- env_t5 %>%
      dplyr::arrange(match(env_t5[[1]], env_t2[[1]]))

    cat("\033[32menv_t6 ---> DONE\n\033[30m")


    # Transpose env
    env6 = t(env_t6)
    # Process the row names and column names of env6
    colnames(env6) = env6[1, ]  # Convert the first row into column names
    env6 = env6[-1, ]           # Remove the first row

    env_t_rownames = as.data.frame(rownames(env6))  # Extract row names into a separate column
    colnames(env_t_rownames) = "sample"      # Rename the column name
    env6 = cbind(env_t_rownames, env6)  # Combine
    rownames(env6) <- NULL  # Reset row names
    cat("\033[32menv6 ---> DONE\n\033[30m")


    #---------------------------------------------------------------------------
    ## Synchronize metadata
    metadata3 = metadata2 %>%
      # Keep the columns "parallel" and the column represented by group
      dplyr::select(dplyr::all_of("parallel"), dplyr::all_of(group)) %>%
      # Remove duplicates in the parallel column
      dplyr::distinct(parallel, .keep_all = TRUE)
    # Rename the first column to "sample"
    colnames(metadata3)[1] <- "sample"

  } else {
    # Do not process parallel samples

    otu6 = otu2
    env6 = env2
    metadata3 = metadata2

    cat("\033[32motu6 ---> DONE2\n\033[30m")
    cat("\033[32menv6 ---> DONE2\n\033[30m")
    cat("\033[32mmetadata3 ---> DONE2\n\033[30m")
  }


  ##
  # Obtain the column number of the OTU_ID.
  if(id_col > 0) {
    cat("The column number for 'OTU ID Column' is: ", id_col, "\n", sep = "")

    otu6 <- as.data.frame(otu6)            # Convert to data.frame
    env6 <- as.data.frame(env6)            # Convert to data.frame

    rownames(otu6) <- otu6[, id_col]       # Rename row names
    rownames(env6) <- env6[, id_col]       # Rename row names

    otu6 <- otu6[, -id_col, drop = FALSE]  # Remove the OTU ID column
    env6 <- env6[, -id_col, drop = FALSE]  # Remove the OTU ID column

    cat("\033[32mid_col ---> DONE\n\033[30m")
  } else {}  # No OTU_ID column


  ## Convert to numeric
  otu6[] <- lapply(otu6, function(x) as.numeric(trimws(x)))
  env6[] <- lapply(env6, function(x) as.numeric(trimws(x)))


  ##
  # DCA Detection
  # Using the decorana function to check if the data is suitable for RDA analysis
  # Selection based on the Axis Lengths value of DCA1
  # If >4.0, choose CCA
  # If between 3.0-4.0, both RDA and CCA are acceptable
  # If <3.0, choose RDA
  DCA <- vegan::decorana(t(otu6))
  print(DCA)
  cat("\033[31mWhen `Axis lengths` value > 4.0, choose CCA;
When `Axis lengths` value is between 3.0 and 4.0, either CCA or RDA can be chosen;
When `Axis lengths` value < 3.0, choose RDA;\n\033[0m")


  ##
  # The 'cca' function from the vegan package supports CCA analysis
  df_cca <- vegan::cca(t(otu6) ~ ., env6, scale = T)

  # View CCA results, using Type I scaling as an example, as described in the reference paper
  scaling1 <- summary(df_cca, scaling = 1)


  ##
  # Adjusted R2 value
  R2 <- vegan::RsquareAdj(df_cca)
  # R2_noadj <- R2$r.squared     # Original R2
  R2_adj <- R2$adj.r.squared   # Adjusted R2

  # Calculate the explained variance for the constrained axes after adjusting R2
  R2_adj_exp <- R2_adj * df_cca$CCA$eig / sum(df_cca$CCA$eig)
  # Calculate axis label data
  CCA1 <- paste0("CCA1(", round(R2_adj_exp[1] * 100, 1), "%)")
  CCA2 <- paste0("CCA2(", round(R2_adj_exp[2] * 100, 1), "%)")


  ## Permutation test for constrained axes and p-value correction
  # test <- vegan::anova.cca(df_cca, permutations = 999)  # Permutation test for all constrained axes, i.e., global test, based on 999 permutations, see ?anova.cca
  # test_axis <- anova.cca(df_cca, by = 'axis', permutations = 999)  # Test each constrained axis individually, based on 999 permutations
  # test_axis$`Pr(>F)` <- p.adjust(test_axis$`Pr(>F)`, method = 'bonferroni')  # p-value correction (example: bonferroni)


  ##
  # Extract plot data
  sites <- data.frame(scaling1$sites)[1:2]

  # Add group information
  sites <- data.frame(sites, rownames(sites))
  colnames(sites) = c("CCA1", "CCA2", "sample")  # Rename columns
  sites <- dplyr::left_join(sites, metadata3, by = c("sample" = "sample"))  # Merge plot data with group information


  ## Standardize the group column name in the 'sites' data to "group"
  if (group %in% colnames(sites)) {
    colnames(sites)[colnames(sites) == group] <- "group"
  }


  ## Extract other data
  df_env <- data.frame(scaling1$biplot)[1:2]   # Extract environmental factor scores


  ## Analysis of the differences between environmental factors and community structure: significance calculation
  # cca_sum <- summary(df_cca)  # Descriptive statistics
  # Test the significance of environmental factors (Monte Carlo permutation test)
  # df_permutest <- permutest(df_cca, permu = 999)  # permu = 999 indicates the number of permutation cycles
  # Significance test for each environmental factor
  df_envfit <- vegan::envfit(df_cca, env6, permu = 999)
  # Data processing
  # cor_data <- data.frame(cca_sum$constr.chi / cca_sum$tot.chi, cca_sum$unconst.chi / cca_sum$tot.chi)
  cor_com <- data.frame(tax = colnames(env6), r = df_envfit$vectors$r, p = df_envfit$vectors$pvals)
  # Mark p < 0.05 as FALSE, p > 0.05 as TRUE, use this data to plot the bar chart.
  cor_com[1:5, 3] = cor_com[, 3] > 0.05


  # CCA Data
  # Axis Eigenvector 1
  # Axis Eigenvector 2
  # Extract Environmental Factor Scores
  # Analysis of differences between environmental factors and community structure: Significance calculation
  result <- NULL
  result <- list(
    data = sites,      # RDA data
    CCA1 = CCA1,       # Axis eigenvalue 1
    CCA2 = CCA2,       # Axis eigenvalue 2
    env = df_env,      # Extract environmental factor scores
    cor_com = cor_com  # Analysis of differences between environmental factors and community structure: significance calculation
  )

  ##
  cat("\033[32m--- Please use the `CCA_plot()` function for visualization. ---\n\033[0m")

  return(result)
}
