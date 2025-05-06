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
#' @param replicate_method (character) Sample processing methods for the same group:
#' mean, sum, median, none.
#'
#' @return A list (containing plot data and the stress values).
#' @export
#'
#' @examples
#' \dontrun{nmds(otu = otu, metadata = metadata, id_col = 1, group = "group",
#'               replicate_method = "none")}
#'
#' @importFrom vegan vegdist metaMDS stressplot
#'
# tools::showNonASCIIfile(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "nmds.R"))
nmds <- function(otu, metadata, id_col = 1, group = "group", replicate_method = "none")
{
  # Check if the column specified by 'group' exists in metadata
  if (!all(group %in% colnames(metadata))) {
    stop(paste("Some values in", ifelse(
      !all(group %in% colnames(metadata)), "group"),
      "are not present in metadata column names."))
  } else {
    cat("\033[32mgroup: `", group, "`\n", sep = "")
  }

  ## Format check
  # Check if the metadata dataframe contains "sample" and "replicate"
  if ("sample" %in% base::tolower(colnames(metadata)) &&
      "replicate" %in% base::tolower(colnames(metadata))) {
    cat("metadata --> DONE\n")
  } else {
    stop("Please ensure that the metadata table contains the `sample` column and the `replicate` column!",
         "\nsample: Sample ID (unique)",
         "\nreplicate: replicate sample identifier")
  }

  ## Process metadata table
  # Extract columns named "sample", "replicate", group, and group1
  metadata2 <- metadata[, c("sample", "replicate", group)]

  # Remove rows with NA values
  na_rows <- apply(metadata2, 1, function(row) any(is.na(row)))
  # Output and remove rows with NA values
  if (any(na_rows)) {
    cat("The following row numbers contain NA values and have been discarded:\n")
    cat(which(na_rows), "\n")
    metadata2 <- metadata2[!na_rows, ]
  } else {
    cat("No 'NA' values found in the grouping information.\n")
  }

  ## Process OTU table
  # Remove columns in OTU based on the grouping information
  sample_values <- c(names(otu)[1], metadata2[["sample"]])  # Get the "sample" column values from metadata2
  keep_columns <- colnames(otu) %in% sample_values  # Create a logical vector indicating which columns to keep
  otu2 <- otu[, keep_columns]    # Keep the columns in OTU that correspond to TRUE in the logical vector

  ## replicate sample processing
  # Define allowed methods
  allowedMethods <- base::tolower(c("mean", "sum", "median", "none"))

  # Convert to lowercase
  replicate_method <- base::tolower(replicate_method)

  # Check the replicate sample processing method
  if(!replicate_method %in% allowedMethods) {
    stop("Please provide a valid replicate_method parameter:\n",
         "Based on the `replicate` column in `metadata`, samples with the same `replicate` value are treated as replicate samples\n",
         "`mean`  : take the mean\n",
         "`sum`   : sum the values\n",
         "`median`: take the median\n",
         "`none`  : no replicate sample processing\n")
  } else {
    cat("\033[32mreplicate replicate_method: `", replicate_method, "`\n\033[30m", sep = "")
  }

  ##
  # Process replicate samples
  if (replicate_method != "none") {
    ## Convert to long format and left join with metadata2
    otu3 <- otu2 %>%
      # Convert from wide format to long format
      tidyr::gather(key = "sample", value = "abun", -1) %>%
      dplyr::left_join(metadata2, by = c("sample" = "sample"))  # Left join OTU3 with metadata2 by "sample" column
    cat("\033[32motu3 ---> DONE\n\033[30m")

    otu4 <- otu3 %>%
      # Perform grouping
      dplyr::group_by_at(dplyr::vars(names(otu3)[1], dplyr::all_of("replicate"))) %>%
      dplyr::select(names(otu3)[1], dplyr::all_of("replicate"), dplyr::all_of(group), dplyr::all_of("abun")) %>%
      dplyr::summarise_if(is.numeric, ~round(match.fun(replicate_method)(.), 1)) %>%
      dplyr::ungroup()
    cat("\033[32motu4 ---> DONE\n\033[30m")

    ## Convert back to wide format
    # Convert the long-format dataframe otu4 back to wide format
    otu5 <- otu4 %>%
      tidyr::spread(key = names(otu4)[1], value = "abun")

    ##
    otu5 <- data.frame(otu5)        # Convert to data frame
    colnames(otu5)[1] <- "#OTU ID"  # Rename the first column to "#OTU ID"

    # Transpose
    otu5 <- as.data.frame(t(otu5))

    # Convert row names to the first column
    row_names <- row.names(otu5)  # Get row names
    otu5 <- data.frame(sample = row_names, otu5)  # Add row names as the first column in the data frame
    row.names(otu5) <- NULL  # Reset row names

    # Use the first row as column names
    colnames(otu5) <- otu5[1, ]
    # Remove the first row
    otu5 <- otu5[-1, ]

    ## Restore original sorting
    otu6 <- otu5 %>%
      dplyr::arrange(match(otu5[[1]], otu2[[1]]))

    cat("\033[32motu6 ---> DONE\n\033[30m")

    ## Synchronize metadata
    metadata3 = metadata2 %>%
      # Keep columns named "replicate" and group
      dplyr::select(dplyr::all_of("replicate"), dplyr::all_of(group)) %>%
      # Remove duplicates based on the "replicate" column
      dplyr::distinct(replicate, .keep_all = TRUE)
    # Rename the first column to "sample"
    colnames(metadata3)[1] <- "sample"

  } else {
    # No replicate sample processing
    otu6 = otu2
    metadata3 = metadata2
    cat("\033[32motu6 ---> DONE2\n\033[30m")
    cat("\033[32mmetadata3 ---> DONE2\n\033[30m")
  }

  ##
  # Obtain the column number of the OTU_ID.
  if(id_col > 0) {
    cat("The column number for 'OTU ID Column' is: ", id_col, "\n", sep = "")
    otu6 <- as.data.frame(otu6)            # Convert to data.frame
    rownames(otu6) <- otu6[, id_col]       # Rename row names
    otu6 <- otu6[, -id_col, drop = FALSE]  # Remove the OTU ID column
    cat("\033[32mid_col ---> DONE\n\033[30m")
  } else {}  # No OTU_ID column


  ## Convert to numeric
  otu6[] <- lapply(otu6, function(x) as.numeric(trimws(x)))

  ##
  otu_t <- t(otu6)  # Transpose OTU table

  # Remove constant or all-zero columns
  otu_t <- otu_t[, apply(otu_t, 2, function(x) var(x) != 0)]

  otu.distance <- vegan::vegdist(otu_t, method = "bray")        # Calculate Bray-Curtis distance

  df_NMDS <- vegan::metaMDS(otu.distance, k = 2, trymax = 20)   # NMDS ordination, focusing on stress, points, and species indices
  Stress <- df_NMDS$stress
  cat("Generally, a Stress value < 0.2 is usable; < 0.1 is good\n",
      "Stress = ", Stress,
      if(Stress < 0.2 && Stress > 0.1) {
        paste0("\033[34m", "(Usable)", "\033[39m")
      } else if(Stress < 0.1) {
        paste0("\033[32m", "(Good)", "\033[39m")
      } else {
        paste0("\033[31m", "(Poor)", "\033[39m")
      }, "\n", sep = "")
  vegan::stressplot(df_NMDS)  # Check the relationship between dissimilarity and ordination distance; if points are far from the line, NMDS analysis is suitable


  # Plotting data
  NMDS <- as.data.frame(df_NMDS$points)                    # Extract data for plotting
  NMDS <- cbind(sample = rownames(NMDS), NMDS)             # Add row names as a separate column
  NMDS <- merge(NMDS, metadata3, by = "sample", sort = F)  # Merge plotting data with grouping information

  # Rename columns
  colnames(NMDS)[colnames(NMDS) == group] <- "group"

  # Return results
  result <- list(NMDS = NMDS,     # Plotting data
                 stress = Stress  # Stress value
  )
  return(result)
}

