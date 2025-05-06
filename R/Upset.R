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
#' @param replicate_method (character) Sample processing methods for the same group:
#' mean, sum, median, none.
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
                  group = "group", replicate_method = "mean")
{
  ##
  # Allow metadata to be NULL. When metadata is NULL, the group parameter will be invalid
  if(!is.null(metadata)) {

    # Check if the column name specified by the group parameter exists in metadata
    if (!all(group %in% colnames(metadata))) {
      stop(paste("Some values in", ifelse(
        !all(group %in% colnames(metadata)), "group"),
        "are not present in metadata column names."))
    } else {
      cat("\033[32mgroup: `", group, "`\n\033[30m", sep = "")
    }


    ## Format check
    # Check if metadata dataframe contains "sample" and "replicate"
    if ("sample" %in% base::tolower(colnames(metadata)) &&
        "replicate" %in% base::tolower(colnames(metadata))) {
      cat("\033[32mmetadata --> DONE\n\033[30m")
    } else {
      stop("Please ensure that the metadata table contains the `sample` column and the `replicate` column!",
           "\nsample: Sample ID (unique)",
           "\nreplicate: replicate sample identifier")
    }


    ## Process the metadata table
    # Extract columns "sample", "replicate", the value of the group parameter, and other relevant columns from metadata
    metadata2 <- metadata[, c("sample", "replicate", group)]

    # For metadata, discard rows containing NA values
    na_rows <- apply(metadata2, 1, function(row) any(is.na(row)))

    # Output and discard rows containing NA values
    if (any(na_rows)) {
      cat("metadata: The following row numbers contain NA values and have been discarded:\n")
      cat(which(na_rows), "\n")

      # Discard corresponding rows from metadata
      metadata2 <- metadata2[!na_rows, ]

    } else {
      cat("metadata: No 'NA' values found in the grouping information.")
    }

    # For otu, discard columns based on group information
    sample_values <- c(names(otu)[1], metadata2[["sample"]])  # Get the values from the "sample" column in metadata2
    keep_columns <- colnames(otu) %in% sample_values  # Create a logical vector to identify columns to keep
    otu2 <- otu[, keep_columns]             # Keep columns with TRUE in the logical vector


    ## replicate sample handling
    # Define allowed methods
    allowedMethods <- base::tolower(c("mean", "sum", "median", "none"))

    # Convert to lowercase
    replicate_method <- base::tolower(replicate_method)

    # Check the replicate sample handling method
    if(!replicate_method %in% allowedMethods) {
      stop("Please input a valid replicate_method parameter:\n",
           "For handling replicate samples based on the `replicate` column in the metadata table, where samples with the same `replicate` value are considered replicate samples\n",
           "`mean`  : Take the mean\n",
           "`sum`   : Take the sum\n",
           "`median`: Take the median\n",
           "`none`  : Do not process replicate samples\n")
    } else {
      cat("\033[32mreplicate replicate_method: `", replicate_method, "`\n\033[30m", sep = "")
    }



    ##
    # Process replicate samples
    if (replicate_method != "none") {
      ## Convert to long format and left join with metadata2
      otu3 <- otu2 %>%
        # Convert the dataframe from wide to long format
        tidyr::gather(key = "sample", value = "abun", -1) %>%
        dplyr::left_join(metadata2, by = c("sample" = "sample"))  # Left join otu3 with metadata2 on sample column
      cat("\033[32motu3 ---> DONE\n\033[30m")


      otu4 <- otu3 %>%
        # Group by replicate and the first column
        dplyr::group_by_at(dplyr::vars(names(otu3)[1], dplyr::all_of("replicate"))) %>%
        dplyr::select(names(otu3)[1], dplyr::all_of("replicate"), dplyr::all_of(group), dplyr::all_of("abun")) %>%
        dplyr::summarise_if(is.numeric, ~round(match.fun(replicate_method)(.), 1)) %>%
        dplyr::ungroup()
      cat("\033[32motu4 ---> DONE\n\033[30m")

      ## Convert back to wide format
      # Convert the long format otu4 dataframe back to wide format
      otu5 <- otu4 %>%
        tidyr::spread(key = names(otu4)[1], value = "abun")

      ##
      otu5 <- data.frame(otu5)        # Convert to dataframe
      colnames(otu5)[1] <- "#OTU ID"  # Rename the first column to "#OTU ID"

      # Transpose
      otu5 <- as.data.frame(t(otu5))

      # Convert row names to the first column
      row_names <- row.names(otu5)  # Get row names
      otu5 <- data.frame(sample = row_names, otu5)  # Add row names as the first column in the dataframe
      row.names(otu5) <- NULL  # Reset row names

      # Use the first row as column names
      colnames(otu5) <- otu5[1, ]
      # Remove the first row
      otu5 <- otu5[-1, ]


      ## Restore original order
      otu6 <- otu5 %>%
        dplyr::arrange(match(otu5[[1]], otu2[[1]]))

      cat("\033[32motu6 ---> DONE\n\033[30m")


      ## Sync metadata
      metadata3 = metadata2 %>%
        # Keep columns "replicate" and the group column
        dplyr::select(dplyr::all_of("replicate"), dplyr::all_of(group)) %>%
        # Remove duplicates in the replicate column
        dplyr::distinct(replicate, .keep_all = TRUE)
      # Rename the first column to "sample"
      colnames(metadata3)[1] <- "sample"



    } else {
      # Do not process replicate samples

      otu6 = otu2
      metadata3 = metadata2
      cat("\033[32motu6 ---> DONE2\n\033[30m")
      cat("\033[32mmetadata3 ---> DONE2\n\033[30m")
    }



    ###
  } else {
    # When metadata is NULL
    otu6 <- otu
    cat("\033[32m metadata is NULL\n\033[30m")
  }


  ##
  # Obtain the column number of the OTU_ID.
  if(id_col > 0) {
    cat("The column number for 'OTU ID Column' is: ", id_col, "\n", sep = "")
    otu6 <- as.data.frame(otu6)            # Convert to data.frame
    rownames(otu6) <- otu6[, id_col]       # Rename row names
    otu6 <- otu6[, -id_col, drop = FALSE]  # Remove the OTU ID column

  } else {}  # No OTU_ID column

  # Convert to numeric
  otu6[] <- lapply(otu6, function(x) as.numeric(trimws(x)))

  ##
  return(otu6)
}

