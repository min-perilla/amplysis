#' @title Venn Diagram Analysis
#'
#' @description
#' When the number of samples is less than 5, a Venn diagram can effectively
#' display the intersections of OTUs among the samples. Although the `venn()`
#' function allows the `metadata` parameter to be NULL, we still recommend using
#' the `group` column in the metadata file to control the grouping information.
#'
#' It is important to note that when the number of samples is 5 or more, we
#' recommend using a `upset` diagram to represent the intersections of OTUs.
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
#' venn(otu, metadata, 1, "group", "mean")
#' }
#'
#' @importFrom dplyr group_by_at select ungroup vars summarise_if
#' @importFrom tidyr gather spread
#' @importFrom VennDiagram get.venn.partitions
#'
venn = function(otu, metadata = NULL, id_col = 1,
         group = "group", parallel_method = "mean")
{
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


    ##
    otu3 <- otu2 %>%
      tidyr::gather(key = "sample", value = "abun", -1) %>%
      dplyr::left_join(metadata2, by = c("sample" = "sample"))
    cat("\033[32motu3 ---> DONE\n\033[30m")


    ##
    if (parallel_method != "none") {
      otu4 <- otu3 %>%
        # 进行分组
        dplyr::group_by_at(dplyr::vars(names(otu3)[1], group)) %>%
        dplyr::select(names(otu3)[1], group, "abun") %>%
        dplyr::summarise_if(is.numeric, match.fun(parallel_method)) %>%
        dplyr::ungroup()
      cat("\033[32motu4 ---> DONE\n\033[30m")
    }


    ##
    otu5 <- otu4 %>%
      tidyr::spread(key = names(otu3)[1], value = "abun")


    ##
    otu5 <- data.frame(otu5)
    colnames(otu5)[1] <- "#OTU_ID"
    otu5 <- as.data.frame(t(otu5))
    row_names <- row.names(otu5)
    otu5 <- data.frame(sample = row_names, otu5)
    row.names(otu5) <- NULL
    colnames(otu5) <- otu5[1, ]
    otu5 <- otu5[-1, ]


  } else {
    #
    otu5 <- otu
  }


  ##
  # Obtain the column number of the OTU_ID.
  if(id_col > 0) {
    cat("The column number for 'OTU ID Column' is: ", id_col, "\n", sep = "")
    otu5 <- as.data.frame(otu5)            # Convert to data.frame
    rownames(otu5) <- otu5[, id_col]       # Rename row names
    otu5 <- otu5[, -id_col, drop = FALSE]  # Remove the OTU ID column

  } else {}  # No OTU_ID column


  ##
  df <- NULL
  for (i in 1:length(colnames(otu5))){
    group <- colnames(otu5)[i]
    df[[group]] <- rownames(otu5)[which(otu5[,i]!= 0)]
  }

  inter <- VennDiagram::get.venn.partitions(df)

  for (i in 1:nrow(inter)){
    inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
  }

  inter2 <- inter[, setdiff(names(inter), c("..set..", "..values.."))]
  utils::write.csv(inter2, file = "Venn_inter.csv", row.names = FALSE)


  ##
  return(otu5)
}
