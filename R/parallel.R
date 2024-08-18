#' @title Parallel Sample Processing
#'
#' @description
#' Parallel samples are determined using the "parallel" column in the sample
#' metadata file. Samples with the same value in the "parallel" column are
#' considered parallel samples. The methods for handling parallel samples are:
#' (1) "mean" (calculate the average); (2) "sum" (calculate the sum);
#' (3) "median" (calculate the median); and (4) "none" (no operation).
#' The default method is "mean" (calculate the average).
#'
#' @param otu Feature table
#' @param metadata Sample metadata file
#' @param id_col (integer)The column number of the OTU ID column in the
#' OTU table, by default, is 1.
#' @param group (character) Group 1, please enter the column name of
#' the grouping information in the metadata table.
#' @param parallel_method  (character) Sample processing methods for the same group:
#' average, sum, median, none.
#' @param digits (integer) When the parallel sample method is set to "mean," the
#' number of decimal places will default to 0.
#'
#' @param metadata_out (logical) The returned result will be a list containing
#' both OTU and metadata.
#'
#' @return feature table
#' @export
#'
#' @examples
#' \dontrun{
#' otu_new = parallel(otu, metadata, 1, "group", "mean")
#' otu_metadata_list = parallel(otu, metadata, 1, "group", "mean", 0, T)
#' }
#'
#' @importFrom dplyr group_by group_cols summarize_at
#' @importFrom ggplot2 vars
#'
parallel = function(otu, metadata, id_col = 1, group = "group",
         parallel_method = "mean", digits = 0, metadata_out = F)
{
  {
    sample_flag = 0
    parallel_flag = 0
    group_flag = 0
  }

  {
    if("sample" %in% base::tolower(colnames(metadata))){ sample_flag = 1 }
    if("parallel" %in% base::tolower(colnames(metadata))){ parallel_flag = 1 }
    if(group %in% colnames(metadata)){ group_flag = 1 }
  }

  #
  if(sample_flag == 0 || parallel_flag == 0 || group_flag == 0){
    stop(
      "Please ensure that the metadata table contains the",
      if(sample_flag == 0) { " `sample`" } else { NULL },
      if(parallel_flag == 0) { ", `parallel`" } else { NULL },
      if(group_flag == 0) { ", `group`" } else { NULL },
      " column.\n",

      ##
      if(sample_flag == 0) { "\nsample: Sample ID (unique)" } else { NULL },
      if(parallel_flag == 0) { "\nparallel: Parallel sample identifier" } else { NULL },
      if(group_flag == 0) { "\ngroup: Group information" } else { NULL }, "\n"
    )
  } else { cat("\033[32mmetadata --> DONE\n\033[30m", sep = "") }


  ##
  metadata2 <- metadata[, c("sample", "parallel", group)]
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
    otu2 = as.data.frame(otu)             # Convert to data.frame
    rownames(otu2) <- otu2[, id_col]      # Rename row names
    otu2 = otu2[, -id_col, drop = FALSE]  # Remove the OTU ID column

  } else {
    otu2 = otu
  }  # No OTU_ID column


  ##
  allowedMethods <- base::tolower(c("mean", "sum", "median", "none"))
  parallel_method = base::tolower(parallel_method)
  #
  if(!base::tolower(parallel_method) %in% allowedMethods) {
    stop("Please enter the correct argument for the parameter 'parallel_method':\n",
         "Process according to the 'parallel' column in the `metadata` table, ",
         "\nsamples with the same 'parallel' value are considered parallel samples.\n",
         "`mean`  : Calculate the average\n",
         "`sum`   : Calculate the sum\n",
         "`median`: Calculate the median\n",
         "`none`  : Do not process parallel samples\n")


  } else {
    cat("\033[32mparallel_method: `", parallel_method, "`\n\033[30m", sep = "")
  }


  ##
  if(parallel_method %in% c("mean", "sum", "median")) {
    otu_t <- as.data.frame(t(otu2))
    otu_t_g <- merge(otu_t, metadata2, by.x = "row.names", by.y = "sample",
                     all.x = TRUE, sort = F)
    row.names(otu_t_g) <- otu_t_g[, 1]
    otu_t_g <- otu_t_g[, -1]
    metadata3 <- otu_t_g[-c(1:ncol(otu_t))]
    rownames(metadata3) <- NULL
    colnames(metadata3)[1] <- "sample"

    #
    otu_t2 <- otu_t %>%
      dplyr::group_by(metadata3[["sample"]]) %>%
      dplyr::summarize_at(ggplot2::vars(-dplyr::group_cols()), parallel_method)

    #
    otu_t2 <- as.data.frame(otu_t2)
    row.names(otu_t2) <- otu_t2[, 1]
    otu_t2 <- otu_t2[, -1]

    #
    otu3 <- as.data.frame(t(otu_t2))
    metadata4 <- unique(as.data.frame(metadata3[["sample"]]))
    colnames(metadata4)[1] <- "sample"
    rownames(metadata4) <- NULL

    #
    matched_indices <- match(metadata4[["sample"]], metadata3[["sample"]])
    unique_group <- metadata3[[group]][matched_indices]
    metadata5 <- data.frame(metadata4, unique_group)
    colnames(metadata5) <- c("sample", "group")
    if(parallel_method == "mean") {
      otu3 = round(x = otu3, digits = digits)
    }


    ##
  } else if(parallel_method %in% c("none")) {
    otu3 <- otu2
    metadata5 <- metadata2
    colnames(metadata5) <- c("sample", "parallel", "group")
  }


  ##
  if(isTRUE(metadata_out)) {
    result = NULL
    result = list(otu = otu3, metadata = metadata5)
  } else {

    #
    return(otu3)
  }
}
