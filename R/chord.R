#' @title Chord Diagram
#'
#' @description
#' Process the feature table into the necessary data format for drawing a chord
#' diagram.
#'
#' @param otu otu table
#' @param metadata metadata table
#' @param tax tax table
#'
#' @param id_col (integer) The OTU_ID column is in which column,
#' defaulting to 0 means there is no OTU_ID column, and the data is already
#' numeric.
#' @param tax_cla (character) Taxonomic level. Only column names from the tax
#' table can be entered, for example tax_cla = 'genus'.
#'
#' @param group (Required, character) Grouping information. please enter the column name of
#' the grouping information in the metadata table.
#' @param parallel_method  (character) Sample processing methods for the same group:
#' average, sum, median, none.
#' @param row_n (integer) Preserve the top N taxa (including the Nth) based on
#' abundance, while merging taxa with lower abundance into "others".
#' @param digits (integer) When the parallel sample method is set to "mean," the
#' number of decimal places will default to 0.
#'
#' @return a data table
#' @export
#'
#' @examples
#' \dontrun{
#' test = chord(otu = otu, metadata = metadata, tax = tax, id_col = 1,
#' tax_cla = "genus", group = "group", parallel_method = "mean", row_n = 8,
#' digits = 0)}
#'
#' @importFrom dplyr across arrange desc group_by mutate row_number summarize
#' summarise_all ungroup where
#'
chord = function(otu, metadata, tax, id_col = 1, tax_cla = "genus",
         group = "group", parallel_method = "mean", row_n = 8, digits = 0)
{
  ##
  otu_matadata = amplysis::parallel(
    otu = otu, metadata = metadata, id_col = id_col,
    group = group, parallel_method = parallel_method,
    digits = digits, metadata_out = T)

  otu2 = otu_matadata[["otu"]]
  metadata2 = otu_matadata[["metadata"]]

  otu_tax = base::merge(x = tax[, c(1, which(colnames(tax) == tax_cla))], y = otu2,
                        by.x = 1, by.y = "row.names", sort = F)
  otu_tax = otu_tax[, -1]

  otu_tax2 <- otu_tax %>%
    dplyr::group_by(otu_tax[[tax_cla]]) %>%
    dplyr::summarize(dplyr::across(dplyr::where(is.numeric), sum, na.rm = TRUE)) %>%
    dplyr::arrange(dplyr::desc(rowSums(dplyr::across(dplyr::where(is.numeric))))) %>%
    dplyr::ungroup()

  colnames(otu_tax2)[1] = tax_cla
  otu_tax3 <- otu_tax2 %>%
    dplyr::mutate(ifelse(dplyr::row_number() <= row_n, otu_tax2[[tax_cla]], "others"))

  otu_tax3[, 1] <- otu_tax3[,ncol(otu_tax3)]
  otu_tax3 <- otu_tax3[, -ncol(otu_tax3)]


  ##
  otu_tax4 <- otu_tax3 %>%
    dplyr::group_by(dplyr::select(otu_tax3, tidyr::all_of(tax_cla))) %>%
    dplyr::summarise_all(sum) %>%
    dplyr::ungroup()

  row_numbers_others <- which(otu_tax4[, 1] == "others")
  otu_tax4 <- rbind(otu_tax4[-row_numbers_others, ], otu_tax4[row_numbers_others, ])


  ##
  otu_tax4 = as.data.frame(otu_tax4)
  rownames(otu_tax4) = otu_tax4[, 1]
  otu_tax4 = otu_tax4[, -1]


  ##
  otu_tax4 = as.matrix(otu_tax4)
  df = data.frame(from = rep(rownames(otu_tax4), times = ncol(otu_tax4)),
                  to = rep(colnames(otu_tax4), each = nrow(otu_tax4)),
                  value = as.vector(otu_tax4),
                  stringsAsFactors = FALSE)


  ##
  cat("\033[32m--- Please use the `chord_plot()` function for visualization. ---\n\033[0m")

  # 返回数据
  return(df)
}
