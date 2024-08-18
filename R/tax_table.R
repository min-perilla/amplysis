#' @title Data Separation
#'
#' @description Call the `separate_wider_delim()` function from the tidyr package
#' to separate the data in the taxonomy column of the data table.
#'
#' @param tax The data table containing the "tax" column.
#' @param index (integer) The column number of the taxonomy column in the data table.
#' @param delim (character) The delimiter. Default is ";".
#' @param names (character) Names for each column after data separation.
#' Ensure the number of names matches the number of columns after separation,
#' such as 'names = c("domain", "phylum", "class", "order", "family", "genus", "species")'.
#'
#' @return data.frame
#' @export
#'
#' @examples
#' \dontrun{tax_separate(tax = data, index = 2)}
#' \dontrun{tax_separate(tax = data, index = 2, delim = ";",
#' names = c("domain", "phylum", "class", "order", "family", "genus", "species"))}
#'
#' @importFrom tidyr separate_wider_delim
#'
# Please use tools:: showNonASCIIfile(file. R) to check for the presence of non ASCII characters.
# tools::showNonASCIIfile(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "tax_table.R"))
tax_separate <- function(tax, index, delim = NULL,
         names = c("domain", "phylum", "class", "order", "family", "genus", "species"))
{
  tax_rownames <- rownames(tax)

  tax_separated <- tidyr::separate_wider_delim(data = tax,
    cols = tidyr::all_of(index), delim = delim, names = names,
    too_few = "align_start", too_many = "merge")

  tax_separated <- as.data.frame(tax_separated)
  rownames(tax_separated) <- tax_rownames

  cat("\nData split successfully!\n", "cols: ", index, "\n", "delim: `", delim,
      "`\n", sep = "")

  return(tax_separated)
}
################################################################################










#' @title Remove The prefix From Tax Information
#' @description
#' The function can remove a fixed-length prefix from taxonomic information, such as p__Firmicutes -> Firmicutes.
#'
#' @param tax (data) A table containing classification information.
#' @param index (integer or character) Which columns contain the classification information to be repaired,
#' for example: index = 3, index = c(2:8), index = c("domain", "phylum", "class", "order", "family", "genus", "species").
#' @param position (integer) The length of the prefix.
#'
#' @return data.frame
#' @export
#'
#' @examples \dontrun{tax_dePrefix <- tax_trim_prefix(tax = tax, index = c(2:8), position = 3)}
#'
tax_trim_prefix <- function(tax, index, position)
{
  table2 <- tax

  pattern = paste0("^.{", position, "}")

  remove_pattern <- function(x, pattern) {
    base::gsub(pattern = pattern, replacement = "", x = x)}

  result <- base::apply(
    X = table2[, index, drop = F], MARGIN = 2,
    FUN = function(col)
    { base::sapply(X = col, FUN = remove_pattern, pattern = pattern) }
  )

  result <- as.data.frame(result)
  table2[, index] <- result

  return(table2)
}
################################################################################










#' @title Taxonomic Information Repair
#'
#' @description
#' In taxonomic information, terms such as \"unclassified\", \"unknown\", \"uncultured\",
#' \"norank\", \"unidentified\", \"Unknown_Species\", \"Unknown_Genus\", \"Unknown_Family\",
#' \"Unknown_Order\", \"Unknown_Class\", and \"Unknown_Phylum\" may exist.
#' When plotting a species composition chart, simply merging these categories may result in the loss of some taxonomic information.
#'
#' @description
#' For instance, one entry may be categorized as Firmicutes at the phylum level and as
#' \"uncultured\" at the genus level, while another entry may be categorized as Armatimonadetes at
#' the phylum level and also as \"uncultured\" at the genus level. Merging them at
#' the genus level would result in both entries being classified as \"uncultured\",
#' leading to the loss of some taxonomic information. The correct approach is to refine
#' and repair the taxonomic information. In this case, the entries at the genus level should be
#' specified as \"uncultured_Firmicutes\" and \"uncultured_Armatimonadetes\", respectively,
#' to preserve the detailed taxonomic information.
#'
#' @param tax (data) A table containing classification information.
#' @param column_to_check (integer or character) The columns to be repaired, such as column_to_check = 8, column_to_check = c(4:8),
#' If column names are known, you can also input column_to_check = "column_name"
#' @param column_to_add (integer or character) If column_to_check is unclassified, append the column number as a suffix.
#'
#' @return data.frame
#' @export
#'
#' @examples \dontrun{tax_names_repair(tax = tax, column_to_check = 7, column_to_add = 3)}
#' @examples \dontrun{tax_names_repair(tax = tax, column_to_check = c(4:8), column_to_add = 3)}
#' @examples \dontrun{tax_names_repair(tax = tax, column_to_check = "genus",
#' column_to_add = "phylum")}
#' @examples \dontrun{tax_names_repair(tax = tax, column_to_check = c("genus", "species"),
#' column_to_add = "phylum")}
#' @examples \dontrun{tax_names_repair(tax = tax, column_to_check = c("genus", "species"),
#'  column_to_add = 3)}
#'
#' @importFrom dplyr mutate_at
#' @importFrom dplyr vars
#' @importFrom dplyr all_of
#'

# Please use tools:: showNonASCIIfile(file. R) to check for the presence of non ASCII characters.
# tools::showNonASCIIfile(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "tax_table.R"))

tax_names_repair <- function(tax,              # A table containing classification information.
                             column_to_check,  # The columns to be repaired, such as column_to_check = 8, column_to_check = c(4:8),
                                               # If column names are known, you can also input column_to_check = "column_name"
                             column_to_add     # If column_to_check is unclassified, append the column number as a suffix.
) {
  # Format check: column_to_check
  if (is.numeric(column_to_check) && any(floor(column_to_check) != column_to_check)) {
    stop("column_to_check must be an integer!")
  }

  # Format check: column_to_add
  if (!(is.numeric(column_to_add) && length(column_to_add) == 1) &&
      !(is.character(column_to_add) && length(column_to_add) == 1)) {
    stop("column_to_add must be a single integer or a single string!")
  }

  # Format check: Ensure that column_to_check and column_to_add are not the same.
  if (!is.null(column_to_check) && !is.null(column_to_add)) {
    if (identical(column_to_check, column_to_add)) {
      stop("column_to_check and column_to_add cannot be the same!")
    }
  }

  # Define a pattern list.
  patterns <- base::tolower(c("unclassified", "unknown", "uncultured", "norank", "unidentified",
                              "Unknown_Species", "Unknown_Genus", "Unknown_Family",
                              "Unknown_Order", "Unknown_Class", "Unknown_Phylum"))

  # Apply the pattern and add suffix.
  tax <- dplyr::mutate_at(tax, dplyr::vars(dplyr::all_of(column_to_check)), function(x) {
    ifelse(base::tolower(x) %in% patterns, paste0(x, '_', tax[[column_to_add]]), x)
  })

  return(tax)
}

