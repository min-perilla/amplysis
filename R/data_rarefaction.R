#' @title Data Rarefaction
#'
#' @description
#' Use the "rrarefy" function from the "vegan" package or the "rarefy_even_depth" function from the "phyloseq" package for data rarefaction.
#'
#' @param file (data.frame) OTU data frame
#' @param id_col (integer) The OTU_ID column is in which column, defaulting to 0 means there is no OTU_ID column, and the data is already numeric.
#' @param method (character) The data rarefaction can be performed using the "vegan" package or the "phyloseq" package.
#' @param seed (integer) (Only applicable when using the "phyloseq" package for data rarefaction) Random seed.
#' @param replace (logical) Whether to sample with replacement (TRUE) or without replacement (FALSE).
#' (1) Sampling with replacement is faster and more memory-efficient;
#' (2) Sampling with replacement means that the read counts for some OTUs could
#' be greater than the original count values.
#' @param trimOTUs (logical) (Only applicable when using the "phyloseq" package for data rarefaction) Whether to remove OTUs with zero counts for each sample.
#' @param tax_table (tax table)
#' @param write_file (logical) Whether to save the rarefied file.
#' @param file_name (character) Customize the saved file name.
#'
#' @return Return an OTU table that has been rarefied,
#' or return a list containing a rarefied OTU table and an aligned taxonomy table after rarefaction.
#' @export
#'
#' @examples
#' \dontrun{data_rarefy(file = otu_table)}
#' \dontrun{data_rarefy(file = otu_table, id_col = 0, method = "vegan")}
#' \dontrun{data_rarefy(file = otu_table, id_col = 0, method = "phyloseq", seed = 1231)}
#' \dontrun{data_rarefy(file = otu, id_col = 1, method = "phyloseq",
#' trimOTUs = T, tax_table = tax, write_file = T)}
#'
#' @importFrom vegan rrarefy
#' @importFrom phyloseq otu_table rarefy_even_depth
#' @importFrom readr write_excel_csv
#'

# Please use tools:: showNonASCIIfile(file. R) to check for the presence of non ASCII characters.
# tools::showNonASCIIfile(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "data_rarefaction.R"))
data_rarefy <- function(file, id_col = 0, method = "vegan", seed = 123,
         replace = F, trimOTUs = T, tax_table = NULL, write_file = TRUE,
         file_name = "table")
{
  if(id_col > 0) {
    file <- as.data.frame(file)
    rownames(file) <- file[, id_col]
    file <- file[, -id_col, drop = FALSE]
  } else {}

  is_numeric <- sapply(
    file, function(col) {
      all(grepl("^\\d+$", as.character(col)))
    }
  )

  is_numeric <- as.data.frame(is_numeric)

  if (all(is_numeric)) {

  } else {
    non_numeric_cols <- which(!is_numeric)
    real_col <- non_numeric_cols

    if(isTRUE(!id_col == 0)) {
      real_col <- non_numeric_cols + 1
    }

    message(paste("Columns ", paste(real_col, collapse = ", "), " are not numeric.", sep = ""))

    cat("Previewing first 3 rows:\n", sep = "")
    print(file[1:3, non_numeric_cols])

    stop("Please review the data and try again!")
  }

  ##
  maxColSums <- max(colSums(file))
  minColSums <- min(colSums(file))

  method <- tolower(method)

  if(method == "vegan") {
    otu_rarefy <- t(vegan::rrarefy(t(file), minColSums))
    otu_rarefy <- as.data.frame(otu_rarefy)
  } else if(method == "phyloseq") {
    set.seed(seed = seed)
    sue = 1231
    otu_table <- phyloseq::otu_table(object = file, taxa_are_rows = T)
    otu_table <- phyloseq::rarefy_even_depth(
      physeq = otu_table,
      replace = replace,
      trimOTUs = trimOTUs,
      verbose = T)
    otu_rarefy <- as.data.frame(otu_table@.Data)
  } else {
    stop("\"method\" parameter error, please enter \"method = vegan\"",
         " or \"method = pholoseq\" and try again.\n", sep = "")
  }

  cat("\033[0;32m", "Succeed!\n",
      "Method: ", method, "\n",
      "Colsums_Maximum: ", maxColSums, "\n",
      "Colsums_Minimum: ", minColSums, "\n",
      "All colsums -> ", minColSums, "\n",
      "Data type: ", class(otu_rarefy), "\033[0m\n", sep = "")
  if(method == "phyloseq") {
    cat("\033[0;32m", "Seed: ", seed, "\033[0m\n", sep = "")
  }

  otu_rarefy <- cbind(OTU_ID = rownames(otu_rarefy), otu_rarefy)
  colnames(otu_rarefy) <- c("#OTU ID", colnames(otu_rarefy)[-1])
  rownames(otu_rarefy) <- NULL

  merged_data <- NULL

  if(method == base::tolower("phyloseq") && isTRUE(trimOTUs)) {
    if(is.null(tax_table)) {

    } else {
      if(id_col > 0) {
        merged_data <- align_otu_tax(
          tax_table = tax_table, otu_table = otu_rarefy, by.x = 1, by.y = 1,
          all.x = F, all.y = T, sort = F)

      } else {
        merged_data <- align_otu_tax(
          tax_table = tax_table, otu_table = otu_rarefy, by.x = "row.names",
          by.y = 1, all.x = F, all.y = T, sort = F)
      }

      if(isTRUE(write_file)) {
        file_name_otu = paste0(file_name, "_rarefy_p.csv")
        file_name_tax = "tax_alignment.csv"

        otu_rarefy2 <- merged_data[[1]]
        tax_align2 <- merged_data[[2]]

        if(id_col <= 0) {
          otu_rarefy3 <- cbind(OTU_ID = rownames(otu_rarefy2), otu_rarefy2)
          tax_align3 <- cbind(OTU_ID = rownames(tax_align2), tax_align2)

          colnames(otu_rarefy3) <- c("#OTU ID", colnames(otu_rarefy3)[-1])
          colnames(tax_align3) <- c("#OTU ID", colnames(tax_align3)[-1])

          rownames(otu_rarefy3) <- NULL
          rownames(tax_align3) <- NULL

          readr::write_excel_csv(x = otu_rarefy3, file = file_name_otu)
          readr::write_excel_csv(x = tax_align3, file = file_name_tax)

          cat("\033[0;32m", "\nThe file \"", file_name_otu, " ,", file_name_tax, "\" has been saved to \n",
              getwd(), "/", file_name_otu, "\n",
              getwd(), "/", file_name_tax, "\033[0m\n",
              sep = "")

        } else {
          readr::write_excel_csv(x = otu_rarefy2, file = file_name_otu)
          readr::write_excel_csv(x = tax_align2, file = file_name_tax)

          cat("\033[0;32m", "\nThe file \"", file_name_otu, " ,", file_name_tax, "\" has been saved to \n",
              getwd(), "/", file_name_otu, "\n",
              getwd(), "/", file_name_tax, "\033[0m\n",
              sep = "")
        }
      }

      return(merged_data)
    }
  }




  if(isTRUE(write_file)) {
    if(method == "vegan") {
      file_name2 = paste0(file_name, "_rarefy_v.csv")

    } else if(method == "phyloseq") {
      file_name2 = paste0(file_name, "_rarefy_p.csv")
    }

    readr::write_excel_csv(x = otu_rarefy, file = file_name2)

    cat("\033[0;32m", "\nThe file \"", file_name2, "\" has been saved to \n",
        getwd(), "/", file_name2, "\033[0m\n", sep = "")
  }

  if(id_col <= 0) {
    rownames(otu_rarefy) <- otu_rarefy[, 1]
    otu_rarefy <- otu_rarefy[, -1]
    return(otu_rarefy)

  } else {
    return(otu_rarefy)
  }
}



################################################################################
# When using the rarefaction method from the phyloseq package, with trimOTUs = TRUE, align the tax and otu tables using the following method.
#' @title Align OTU and TAX.
#'
#' @description
#' When using the rarefaction method of the phyloseq package with the parameter
#' trimOTUs = TRUE, some OTUs in the OTU table will be removed. In this case,
#' the number of rows in the OTU table and the TAX table will not be consistent.
#' The align_otu_tax() function can be used to align these two tables.
#'
#' @param tax_table tax table
#' @param otu_table otu table
#' @param by.x (integer or character) Used to specify the column name in the tax data frame,
#' defaulting to 1, indicating the column name of the first column. Alternatively,
#' you can directly enter the column name, such as by.x = "OTU ID".
#' @param by.y (integer or character) Used to specify the column name in the otu data frame,
#' defaulting to 1, indicating the column name of the first column. Alternatively,
#' you can directly enter the column name, such as by.x = "OTU ID".
#' @param all.x (logical) Should all rows of the tax data frame be retained? Default is FALSE.
#' @param all.y (logical) Should all rows of the otu data frame be retained? Default is FALSE.
#' @param sort (logical) Should reordering be done? Default is no.
#'
#' @return Return a list containing the aligned OTU and TAX tables.
#' @export
#'
#' @examples \dontrun{align_otu_tax(tax_table = tax, otu_table = otu_rare_p,
#' by.x = 1, by.y = 1, all.x = F, all.y = T, sort = F)}
#'
align_otu_tax <- function(tax_table, otu_table, by.x = 1, by.y = 1, all.x = F,
         all.y = F, sort = F)
{
  result_list = NULL
  len_tax <- length(tax_table)

  merged_data <- base::merge(x = tax_table, y = otu_table,
                             by.x = by.x, by.y = by.y,
                             all.x = all.x, all.y = all.y, sort = sort)

  if(by.x == "row.names") {
    len_tax <- len_tax + 1
  }

  tax_table2 <- merged_data[, c(1:len_tax)]
  otu_table2 <- merged_data[, -c(2:len_tax)]

  if(by.x == "row.names") {
    rownames(otu_table2) <- otu_table2[, 1]
    rownames(tax_table2) <- tax_table2[, 1]

    otu_table2 <- otu_table2[, -1]
    tax_table2 <- tax_table2[, -1]
  }

  result_list = list(otu = otu_table2, tax = tax_table2)

  return(result_list)
}

