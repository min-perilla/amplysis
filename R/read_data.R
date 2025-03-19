#' @title Read Data
#'
#' @description Can automatically recognize and import CSV, tsv, txt, xls, xlsx, biom, fasta, fastq files.
#'
#' @param file (character) File path, such as "file = ./table. csv".
#' @param col_names (logical) Set the first row as the column header, defaulting to TURE, for example "col_names = TRUE".
#' @param row_names (logical) Use the first column # OTU ID as the row name, defaulting to TURE, for example "otu_rownames = True".
#' @param delim (character) Custom segmentation symbols are only valid for importing txt files, such as "delim = "\\t".
#' @param skip (integer) Ignoring the first n rows of metadata, such as "skip = 1".
#' @param sheet (integer) Only valid for xls or xlsx files, select which sheet to read, such as "sheet = 2".
#' @param range (character) The range of tables to be read, such as "sheet_name!B2:G14".
#' @param fast_to_csv (logical) Whether to convert FAST files to CSV files. The default value is FALSE.
#' @param f_to_c_fileName (character) Customize the CSV file name? If yes, please enter; otherwise, use the default value "fast". For example, if you type "rep_seq", the file name will be saved as "rep_seq.csv".
#'
#' @return data.frame
#' @export
#'
#' @examples
#' \dontrun{read_data()}
#'
#' @importFrom readr read_csv read_tsv read_delim write_excel_csv
#' @importFrom readxl read_excel
#' @importFrom phyloseq import_biom otu_table tax_table
#' @importFrom Biostrings readDNAStringSet
#' @importFrom ape read.tree
#'

# Please use tools:: showNonASCIIfile(file. R) to check for the presence of non ASCII characters.
# tools::showNonASCIIfile(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "read_data.R"))
read_data <- function(
    file = NULL,             # file path
    col_names = TRUE,        # Whether to convert the first row into column names, default TRUE.
    row_names = FALSE,       # Whether to convert the first column into column names, default TRUE.
    delim = NULL,            # Specify the delimiter; when provided, it will call the `readr::read_delim()` function.
    skip = 0,                # Ignore the first n rows when reading the data

    # The following parameters only take effect when reading Excel files.
    sheet = NULL,  # Specify which sheet to read, for example: sheet = 2.
    range = NULL,  # Read a specific range, for example: range = sheet name!B2:G14.

    # The following parameters only apply to reading 'fast' series files.
    fast_to_csv = FALSE,       # Whether to convert the 'fast' file to a 'CSV' file.
    f_to_c_fileName = "fast"   # Customize the CSV file name without entering the extension. For example, entering "rep_seq" will save the file as "rep_seq.csv".
) {

  # Initialization results
  result <- NA

  # Marking for the fast series files
  flag_fast = -1  # Fasta is indicated by 1, and Fastq is indicated by 2.

  # check if the file name is empty.
  if(is.null(file)) {
    stop("File name is empty!")
  } else {

    # Retrieve the file type and convert it to lowerca using the base package's 'tolower()' function.
    fileType <- tolower(tools::file_ext(file))

    # Read the file using different methods based on its file type.
    # 1. csv file
    if(fileType == "csv") {
      # call the 'read_csv()' function from the 'readr' package.
      result <- read_csv(file = file, col_names = col_names, skip = skip)
    }


    # 2. tsv file
    else if(fileType == "tsv") {
      # call the 'read_tsv()' function from the 'readr' package.
      result <- read_tsv(file = file, col_names = col_names, skip = skip)
    }


    # 3. txt file
    # Require user specification of the delimiter, defaulting to `\t`, for a TXT file.
    else if(fileType == "txt") {
      # call the `read_delim()` function from the `readr` package.
      result <- read_delim(file = file, delim = delim, col_names = col_names, skip = skip)
    }


    # 4. excel file
    else if(fileType == "xls" | fileType == "xlsx") {
      # call the `read_excel()` function from the `readxl` package.
      result <- read_excel(path = file, sheet = sheet, range = range, skip = skip)
    }


    # 5. biom
    else if(fileType == "biom") {
      # call the `import_biom()` function from the `phyloseq` package.
      biom_data <- import_biom(BIOMfilename = file)

      # otu table
      otu_df <- as.data.frame(otu_table(biom_data))    # Extract the OTU table.
      # colnames(otu_df)[1] <- "#OTU ID"               # Change the column name of the first column to "#OTU ID".


      # tax table
      tax_df <- as.data.frame(tax_table(biom_data))  # Extract the tax table.
      # Rename column names.
      colnames(tax_df) <- c("domain", "phylum", "class", "order", "family", "genus", "species")

      # Left join `otu_df` and `tax_df` with `all = FALSE` to discard unmatched rows, and `sort = FALSE` to avoid re-sorting.
      result <- merge(otu_df, tax_df, by.x = "row.names", by.y = "row.names", all = FALSE, sort = FALSE)

      # Rename the first column as `#OTU ID`.
      colnames(result)[1] <- "#OTU ID"
    }


    # 6. fasta (rep_seq)
    # Read a FASTA format file and convert it to a data frame and CSV file
    # .fna: FASTA nucleotide file
    # .faa: FASTA amino acid file
    # .ffn: FASTA nucleotide coding region file
    else if(fileType %in% c("fasta", "fa", "fna", "faa", "ffn")) {
      # Call the `readDNAStringSet` function from the R package Biostrings to read a FASTA file.
      result <- readDNAStringSet(filepath = file,
                                 format = "fasta",
                                 skip = skip
      )

      # Output summary information.
      # print(result@ranges)
      print(result)

      result <- as.data.frame(result)  # Convert it to the data.frame format.
      colnames(result) <- c("rep_seq") # Rename the first column.


      # If `row_names = FALSE`, set the row names as the first column and initialize row names.
      if(isFALSE(row_names)){
        cat("row_names is FALSE")
        # Add row names as the first column and set the column name to "#OTU ID".
        result <- cbind("#OTU ID" = rownames(result), result)
        rownames(result) <- NULL  # Initialize row names.
      }


      # Should the FASTA file be converted to a CSV file?
      # yes
      if(isTRUE(fast_to_csv)) {
        cat("fast_to_csv is TRUE")
        result2 <- result

        # If `row_names` is set to TRUE, the following operations will add row names as the first column.
        if(isTRUE(row_names)) {
          # Add row names as the first column and set the column name to "#OTU ID".
          result2 <- cbind("#OTU ID" = rownames(result2), result2)
          rownames(result2) <- NULL  # Initialize row names.
        }


        # Generate a file name.
        file_name_csv <- paste0(f_to_c_fileName, ".csv")

        # Call the `write_excel_csv` function from the readr package to save the data as a CSV file.
        readr::write_excel_csv(x = result2, file = file_name_csv)

        cat("\033[0;32m", "\nThe file \"", file_name_csv, "\" has been saved to \n",
            getwd(), "/", file_name_csv, "\033[0m\n", sep = "")
        # cat("\033[1;32m", "", "\033[0m\n", sep = "")  # Green text.
      }

      # Set the identifier to 1 to indicate that the input is a FASTA file.
      flag_fast = 1
    }


    # 7. fastq (rep_seq)
    # Read a FASTQ format file and convert it into a data frame and CSV file.
    else if(fileType %in% c("fastq", "fq")) {
      # Call the `readDNAStringSet` function from the R package Biostrings to read a FASTA file.
      result <- readDNAStringSet(filepath = file,
                                 format = "fastq",
                                 skip = skip
      )

      # Output summary information.
      # print(result@ranges)
      print(result)

      result <- as.data.frame(result)  # Convert it to the data.frame format.
      colnames(result) <- c("rep_seq") # Rename the first column.


      # If `row_names = FALSE`, set the row names as the first column and initialize row names.
      if(isFALSE(row_names)){
        cat("row_names is FALSE")
        # Add row names as the first column and set the column name to "#OTU ID".
        result <- cbind("#OTU ID" = rownames(result), result)
        rownames(result) <- NULL  # Initialize row names.
      }


      # Should the FASTA file be converted to a CSV file?
      # yes
      if(isTRUE(fast_to_csv)) {
        cat("fast_to_csv is TRUE")
        result2 <- result

        # If `row_names` is set to TRUE, the following operations will add row names as the first column.
        if(isTRUE(row_names)) {
          # Add row names as the first column and set the column name to "#OTU ID".
          result2 <- cbind("#OTU ID" = rownames(result2), result2)
          rownames(result2) <- NULL  # Initialize row names.
        }


        # Generate a file name.
        file_name_csv <- paste0(f_to_c_fileName, ".csv")

        # Call the `write_excel_csv` function from the readr package to save the data as a CSV file.
        readr::write_excel_csv(x = result2, file = file_name_csv)

        cat("\033[0;32m", "\nThe file \"", file_name_csv, "\" has been saved to \n",
            getwd(), "/", file_name_csv, "\033[0m\n", sep = "")
        # cat("\033[1;32m", "", "\033[0m\n", sep = "")  # Green text.

      }

      # Set the identifier to 2 to indicate that the input is a FASTQ file.
      flag_fast = 2
    }


    # 8. Read the tree file.
    # nwk format: Newick
    # nhx format: New Hampshire eXtended
    # tree format: tree
    else if(fileType %in% c("nwk", "nhx", "tree")) {
      tree <- ape::read.tree(file)

      cat("\033[0;32m",  # Green text.
          "The tree file has been successfully read!\n",
          "File type: ", fileType, "\n",
          "\033[0m\n",   # Restore the font color.
          sep = "")

      return(tree)
    }



    else {
      stop("The file type is incorrect!")
    }



    # Set the # OTU ID in the first column to the row name and remove the first column
    # No need to convert Fasta and Fastq format files
    if(isTRUE(row_names) && flag_fast == -1) {
      result <- as.data.frame(result)
      rownames(result) <- result[, 1]       # Set the first column as the row name
      result <- result[, -1, drop = FALSE]  # Remove the first column, parameter drop=FALSE means that when the data is one-dimensional, it will not be converted to a vector
    }


    # Output result information
    cat("\033[0;32m", "\nSuccess: Data imported!\n",
        "File type: ", fileType, "\n",
        "Otu rownames: ", row_names, "\n",
        "Skip: ", skip, "\033[0m\n", sep = "")

    return(result)
  }
}

