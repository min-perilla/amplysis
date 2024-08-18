#' @title Save as Rdata
#' @description
#' The use of the `save()` function must ensure that the directory exists.
#' In contrast, the `save_Rdata()` function can automatically check whether the
#' directory exists (and create it if it does not), making it easier to save
#' data as RData.
#'
#' @param ... (Parameter) Parameter
#' @param file (character) File path
#'
#' @return Rdata
#' @export
#'
#' @examples
#' \dontrun{
#' otu <- read_data("otu.csv")
#' save_Rdata(otu, "data/otu.Rdata")
#' }
save_Rdata <- function(..., file) {
  ##
  folders <- unlist(strsplit(file, "/"))
  filename = folders[length(folders)]
  folders <- folders[-length(folders)]

  #
  path <- getwd()

  #
  for (folder in folders) {
    path <- file.path(path, folder)

    #
    if (!dir.exists(path)) {
      #
      dir.create(path)
      cat("Created folder: ", path, "\n")
    } else {
      cat("Folder already exists:", path, "\n")
    }
  }


  ##
  save(... = ..., file = file)

  return("success")
}
