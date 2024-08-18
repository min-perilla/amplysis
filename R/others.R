#' @title Set Working Directory
#'
#' @description
#' By calling the setwd() function, you can quickly set the location of the
#' current R file as the working directory.
#'
#' @return character
#' @export
#'
#' @examples \dontrun{set_wd()}
#'
#' @importFrom rstudioapi getActiveDocumentContext
#'

set_wd <- function() {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  cat("Set Working Directory: \n\"", getwd(), "\"\n", sep = "")
}
