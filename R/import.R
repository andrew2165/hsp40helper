#' AnalysisData S4 Class
#'
#' An S4 class to store file path and loaded CSV data for the analysis pipeline.
#'
#' @slot filepath Character. The file path to the CSV data file.
#' @slot dataframe Data.frame. The loaded data.

#' @export
setClass(
  "AnalysisData",
  slots = list(
    filepath = "character",
    dataframe = "data.frame"
  )
)

#' @describeIn AnalysisData S4 initializer that loads data from a CSV file path.
#' @import dplyr
setMethod(
  f = "initialize",
  signature = "AnalysisData",
  definition = function(.Object, filepath = NA_character_, ...) {
    if (!is.na(filepath) && file.exists(filepath)) {
      df <- read.csv(filepath, ...)
      .Object@filepath <- filepath
      .Object@dataframe <- df
    } else {
      stop("A valid file path must be provided and the file must exist.")
    }
    .Object
  }
)

#' Create a new AnalysisData object
#'
#' @param filepath Path to the CSV data file.
#' @param ... Additional arguments passed to read.csv.
#' @return A new AnalysisData object with data loaded from the CSV file.
#' @export
AnalysisData <- function(filepath, ...) {
  new("AnalysisData", filepath = filepath, ...)
}
