#' ProteinAnalysisData S4 Class
#'
#' An S4 class to store file path and loaded CSV data for the analysis pipeline.
#'
#' @slot filepath Character. The file path to the CSV data file.
#' @slot dataframe Data.frame. The loaded data.

#' @export
setClass(
  "ProteinAnalysisData",
  slots = list(
    filepath = "character",
    dataframe = "data.frame",
    vehicleSamples = "character",
    treatmentSamples = "character",
    baitID = "character"
  )
)

#' @describeIn ProteinAnalysisData S4 initializer that loads data from a CSV file path.
#' @import dplyr
#' @import tools
setMethod(
  f = "initialize",
  signature = "ProteinAnalysisData",
  definition = function(.Object, filepath = NA_character_, ...) {
    if (!is.na(filepath) && file.exists(filepath) && (tools::file_ext(filepath) == "tsv")) {
      df <- read.csv(filepath, ...)
      .Object@filepath <- filepath
      .Object@dataframe <- df
    } else {
      stop("A valid file path must be provided and the file must exist.")
    }
    .Object
  }
)


#' Create a new ProteinAnalysisData object
#'
#' @param filepath Path to the CSV data file.
#' @param ... Additional arguments passed to read.csv.
#' @return A new ProteinAnalysisData object with data loaded from the CSV file.
#' @export
read.fragpipe <- function(filepath, ...) {
  new(Class = "ProteinAnalysisData", filepath = filepath, ...)
}


#' Set Samples
#'
#' @param object A ProteinAnalysisData object.
#' @param vehicleSamples Character vector of column names for vehicle samples.
#' @param treatmentSamples Character vector of column names for treatment samples.
#' @return The updated ProteinAnalysisData object.
#' @import dplyr
#' @export
setGeneric(
  "setSamples",
  function(object, vehicleSamples, treatmentSamples, ...) standardGeneric("setSamples")
)
setMethod(
  f = "setSamples",
  signature = "ProteinAnalysisData",
  definition = function(object, vehicleSamples, treatmentSamples, ...) {
    # Store the vectors as slots
    object@vehicleSamples <- vehicleSamples
    object@treatmentSamples <- treatmentSamples
    return(object)
  }
)


# TODO: set a standard set of protein IDs to be removed
# e.g., Immunoglobulins, non-human proteins, etc.
#' Remove unwanted ProteinIDs in ProteinAnalysisData
#'
#' @param object A ProteinAnalysisData object.
#' @param idColumn The Protein ID column name (character) to normalize.
#' @param cleanUpIDs The Protein IDs to be removed (list).
#' @return The updated ProteinAnalysisData object.
#' @import dplyr
#' @export
setGeneric(
  "cleanData",
  function(object, idColumn, cleanUpIDs, ...) standardGeneric("cleanData")
)
setMethod(
  f = "cleanData",
  signature = "ProteinAnalysisData",
  definition = function(object, idColumn, cleanUpIDs, ...) {
    df <- object@dataframe
    if (idColumn %in% names(df)) {
      # df[[column]] <- (df[[column]] - min(df[[column]], na.rm = TRUE)) /
      #                 (max(df[[column]], na.rm = TRUE) - min(df[[column]], na.rm = TRUE))
      # object@dataframe <- df
      df <- dplyr::filter(df, !(idColumn %in% cleanUpIDs))
      object@dataframe <- df
    } else {
      warning("Protein ID Column not found.")
    }
    return(object)
  }
)

# TODO: write a method to normalize to bait (normally DNAJB8)
#' Normalize to Bait
#'
#' @param object A ProteinAnalysisData object.
#' @param idColumn The Protein ID column name (character) to normalize.
#' @param baitID The UniProt Protein ID of the Bait used (character).
#' @param removeBait Remove the Bait entry from the object (boolean).
#' @return The updated ProteinAnalysisData object.
#' @import dplyr
#' @export
setGeneric(
  "baitNormalize",
  function(object, idColumn, cleanUpIDs, removeBait = TRUE, ...) standardGeneric("baitNormalize")
)
setMethod(
  f = "baitNormalize",
  signature = "ProteinAnalysisData",
  definition = function(object, idColumn, baitID, removeBait = TRUE, ...) {
    df <- object@dataframe
    if (baitID %in% idColumn) {
      # df[[column]] <- (df[[column]] - min(df[[column]], na.rm = TRUE)) /
      #                 (max(df[[column]], na.rm = TRUE) - min(df[[column]], na.rm = TRUE))
      # object@dataframe <- df
      # df = dplyr::filter(df, !(idColumn %in% cleanUpIDs))
      # object@dataframe = df

      # TODO: write this part of the function
    } else {
      warning("Bait ProteinID not found.")
    }
    return(object)
  }
)

# TODO: write a method to calculate fold change
#' Calculate Fold Change
#'
#' @param object A ProteinAnalysisData object.
#' @param idColumn The Protein ID column name (character).
#' @param vehicleSamples The columns containing the normalized vehicle intensities (list).
#' @param treatmentSamples The columns containing the normalize treatment intensities (list).
#' @return The updated ProteinAnalysisData object.
#' @import dplyr
#' @export
setGeneric(
  "calculateFoldChange",
  function(object, idColumn, vehicleSamples, treatmentSamples, ...) standardGeneric("calculateFoldChange")
)
setMethod(
  f = "calculateFoldChange",
  signature = "ProteinAnalysisData",
  definition = function(object, idColumn, vehicleSamples, treatmentSamples, ...) {
    df <- object@dataframe
    if (baitID %in% idColumn) {
      # df[[column]] <- (df[[column]] - min(df[[column]], na.rm = TRUE)) /
      #                 (max(df[[column]], na.rm = TRUE) - min(df[[column]], na.rm = TRUE))
      # object@dataframe <- df
      # df = dplyr::filter(df, !(idColumn %in% cleanUpIDs))
      # object@dataframe = df

      # TODO: write this part of the function
    } else {
      warning("Bait ProteinID not found.")
    }
    return(object)
  }
)


# TODO: write a method to generate t-test p values
#' Run t-test
#'
#' @param object A ProteinAnalysisData object.
#' @param vehicleSamples The columns containing the normalized vehicle intensities (list).
#' @param treatmentSamples The columns containing the normalize treatment intensities (list).
#' @return The updated ProteinAnalysisData object.
#' @import dplyr
#' @export
setGeneric(
  "runCorrectedTTest",
  function(object, vehicleSamples, treatmentSamples, method = "Benjamini-Hochberg", ...) standardGeneric("runCorrectedTTest")
)
setMethod(
  f = "runCorrectedTTest",
  signature = "ProteinAnalysisData",
  definition = function(object, vehicleSamples, treatmentSamples, method = "Benjamini-Hochberg", ...) {
    df <- object@dataframe
    if (baitID %in% idColumn) {
      # df[[column]] <- (df[[column]] - min(df[[column]], na.rm = TRUE)) /
      #                 (max(df[[column]], na.rm = TRUE) - min(df[[column]], na.rm = TRUE))
      # object@dataframe <- df
      # df = dplyr::filter(df, !(idColumn %in% cleanUpIDs))
      # object@dataframe = df

      # TODO: write this part of the function
    } else {
      warning("Bait ProteinID not found.")
    }
    return(object)
  }
)


# TODO: write a method to generate volcano plots
#' Generate Volcano Plots
#'
#' @param object A ProteinAnalysisData object.
#' @param idColumn Column containing ProteinIDs
#' @param foldChange Column containing fold change values
#' @return A ggplot2 object
#' @import ggplot2
#' @export
setGeneric(
  "volcanoPlot",
  function(object, idColumn, foldChange, ...) standardGeneric("volcanoPlot")
)
setMethod(
  f = "volcanoPlot",
  signature = "ProteinAnalysisData",
  definition = function(object, idColumn, foldChange, ...) {
    df <- object@dataframe
    if (baitID %in% idColumn) {
      # df[[column]] <- (df[[column]] - min(df[[column]], na.rm = TRUE)) /
      #                 (max(df[[column]], na.rm = TRUE) - min(df[[column]], na.rm = TRUE))
      # object@dataframe <- df
      # df = dplyr::filter(df, !(idColumn %in% cleanUpIDs))
      # object@dataframe = df

      # TODO: write this part of the function
    } else {
      warning("Bait ProteinID not found.")
    }
    return(object)
  }
)
