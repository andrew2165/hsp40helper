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
#' Standard set of ProteinIDs to remove (e.g., contaminants or unwanted proteins)
#'
#' @export
commonContaminantProteinIDs <- c(
  "InsertCommonContaminantsHere"
)


#' Remove unwanted ProteinIDs in ProteinAnalysisData
#'
#' Removes rows where the value in idColumn matches any of the provided cleanUpIDs.
#'
#' @param object A ProteinAnalysisData object.
#' @param idColumn Character. The column name containing Protein IDs.
#' @param cleanUpIDs Character vector of Protein IDs to remove.
#' @return The updated ProteinAnalysisData object.
#' @import dplyr
#' @export
setGeneric(
  "removeContaminants",
  function(object, idColumn, cleanUpIDs, ...) standardGeneric("removeContaminants")
)
setMethod(
  f = "removeContaminants",
  signature = "ProteinAnalysisData",
  definition = function(object, idColumn, cleanUpIDs, ...) {
    df <- object@dataframe
    if (idColumn %in% names(df)) {
      # Remove rows where the value in idColumn matches any cleanUpIDs
      df <- dplyr::filter(df, !(!!rlang::sym(idColumn) %in% cleanUpIDs))
      object@dataframe <- df
    } else {
      warning("Protein ID column not found in dataframe.")
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
  function(object, idColumn, baitID, removeBait = TRUE, ...) standardGeneric("baitNormalize")
)
setMethod(
  f = "baitNormalize",
  signature = "ProteinAnalysisData",
  definition = function(object, idColumn, baitID, removeBait = TRUE, ...) {
    df <- object@dataframe
    # Find the bait row
    bait_row <- df[[idColumn]] == baitID
    if (any(bait_row)) {
      # Identify sample columns (all columns except the idColumn)
      sample_cols <- setdiff(colnames(df), idColumn)
      # Get the bait intensities for all samples
      bait_intensities <- df[bait_row, sample_cols, drop = FALSE]
      if (nrow(bait_intensities) != 1) {
        warning("More than one row matched for bait protein ID; using the first.")
        bait_intensities <- bait_intensities[1, , drop = FALSE]
      }
      # Normalize each sample column by the bait intensity in that sample
      df[ , sample_cols] <- sweep(df[ , sample_cols, drop = FALSE], 2, as.numeric(bait_intensities), `/`)
      # Optionally remove the bait row
      if (removeBait) {
        df <- df[!bait_row, , drop = FALSE]
      }
      object@dataframe <- df
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
#' @param vehicleSamples The columns containing the normalized vehicle intensities (character vector; optional).
#' @param treatmentSamples The columns containing the normalized treatment intensities (character vector; optional).
#' @return The updated ProteinAnalysisData object.
#' @import dplyr
#' @export
setGeneric(
  "calculateFoldChange",
  function(object, idColumn, vehicleSamples = NULL, treatmentSamples = NULL, ...) standardGeneric("calculateFoldChange")
)
setMethod(
  f = "calculateFoldChange",
  signature = "ProteinAnalysisData",
  definition = function(object, idColumn, vehicleSamples = NULL, treatmentSamples = NULL, ...) {
    df <- object@dataframe
    # Use class slots if input not specified
    if (is.null(vehicleSamples) && !is.null(slot(object, "vehicleSamples"))) {
      vehicleSamples <- slot(object, "vehicleSamples")
    }
    if (is.null(treatmentSamples) && !is.null(slot(object, "treatmentSamples"))) {
      treatmentSamples <- slot(object, "treatmentSamples")
    }
    if (is.null(vehicleSamples) || is.null(treatmentSamples)) {
      stop("vehicleSamples and treatmentSamples must be provided (either as arguments or set with setSamples()).")
    }
    # Input checks
    if (!all(vehicleSamples %in% colnames(df))) {
      stop("Some vehicleSamples columns not found in dataframe.")
    }
    if (!all(treatmentSamples %in% colnames(df))) {
      stop("Some treatmentSamples columns not found in dataframe.")
    }

    # TODO: check for potential issue here where the vehicleSamples and treatmentSamples are either not
    # the normalized intensities bc of how the baitNormalize function works or check if the raw intensities
    # are being dropped from the dataframe entirely

    # Calculate means and fold change
    vehicle_means <- rowMeans(df[, vehicleSamples, drop=FALSE], na.rm=TRUE)
    treatment_means <- rowMeans(df[, treatmentSamples, drop=FALSE], na.rm=TRUE)
    df$FoldChange <- treatment_means / vehicle_means
    object@dataframe <- df
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


# TODO: check all of the functions in this file for completeness and correctness