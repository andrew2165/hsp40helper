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
#' @param vehicleSamples The columns containing the normalized vehicle intensities (character vector; optional).
#' @param treatmentSamples The columns containing the normalized treatment intensities (character vector; optional).
#' @param method Correction method for multiple comparison (default: "Benjamini-Hochberg").
#' @return The updated ProteinAnalysisData object.
#' @import dplyr
#' @export
setGeneric(
  "runCorrectedTTest",
  function(object, vehicleSamples = NULL, treatmentSamples = NULL, method = "Benjamini-Hochberg", ...) standardGeneric("runCorrectedTTest")
)
setMethod(
  f = "runCorrectedTTest",
  signature = "ProteinAnalysisData",
  definition = function(object, vehicleSamples = NULL, treatmentSamples = NULL, method = "Benjamini-Hochberg", ...) {
    df <- object@dataframe

    # Use class slots if not provided as arguments
    if (is.null(vehicleSamples) && !is.null(slot(object, "vehicleSamples"))) {
      vehicleSamples <- slot(object, "vehicleSamples")
    }
    if (is.null(treatmentSamples) && !is.null(slot(object, "treatmentSamples"))) {
      treatmentSamples <- slot(object, "treatmentSamples")
    }
    if (is.null(vehicleSamples) || is.null(treatmentSamples)) {
      stop("vehicleSamples and treatmentSamples must be provided (either directly or via setSamples()).")
    }
    if (!all(vehicleSamples %in% colnames(df))) {
      stop("Some vehicleSamples columns not found in dataframe.")
    }
    if (!all(treatmentSamples %in% colnames(df))) {
      stop("Some treatmentSamples columns not found in dataframe.")
    }

    # Run t-test for each row/protein
    p.values <- apply(df, 1, function(row) {
      vehicle <- as.numeric(row[vehicleSamples])
      treatment <- as.numeric(row[treatmentSamples])
      if (length(na.omit(vehicle)) < 2 || length(na.omit(treatment)) < 2) return(NA)
      tryCatch({
        t.test(vehicle, treatment)$p.value
      }, error = function(e) NA)
    })

    # TODO: only apply this if that is specified in the function call
    # Benjamini-Hochberg correction
    ranks <- rank(p.values, ties.method = "min")
    m <- length(p.values)
    q.values <- p.values * m / ranks
    q.values[q.values > 1] <- 1

    df$TTestPValue <- p.values
    df$TTestQValue <- q.values
    object@dataframe <- df
    return(object)
  }
)


# TODO: write a method to generate volcano plots
#' Generate Volcano Plot
#'
#' @param object A ProteinAnalysisData object.
#' @param idColumn Column containing ProteinIDs (character scalar)
#' @param foldChange Column containing fold change values (character scalar)
#' @return A ggplot2 object
#' @import ggplot2
#' @export
setGeneric(
  "volcanoPlot",
  function(object, idColumn, foldChange = "FoldChange", pValueCol = "TTestPValue", ...) standardGeneric("volcanoPlot")
)
setMethod(
  f = "volcanoPlot",
  signature = "ProteinAnalysisData",
  definition = function(object, idColumn, foldChange = "FoldChange", pValueCol = "TTestPValue", ...) {
    df <- object@dataframe
    # Check necessary columns
    if (!(idColumn %in% colnames(df))) stop(paste("Column", idColumn, "not in dataframe"))
    if (!(foldChange %in% colnames(df))) stop(paste("Column", foldChange, "not in dataframe"))
    if (!(pValueCol %in% colnames(df))) stop(paste("Column", pValueCol, "not in dataframe"))

    # Prepare data, handling zeros and NAs
    plot_df <- df[!is.na(df[[foldChange]]) & !is.na(df[[pValueCol]]), ]
    plot_df$log2FC <- log2(plot_df[[foldChange]])
    # Avoid -Inf when p = 0; set very small value
    plot_df[[pValueCol]][plot_df[[pValueCol]] <= 0] <- .Machine$double.xmin
    plot_df$negLog10P <- -log10(plot_df[[pValueCol]])

    library(ggplot2)
    p <- ggplot(plot_df, aes(x = log2FC, y = negLog10P, label = plot_df[[idColumn]])) +
      geom_point(alpha = 0.6) +
      labs(
        x = "log2(Fold Change)",
        y = "-log10(p-value)",
        title = "Volcano Plot"
      ) +
      theme_minimal()

    return(p)
  }
)


# TODO: check all of the functions in this file for completeness and correctness