#' Encode a 5-mer DNA Sequence as Categorical Position Columns
#'
#' Splits each DNA sequence into its individual nucleotide positions and
#' returns one factor column per position. This is the encoding used by the
#' random forest model shipped with the package, so the resulting columns can
#' be appended directly to a feature table before prediction.
#'
#' @param dna_strings Character vector of DNA sequences, all of the same
#'   length (five nucleotides for the m6A model, for example `"GGACA"`). The
#'   central nucleotide must match the context the model was trained on.
#'
#' @return A data frame with one row per input sequence and one factor column
#'   per nucleotide position (`nt_pos1`, `nt_pos2`, ...). Every factor has the
#'   levels `"A"`, `"T"`, `"C"` and `"G"`.
#'
#' @keywords internal
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A Sites for Multiple Candidate Sites
#'
#' Applies the fitted random forest model to a table of sequence and
#' structural features and returns the same table with one predicted
#' probability and one binary status per row. Candidate sites with a predicted
#' probability strictly greater than `positive_threshold` are labelled
#' `"Positive"`; a probability equal to the threshold is labelled
#' `"Negative"`.
#'
#' @param ml_fit A fitted random forest model, as supplied in the package file
#'   `inst/extdata/rf_fit.rds` and read with [readRDS()]. The model must have
#'   been trained with the features listed under `feature_df`.
#' @param feature_df A data frame with one row per candidate site and the
#'   columns `gc_content`, `RNA_type`, `RNA_region`, `exon_length`,
#'   `distance_to_junction`, `evolutionary_conservation` and `DNA_5mer`. The
#'   length and distance features are already log2-transformed in the course
#'   data, and `RNA_type` and `RNA_region` must use the labels `"mRNA"`,
#'   `"lincRNA"`, `"lncRNA"`, `"pseudogene"` and `"CDS"`, `"intron"`,
#'   `"3'UTR"`, `"5'UTR"` respectively.
#' @param positive_threshold Numeric between 0 and 1. The decision cut-off on
#'   the predicted probability of the Positive class. Defaults to `0.5`.
#'
#' @return `feature_df` with two added columns: `predicted_m6A_prob` (numeric,
#'   probability of the site being an m6A site) and `predicted_m6A_status`
#'   (character, either `"Positive"` or `"Negative"`).
#'
#' @seealso [prediction_single()] for a single site.
#'
#' @import randomForest
#' @importFrom stats predict
#'
#' @examples
#' # The model and the example feature table are bundled with the package.
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction", mustWork = TRUE))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                    package = "m6APrediction", mustWork = TRUE))
#'
#' # Predict the first five candidate sites.
#' batch_result <- prediction_multiple(ml_fit, example_df[1:5, ],
#'                                     positive_threshold = 0.6)
#' batch_result[, c("DNA_5mer", "predicted_m6A_prob", "predicted_m6A_status")]
#'
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}

#' Predict the m6A Status of a Single Site
#'
#' A convenience wrapper around [prediction_multiple()] for one candidate site.
#' The seven feature arguments are collected into a one-row data frame, scored
#' by the supplied model, and returned as a short named vector.
#'
#' @param ml_fit A fitted random forest model, as supplied in the package file
#'   `inst/extdata/rf_fit.rds`.
#' @param gc_content Numeric between 0 and 1, the GC content of the sequence
#'   window around the candidate site.
#' @param RNA_type Character, one of `"mRNA"`, `"lincRNA"`, `"lncRNA"` or
#'   `"pseudogene"`.
#' @param RNA_region Character, one of `"CDS"`, `"intron"`, `"3'UTR"` or
#'   `"5'UTR"`.
#' @param exon_length Numeric, the log2-transformed length of the hosting exon.
#' @param distance_to_junction Numeric, the log2-transformed distance from the
#'   candidate site to the nearest exon-exon junction.
#' @param evolutionary_conservation Numeric between 0 and 1, the conservation
#'   score of the site across species.
#' @param DNA_5mer Character, the five-nucleotide sequence context. The model
#'   expects an adenine in the central position, for example `"GGACA"`.
#' @param positive_threshold Numeric between 0 and 1, the decision cut-off on
#'   the predicted probability. Defaults to `0.5`; probabilities strictly
#'   greater than the threshold are labelled `"Positive"`.
#'
#' @return A named character vector with two elements: `predicted_m6A_prob`,
#'   the rounded probability of the site being an m6A site as a character
#'   string, and `predicted_m6A_status`, either `"Positive"` or `"Negative"`.
#'   The probability is coerced to character when the two elements are
#'   combined, so convert it with [as.numeric()] before doing arithmetic.
#'
#' @seealso [prediction_multiple()] for a table of sites.
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction", mustWork = TRUE))
#'
#' prediction_single(ml_fit,
#'                   gc_content = 0.6,
#'                   RNA_type = "mRNA",
#'                   RNA_region = "CDS",
#'                   exon_length = 12,
#'                   distance_to_junction = 5,
#'                   evolutionary_conservation = 0.8,
#'                   DNA_5mer = "GGACA",
#'                   positive_threshold = 0.5)
#'
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
