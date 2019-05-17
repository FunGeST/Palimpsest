#' signature_colour_generator
#'
#' Generates colours for de novo signatures.
#' @param signature_names Character vector of the names of the mutational signatures for which colours are to be generated.
#' @keywords Signatures
#' @export
#' @import RColorBrewer
#' @examples
#' SBS_colours <- signature_colour_generator(signature_names = rownames(SBS_denovo_sigs))



signature_colour_generator <- function(signature_names = NULL){
  requireNamespace("RColorBrewer",quietly = T)
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  sig_cols <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  sig_cols <- sig_cols[sample.int(length(sig_cols),length(signature_names))];names(sig_cols) <- signature_names
  return(sig_cols)
}