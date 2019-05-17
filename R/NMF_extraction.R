#' NMF_Extraction
#'
#' Extracts mutational signatures de novo using NMF. Also estimates the optimal number of mutational signatures in the input. 
#' @param input_matrices Palimpsest input list of mutation number and proportion matrices.
#' @param Type Mutation type (SBS, DBS, ID or SV).
#' @param range_of_sigs Numerical range of signatures. If "num_of_sigs" is set to auto, NMF will estiamte the optimal number of signatures in the input within this range before extracting with this number.
#' @param num_of_sigs The number of mutational signatures to extract. If left to default "auto" value, the appropriate number of signatures will be estimated from NMF metrics (which are plotted in resdir). "auto" is recommended on the first run, but once the optimal number of signatures is known, setting this argument to that value will make the extraction quicker.
#' @param method Specification of the NMF algorithm. ‘brunet’ method corresponds to the standard NMF algorithm from Brunet et al. (2004, PNAS).
#' @param plot_sigs Logical, if TRUE (default), the extracted mutational signatures will be plotted in the resdir. 
#' @param resdir Results directory.

#' @keywords Signatures
#' @export
#' @import NMF
#' @examples
#' SBS_denovo_signatures <- NMF_Extraction(input_matrices = SBS_input,Type = "SBS",range_of_sigs = 1:20,nrun = 10, num_of_sigs = "auto", method = "brunet",plot_sigs=TRUE,resdir = resdir)


NMF_Extraction <- function (input_matrices = NULL, Type = NULL, range_of_sigs = NULL, 
          num_of_sigs = "auto", nrun = 10, method = "brunet", plot_sigs = TRUE, 
          resdir = NA) 
{
  input_data <- as.matrix(input_matrices$mut_props)
  requireNamespace("NMF", quietly = TRUE)
  sumRows <- rowSums(input_data)
  sort(sumRows)
  zeroes <- which(sumRows == 0)
  if (length(zeroes)) {
    input_data[zeroes, 1] <- 1e-10
  }
  if (num_of_sigs == "auto") {
    print(paste("Estimating the optimal number of mutational signatures in the",ncol(input_data),"input samples.. (be patient!)"),quote = F)
    estimate <- nmfEstimateRank(x = input_data, range_of_sigs, 
                                method = method, nrun = nrun, seed = 123456)
    pdf(paste0(resdir, "NMF_Rank_Estimates.pdf"), width = 8, 
        height = 6)
    p <- plot(estimate, y = NULL, what = "all", na.rm = FALSE, 
         xname = "x", yname = "y", xlab = "Factorization rank", 
         ylab = "", main = "NMF rank survey")
    print(p)
    dev.off()
    z <- estimate$measures$cophenetic[which(diff(estimate$measures$cophenetic) < 
                                              0)]
    steep_index <- which(estimate$measures$cophenetic == 
                           z[which(abs(diff(z)) == max(abs(diff(z))))])
    steep_index <- steep_index[length(steep_index)]
    estimated_rank <- estimate$measures$rank[steep_index]
  }
  else {
    estimated_rank <- num_of_sigs
  }
  print(paste("Extracting", estimated_rank, Type,"signatures from the input samples.."),quote = F)
  res <- nmf(input_data, rank = estimated_rank, method = method, 
             nrun = nrun, seed = 123456)
  sigs = t(basis(res))
  rownames(sigs) <- paste0(Type,"_denovo_", 1:nrow(sigs))
  spec <- sigs/rowSums(sigs)
  if (plot_sigs == TRUE) {
    if (Type == "SV") {
      pdf(file.path(resdir, "SV_Denovo_Signature_Profiles.pdf"), 
          width = 24, height = 5)
      plot.SV.sigs(spec)
      dev.off()
    }
    if (Type %in% c("SBS","DBS","ID")) {
      pdf(file=paste0(resdir,Type,"_Denovo_Signature_Profiles.pdf"),width=24,height=7)
      plot_signatures(input_data = spec,Title = rownames(spec),Type = Type,Individual = F,label = "Full")
      dev.off()
    }
  }
  denovo_signatures <- spec
  save(denovo_signatures,file = paste0(resdir,Type,"_denovo_signatures.RData"))
  return(denovo_signatures)
}
