#' deconvolution_fit2
#'
#' Function to calculate the number and proportion of each signature in each sample, in addition to plotting each sample's mutational profile and its signature contribution. 
#' @param input_matrices Palimpsest input list of mutation number and proportion matrices.
#' @param input_signatures Matrix of the mutational signatures to fit within the provided cohort of samples.
#' @param threshold Signatures contributing less then this percentage of total mutations within a sample will be discarded (e.g. if set to 6 and signature X contributes 5% of a sample's mutations, signature X will not be reported as being present in this sample).
#' @param sig_cols Character vector of R-compatable colours representing each signature to be used graphical outputs. Each signature in input_signatures must have named colour in this vector for grpahical outputs to work. 
#' @param doplot Logical indicating whether graphical outputs should be generated (defaults to TRUE). 
#' @param resdir Results directory.
#' @param save_signatures_exp Logical indicating whether or not signatures_exp object should be saved in the redsir (defaults to TRUE).
#' @keywords Signatures
#' @export
#' @import tibble
#' @import NMF
#' @examples
#' signatures_exp <- deconvolution_fit2(input_matrices = SBS_input,
#'                                          threshold = 6,input_signatures = SBS_liver,
#'                                            sig_cols = sig_cols,plot = T,resdir = resdir)
#' 

deconvolution_fit2 <- function (input_matrices = NULL,
                                 input_signatures = NULL, threshold = 6, sig_cols = NA,
                                 doplot = TRUE, save_signatures_exp = TRUE, resdir = resdir) {
  requireNamespace("tibble", quietly = TRUE);requireNamespace("NMF", quietly = TRUE)
  prop_matrix <- input_matrices$mut_props; num_matrix <- input_matrices$mut_nums
  if(nrow(prop_matrix) %!in% c(38,78,83,96)) stop("input_matrices format incorrect")
  if(nrow(prop_matrix)==96) Type <- "SBS"; if(nrow(prop_matrix)==78) Type <- "DBS"; if(nrow(prop_matrix)==83) Type <- "ID"; if(nrow(prop_matrix)==38) Type <- "SV" 
  
  resdir_parent <- resdir
  if (doplot == TRUE) {
    print(paste("Plotting the contribution of",Type,"signatures in each sample.."),quote = F)
    resdir <- paste0(resdir,"Samples/");if (!file.exists(resdir)) dir.create(resdir)
  }

  
  mutSign_props <- c()
  for (s in unique(colnames(prop_matrix))) {
    if (doplot == TRUE) {
      resdir.. <- file.path(resdir, s)
      if (!file.exists(resdir..)) {
        dir.create(resdir..)
      }
      pdf(file=paste0(resdir..,"/",s,"_",Type,"_profile.pdf"),width=24,height=7)
      Mean_plot_input <- as.matrix(prop_matrix[,s]); rownames(Mean_plot_input) <- rownames(prop_matrix)
      plot_signatures(input_data = Mean_plot_input, Title = paste(s))
      dev.off()
    }
    res <- fcnnls(as.matrix(t(input_signatures)), prop_matrix[,s], verbose = TRUE, pseudo = FALSE)
    sig.tot <- margin.table(res$x, 2)
    num.vec <- as.numeric(res$x)
    sig.fit <- data.frame(res$x)
    sig.fit$percent.fit <- c(sig.fit$res.x)/c(sig.tot) * 100
    sig.fit$sig_fit <- sig.fit$percent.fit
    sig.fit$sig_fit[sig.fit$sig_fit < threshold] <- 0
    mat.sigs <- data.frame(t(sig.fit$sig_fit))
    colnames(mat.sigs) <- rownames(res$x)
    mat.sigs$sums <- rowSums(mat.sigs)
    prop.sigs <- (mat.sigs/mat.sigs$sum) * 100
    prop.sigs <- prop.sigs[, -dim(prop.sigs)[2]]
    rownames(prop.sigs) <- s
    signature_content <- prop.sigs/100
    prop.sigs. <- data.frame(t(prop.sigs))
    prop.sigs. <- subset(prop.sigs., prop.sigs.[, 1] > 0)
    if (doplot == TRUE) {
      pdf(file.path(resdir.., paste0(s,"_",Type,"_Signature_Contribution.pdf")), width = 12, height = 10)
      pie(t(prop.sigs.), labels = rownames(prop.sigs.), 
          main = paste(Type,"Mutational Signatures Contribution in:", colnames(prop.sigs.)), 
          col = sig_cols[rownames(prop.sigs.)], 
          border = sig_cols[rownames(prop.sigs.)])
      dev.off()
    }
    mutSign_props <- rbind(mutSign_props, signature_content)
  }
  print("creating signatures_exp object",quote = F)
  
  mutSign_nums <- mutSign_props
  tot.muts <- as.data.frame(colSums(num_matrix)) %>%              ## similar method to Signature Analyzer
    rownames_to_column()
  colnames(tot.muts) <- c("Var1","Freq")
  mutSign_nums$total_mutations <- tot.muts[match(rownames(mutSign_nums), tot.muts$Var1), "Freq"]
  mutSign_nums <- round(mutSign_nums * mutSign_nums$total_mutations)
  mutSign_nums <- mutSign_nums[-dim(mutSign_nums)[2]]
  signatures_exp <- list(sig_props = mutSign_props, sig_nums = mutSign_nums)
  if(save_signatures_exp == TRUE) save(signatures_exp,file = paste0(resdir_parent,Type,"_signatures_exposure.RData"))
  return(signatures_exp)
  dev.off()
  
}