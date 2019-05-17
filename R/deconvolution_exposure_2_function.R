#' deconvolution_exposure_2
#'
#' Plots the exposure of signatures across the series (bargraphs of the number and proportion of each signature in each sample).
#' @param signature_contribution List of signatures exposure numbers and proportions matrices (output from deconvolution_fit function).
#' @param rm_samples List of names of samples to remove (if any). Once the signature contribution of hypermutated samples is known, removing these from the grpahical output of this function can make the comparison of the rest of the samples easier. 
#' @param sig_cols Character vector of R-compatable colours representing each signature to be used graphical outputs. Each signature in input_signatures must have named colour in this vector for grpahical outputs to work. 
#' @keywords Signatures
#' @export
#' @import ggpubr
#' @import reshape2
#' @examples
#' pdf(file.path(resdir.,"signature_content_plot.pdf"),width=12,height=7)
#' signature_content_plot <- deconvolution_exposure_2(signature_contribution = SBS_signatures_exp, rm_samples = "CHC892T", sig_cols = sig_cols)
#' dev.off()

deconvolution_exposure_2 <- function(signature_contribution = signatures_exp, rm_samples = NA,sig_cols = NA){
  
  requireNamespace("ggpubr", quietly = TRUE); requireNamespace("reshape2", quietly = TRUE)
  sig_nums <- signature_contribution$sig_nums; sig_props <- signature_contribution$sig_props
  if(!is.na(rm_samples)){
    sig_nums <- sig_nums[which(rownames(sig_nums) %!in% rm_samples),]
    sig_props <- sig_props[which(rownames(sig_props) %!in% rm_samples),]
  }
  
  nsamp <- nrow(sig_nums)
  nrows <- ceiling(nsamp/200)
  signature_content_plot <- c()
  if(nrows == 1){
    content_plot <- signature_exposure_plot(sig_nums, sig_props, sig_cols = sig_cols)
  }
  if(nrows>1){  ## for >200 samples we find it best to sort the plot into several rows to make it easier to view, the rest of the function does this
    
    nsamp_new <- nsamp
    for(j in 1:20){    ## making number of samples divisible by the number of rows 
      if( (nsamp_new/nrows)%%1 == 0){
        break
      }
      if( (nsamp_new/nrows)%%1 != 0){
        nsamp_new <- nsamp_new+1
      }
    }
    if(nsamp != nsamp_new){    ## adding "blank" samples to make number of samples up to number calculated in the previous step
      for(i in (1+nsamp):nsamp_new){
        sig_nums[i,] <- 0
        sig_props[i,] <- 0
        rownames(sig_nums)[i] <- paste(i)
        rownames(sig_props)[i] <- paste(i)
      }
    }
    
    for(i in 1:nrows){    ## sort samples into groups, one for each row
      n4me <- paste0("row_",i)
    
      if(i == 1) n4mes <- names(sort(rowSums(sig_nums)))[1:(nsamp_new/(nrows/i))]
      if(i > 1) n4mes <- names(sort(rowSums(sig_nums)))[((nsamp_new/(nrows/(i-1)))+1):(nsamp_new*i/nrows)]
      sig_nums_tmp <- sig_nums[n4mes,];sig_props_tmp <- sig_props[n4mes,]
      
      signature_content_plot[[i]] <- signature_exposure_plot(sig_nums_tmp,sig_props_tmp,sig_cols = sig_cols)

    }
    content_plot <-  ggarrange(plotlist =  signature_content_plot[nrows:1], heights = rep(7,nrows), widths = rep(10,nrows),
              ncol = 1, nrow = nrows)
  }
  print(content_plot)
}


#' signature_exposure_plot
#'
#' Function to plot the exposure of mutational signatures in samples across the entire series
#' @param mutSign_nums Matrix in sample x mutational signature exposure format in numbers
#' @param mutSign_props Matrix in sample x mutational signature exposure format in proportions
#' @param sig_cols Character vector indicating the colors representing each signature in graphical outputs. Must match to the total number of provided signatures
#'
#' @export
#' @import ggplot2
#' @importFrom ggplot2 theme_bw
#' @import gplots
#' @import reshape2
#' @importFrom reshape2 melt


signature_exposure_plot <- function (mutSign_nums, mutSign_props, sig_cols = NULL) 
{
  scale <- 1
  .theme_ss <- theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 90, 
                                                                           vjust = 0.5, size = 8 * scale, family = "mono"), axis.text.y = element_text(hjust = 0.5, 
                                                                                                                                                       size = 12 * scale, family = "mono"), axis.text = element_text(size = 12 * 
                                                                                                                                                                                                                       scale, family = "mono"))
  ordering <- order(colSums(t(mutSign_nums)), decreasing = T)
  mutSign_nums <- t(mutSign_nums)
  mutSign_props <- t(mutSign_props)
  mutSign_nums <- mutSign_nums[, ordering]
  mutSign_props <- mutSign_props[, ordering]
  sample.ordering <- colnames(mutSign_nums)
  x1 <- melt(mutSign_nums)
  x2 <- melt(mutSign_props)
  colnames(x1) <- c("Signature", "Sample", "Activity")
  colnames(x2) <- c("Signature", "Sample", "Activity")
  x1[, "class0"] <- c("Counts")
  x2[, "class0"] <- c("Proportions")
  df2 <- rbind(x1, x2)
  df2$class0 <- factor(df2$class0, c("Counts", "Proportions"))
  df2$Sample <- factor(df2$Sample, sample.ordering)
  p = ggplot(df2, aes(x = factor(Sample), y = Activity, fill = Signature))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + scale_fill_manual(values = sig_cols)
  p = p + ggtitle("Mutational Signature Exposures")
  p = p + facet_grid(class0 ~ ., scale = "free_y")
  p = p + theme(plot.title = element_text(lineheight = 1, face = "bold", 
                                          size = 15 * scale))
  p = p + xlab("Samples") + ylab("Mutational Signature Content")
  p = p + theme(axis.title.x = element_text(face = "bold", 
                                            colour = "black", size = 15 * scale))
  p = p + theme(axis.title.y = element_text(face = "bold", 
                                            colour = "black", size = 15 * scale))
  p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                           size = 12 * scale, face = "bold", colour = "black"))
  p = p + theme(axis.text.y = element_text(size = 10 * scale, 
                                           face = "bold", colour = "black"))
  p = p + theme(legend.title = element_blank())
  p = p + .theme_ss
  p = p + theme(legend.position = "top")
  return(p)
}
