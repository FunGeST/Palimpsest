#' NMF_Extraction
#'
#' Extracts mutational signatures de novo using NMF. Also estimates the optimal number of mutational signatures in the input. 
#' @param input_matrices Palimpsest input list of mutation number and proportion matrices.
#' @param range_of_sigs Numerical range of signatures. If "num_of_sigs" is set to auto, NMF will estiamte the optimal number of signatures in the input within this range before extracting with this number.
#' @param num_of_sigs The number of mutational signatures to extract. If left to the default "auto" value, the appropriate number of signatures will be estimated from NMF metrics (which are plotted in resdir). 
#' "auto" is recommended on the first run, but once the optimal number of signatures is known, setting this argument to that value will make the extraction quicker.
#' @param method Specification of the NMF algorithm. ‘brunet’ method corresponds to the standard NMF algorithm from Brunet et al. (2004, PNAS).
#' @param plot_sigs Logical, if TRUE (default), the extracted mutational signatures will be plotted in the resdir. 
#' @param resdir Results directory.

#' @keywords Signatures
#' @export
#' @import NMF
#' @examples
#' SBS_denovo_signatures <- NMF_Extraction(input_matrices = SBS_input, range_of_sigs = 1:20,nrun = 10,resdir = resdir)


NMF_Extraction <- function (input_matrices = NULL, range_of_sigs = 1:20, 
          num_of_sigs = "auto", nrun = 10, method = "brunet", plot_sigs = TRUE, 
          resdir = NA) 
{
  input_data <- as.matrix(input_matrices$mut_props)
  if(nrow(input_data) %!in% c(38,78,83,96)) stop("input_matrices format incorrect")
  if(nrow(input_data)==96) Type <- "SBS"; if(nrow(input_data)==78) Type <- "DBS"; if(nrow(input_data)==83) Type <- "ID"; if(nrow(input_data)==38) Type <- "SV"  
  requireNamespace("NMF", quietly = TRUE)
  sumRows <- rowSums(input_data)
  sort(sumRows)
  zeroes <- which(sumRows == 0)
  if (length(zeroes)) {
    input_data[zeroes, 1] <- 1e-10
  }
  if (num_of_sigs == "auto") {
    print(paste("Estimating the optimal number of",Type,"mutational signatures in the",ncol(input_data),"input samples.. (be patient!)"),quote = F)
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
    print(paste0("The optimal number of ",Type," signatures is ",estimated_rank,"."), quote = F)
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
      plot_signatures(input_data = spec,Title = rownames(spec))
      dev.off()
    }
  }
  denovo_signatures <- spec
  save(denovo_signatures,file = paste0(resdir,Type,"_denovo_signatures.RData"))
  return(denovo_signatures)
}


#' deconvolution_compare
#'
#' Function to calculate cosine similarity scores between two sets of mutational signatures (e.g. SBS COSMIC signatures vs Palimpsest de novo signatures).
#' @param new_signatures Data frame of de novo extracted mutational signatures
#' @param COSMIC_Signatures Data frame of reference mutational signatures (The reference SBS_cosmic, DBS_cosmic & ID_cosmic matrices are taken from Alexandrov et al. (2018)).
#'
#' @export
#' @import lsa
#' @examples
#' SBS_cosine_similarities <- deconvolution_compare(SBS_denovo_sigs,SBS_cosmic)

deconvolution_compare <- function(new_signatures, COSMIC_Signatures) {
  if(any(!grepl("denovo",rownames(new_signatures),ignore.case = T))) rownames(new_signatures) <- rep(paste("DeNovo", rownames(new_signatures), sep = "_"))
  mutmat <- t(rbind(new_signatures, COSMIC_Signatures))
  m <- as.matrix(palimpsest_distCosine(t(mutmat)))
  row_distance = as.dist(palimpsest_distCosine(t(m)))
  row_cluster = hclust(row_distance, method = "ward.D")
  col_distance = as.dist(palimpsest_distCosine(t(m)))
  col_cluster = hclust(col_distance, method = "ward.D")
  heatmap.2(1 - m, key.title = "Cosine similarity", keysize = 0.8, 
            main = " Cosine similarity Matrix", notecol = "black", 
            density.info = "none", trace = "none", margins = c(12, 
                                                               9), col = colorRampPalette(c("white", "yellow", "red"))(n = 299), 
            Rowv = as.dendrogram(row_cluster), Colv = as.dendrogram(col_cluster))
  return(1 - m)
}


#' signature_origins
#'
#' Annotates each mutation in a VCF with the signature with which it is most likely associated. E.g. when Type is set to "DBS", each line in the VCF corresponding to a DBS mutation is annotated with the probility that each DBS signature caused it.
#' @param vcf The input VCF file to which signature origin annotations are to be added. 
#' @param Type Mutation type (SBS, DBS, ID or SV).
#' @param signature_contribution List of signatures exposure numbers and proportions matrices (output from deconvolution_fit function).
#' @param input_signatures Matrix of the input signatures with which the VCF is to be annotated .
#' @keywords Signatures
#' @export
#' @examples
#' vcf <- signature_origins(vcf=vcf, Type = "SBS", signature_contribution = signatures_exp, input_signatures = SBS_liver_signatures)


signature_origins <- function (input = NULL, Type = Type,  
                               signature_contribution = signatures_exp, input_signatures = NULL){
  signature_contribution <- signature_contribution$sig_nums
  
  if(Type == "SBS") mutcat.col <- "SBS_cat3"; if(Type == "DBS") mutcat.col <- "DBS_cat"; if(Type == "ID") mutcat.col <- "ID_cat"; if(Type == "SV") mutcat.col <- "Category1"
  if(mutcat.col %!in% colnames(input)) stop(paste("vcf is missing the <<",mutcat.col,">> mutation category column. Use the << annotate_VCF >> function to add the appropriate column."))
  
  Sig.max.col <- paste0(Type,".Sig.max"); Sig.max.prob.col <- paste0(Type,".Sig.max.prob") 
  contrib <- t(signature_contribution)
  input[, paste(rownames(contrib), "prob", sep = ".")] <- NA
  input[,Sig.max.col] <- NA; input[,Sig.max.prob.col] <- NA
  input$unique <- paste0(input$Sample,input$CHROM,"_",input$POS); ordering <- input$unique
  vcf <- input[!is.na(input[,mutcat.col]),]
  output <- input[is.na(input[,mutcat.col]),]
  
  if (Type == "SBS") {
    bases <- c("A", "C", "G", "T")
    ctxt16 <- paste(rep(bases, each = 4), rep(bases, 4), 
                    sep = ".")
    mt <- c("CA", "CG", "CT", "TA", "TC", "TG")
    Types <- paste(rep(mt, each = 16), rep(ctxt16, 6), sep = "_")
    Types <- sapply(Types, function(z) {sub("\\.", substr(z, 1, 1), z)}) 
  }
  if (Type == "SV") {
    Types = c("BND_clust_0", "BND_clust_1", "DEL_0-1kb_clust_0", 
              "DEL_0-1kb_clust_1", "DEL_100kb-1Mb_clust_0", "DEL_100kb-1Mb_clust_1", 
              "DEL_10-100kb_clust_0", "DEL_10-100kb_clust_1", 
              "DEL_>10Mb_clust_0", "DEL_>10Mb_clust_1", "DEL_1-10kb_clust_0", 
              "DEL_1-10kb_clust_1", "DEL_1-10Mb_clust_0", "DEL_1-10Mb_clust_1", 
              "DUP_0-1kb_clust_0", "DUP_0-1kb_clust_1", "DUP_100kb-1Mb_clust_0", 
              "DUP_100kb-1Mb_clust_1", "DUP_10-100kb_clust_0", 
              "DUP_10-100kb_clust_1", "DUP_>10Mb_clust_0", "DUP_>10Mb_clust_1", 
              "DUP_1-10kb_clust_0", "DUP_1-10kb_clust_1", "DUP_1-10Mb_clust_0", 
              "DUP_1-10Mb_clust_1", "INV_0-1kb_clust_0", "INV_0-1kb_clust_1", 
              "INV_100kb-1Mb_clust_0", "INV_100kb-1Mb_clust_1", 
              "INV_10-100kb_clust_0", "INV_10-100kb_clust_1", 
              "INV_>10Mb_clust_0", "INV_>10Mb_clust_1", "INV_1-10kb_clust_0", 
              "INV_1-10kb_clust_1", "INV_1-10Mb_clust_0", "INV_1-10Mb_clust_1")
  }
  if(Type == "DBS"){
    Types = c("AC>CA", "AC>CG", "AC>CT", "AC>GA", "AC>GG", "AC>GT", "AC>TA", "AC>TG", "AC>TT", "AT>CA", "AT>CC", "AT>CG", "AT>GA", "AT>GC", "AT>TA", "CC>AA", "CC>AG",
              "CC>AT", "CC>GA" ,"CC>GG" ,"CC>GT", "CC>TA", "CC>TG", "CC>TT", "CG>AA", "CG>AC", "CG>AT", "CG>GA" ,"CG>GC" ,"CG>TA", "CT>AA", "CT>AC", "CT>AG", "CT>GA",
              "CT>GC" ,"CT>GG" ,"CT>TA" ,"CT>TC", "CT>TG", "GC>AA", "GC>AG", "GC>AT" ,"GC>CA" ,"GC>CG" ,"GC>TA" ,"TA>AC", "TA>AG", "TA>AT", "TA>CC", "TA>CG", "TA>GC",
              "TC>AA", "TC>AG", "TC>AT", "TC>CA", "TC>CG", "TC>CT", "TC>GA", "TC>GG", "TC>GT", "TG>AA", "TG>AC", "TG>AT", "TG>CA", "TG>CC", "TG>CT", "TG>GA", "TG>GC",
              "TG>GT", "TT>AA", "TT>AC", "TT>AG", "TT>CA", "TT>CC", "TT>CG", "TT>GA", "TT>GC", "TT>GG")
  }
  if(Type == "ID"){
    Types = c("DEL_C_1_0","DEL_C_1_1","DEL_C_1_2","DEL_C_1_3","DEL_C_1_4","DEL_C_1_5+", "DEL_T_1_0","DEL_T_1_1","DEL_T_1_2",
              "DEL_T_1_3","DEL_T_1_4","DEL_T_1_5+", "INS_C_1_0","INS_C_1_1","INS_C_1_2","INS_C_1_3","INS_C_1_4","INS_C_1_5+", 
              "INS_T_1_0","INS_T_1_1","INS_T_1_2","INS_T_1_3" ,"INS_T_1_4","INS_T_1_5+","DEL_repeats_2_0","DEL_repeats_2_1",
              "DEL_repeats_2_2","DEL_repeats_2_3","DEL_repeats_2_4","DEL_repeats_2_5+","DEL_repeats_3_0",
              "DEL_repeats_3_1","DEL_repeats_3_2","DEL_repeats_3_3","DEL_repeats_3_4","DEL_repeats_3_5+", 
              "DEL_repeats_4_0","DEL_repeats_4_1","DEL_repeats_4_2","DEL_repeats_4_3","DEL_repeats_4_4","DEL_repeats_4_5+",
              "DEL_repeats_5+_0","DEL_repeats_5+_1","DEL_repeats_5+_2","DEL_repeats_5+_3","DEL_repeats_5+_4","DEL_repeats_5+_5+",
              "INS_repeats_2_0","INS_repeats_2_1","INS_repeats_2_2","INS_repeats_2_3","INS_repeats_2_4","INS_repeats_2_5+" ,
              "INS_repeats_3_0","INS_repeats_3_1","INS_repeats_3_2","INS_repeats_3_3","INS_repeats_3_4","INS_repeats_3_5+" ,
              "INS_repeats_4_0","INS_repeats_4_1","INS_repeats_4_2","INS_repeats_4_3","INS_repeats_4_4","INS_repeats_4_5+", 
              "INS_repeats_5+_0","INS_repeats_5+_1","INS_repeats_5+_2","INS_repeats_5+_3","INS_repeats_5+_4","INS_repeats_5+_5+",
              "DEL_MH_2_1","DEL_MH_3_1","DEL_MH_3_2","DEL_MH_4_1","DEL_MH_4_2","DEL_MH_4_3",      
              "DEL_MH_5+_1","DEL_MH_5+_2","DEL_MH_5+_3","DEL_MH_5+_4","DEL_MH_5+_5+")
      Types <- gsub("_",",",Types)
  }
  if(Type == "SV"){
    Types = c("BND_clust_0","BND_clust_1","DEL_0-1kb_clust_0","DEL_0-1kb_clust_1","DEL_100kb-1Mb_clust_0",
              "DEL_100kb-1Mb_clust_1","DEL_10-100kb_clust_0","DEL_10-100kb_clust_1","DEL_>10Mb_clust_0",
              "DEL_>10Mb_clust_1","DEL_1-10kb_clust_0","DEL_1-10kb_clust_1","DEL_1-10Mb_clust_0",
              "DEL_1-10Mb_clust_1","DUP_0-1kb_clust_0","DUP_0-1kb_clust_1","DUP_100kb-1Mb_clust_0",
              "DUP_100kb-1Mb_clust_1","DUP_10-100kb_clust_0","DUP_10-100kb_clust_1","DUP_>10Mb_clust_0",
              "DUP_>10Mb_clust_1","DUP_1-10kb_clust_0","DUP_1-10kb_clust_1","DUP_1-10Mb_clust_0",
              "DUP_1-10Mb_clust_1","INV_0-1kb_clust_0","INV_0-1kb_clust_1","INV_100kb-1Mb_clust_0",
              "INV_100kb-1Mb_clust_1","INV_10-100kb_clust_0","INV_10-100kb_clust_1","INV_>10Mb_clust_0",
              "INV_>10Mb_clust_1","INV_1-10kb_clust_0","INV_1-10kb_clust_1","INV_1-10Mb_clust_0","INV_1-10Mb_clust_1")
  }
  
  sigs <- t(input_signatures)
  rownames(sigs) <- Types
  print("Calculating the probability of each mutation being due to each signature..",quote = F)
  total <- length(unique(vcf$Sample))
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  counter_start <- c()
  for (samp in unique(vcf$Sample)) {
    counter_start <- match(samp, unique(vcf$Sample))
    Sys.sleep(0.1)
    setTxtProgressBar(pb, counter_start)
    for (mut_Type in Types) {
      probs <- round(prop.table(sigs[mut_Type, ] * contrib[, 
                                                           samp]), digits = 4)
      if (all(is.na(probs))) {
        probs <- prop.table(contrib[, samp])
      }
      ind <- which(vcf$Sample == samp & vcf[, mutcat.col] == 
                     mut_Type)
      length(ind)
      for (sig in names(probs)) {
        vcf[ind, paste(sig, "prob", sep = ".")] <- probs[sig]
      }
    }
  }
  close(pb)
  print("Assigning the origin of each mutation..",quote = F)

  vcf[,Sig.max.col] <- sapply(1:nrow(vcf), function(i) {
    names(which.max(vcf[i, paste(rownames(contrib), "prob",
                                 sep = ".")]))
  })

  vcf[,Sig.max.col] <- substr(vcf[,Sig.max.col], 1, nchar(vcf[,Sig.max.col]) - 5)
  vcf[,Sig.max.prob.col] <- sapply(1:nrow(vcf), function(i){
    vcf[i,paste(names(which.max(vcf[i, paste(rownames(contrib), "prob", sep = ".")])))]
  })
  
  output <- rbind(output,vcf)
  output <- output[order(match(output$unique, ordering)), ] %>% 
    dplyr::select(-unique)
  return(output)
}



#' deconvolution_fit
#'
#' Function to calculate the number and proportion of each SBS, DBS or indel signature in each sample, in addition to plotting each sample's mutational profile and its signature contribution. 
#' @param input_matrices Palimpsest input list of mutation number and proportion matrices.
#' @param input_signatures Matrix of the mutational signatures to fit within the provided cohort of samples.
#' @param input_vcf The VCF used in the current analysis. Only required when Type == "SBS" so that the strand bias of SBS mutations can be plotted. 
#' @param threshold Signatures contributing less then this percentage of total mutations within a sample will be discarded (e.g. if set to 6 and signature X contributes 5 per cent of a sample's mutations, signature X will not be reported as present in this sample).
#' @param signature_colours Character vector of R-compatible colours representing each signature to be used graphical outputs. Each signature in input_signatures must have named colour in this vector for grpahical outputs to work. Use the "signature_colour_generator" function to generate colours for new signatures.
#' @param doplot Logical indicating whether graphical outputs should be generated (defaults to TRUE). 
#' @param resdir Results directory.
#' @param save_signatures_exp Logical indicating whether or not signatures_exp object should be saved in the redsir (defaults to TRUE).
#' @keywords Signatures
#' @export
#' @import tibble
#' @import NMF
#' @examples
#' signatures_exp <- deconvolution_fit(input_matrices = SBS_input, threshold = 8,input_signatures = SBS_liver,signature_colours = sig_cols,resdir = resdir)

deconvolution_fit <- function (input_matrices = NULL,
                                 input_signatures = NULL, input_vcf = vcf, threshold = 6, signature_colours = NA,
                                 doplot = TRUE, save_signatures_exp = TRUE, resdir = resdir) {
  requireNamespace("tibble", quietly = TRUE);requireNamespace("NMF", quietly = TRUE)
  prop_matrix <- input_matrices$mut_props; num_matrix <- input_matrices$mut_nums
  if(nrow(prop_matrix) %!in% c(38,78,83,96)) stop("input_matrices format incorrect")
  if(nrow(prop_matrix)==96) Type <- "SBS"; if(nrow(prop_matrix)==78) Type <- "DBS"; if(nrow(prop_matrix)==83) Type <- "ID"; if(nrow(prop_matrix)==38) stop("this isn't the SV function, please use 'deconvolution_fit_SV()' ")
  
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
      if(Type  == "SBS"){
        vcf_filt = filter(input_vcf, Sample == s & Type == "SNV")
        plotStrandBias96types(vcf_filt, mutcat3.col = "SBS_cat3",
          plot.file = paste0(resdir..,"/",s,"_Strand_bias_96_substitution_types.pdf"))
      }
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
          col = signature_colours[rownames(prop.sigs.)], 
          border = signature_colours[rownames(prop.sigs.)])
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


deconvolution_fit_SV <- function (vcf = vcf, input_data = data,
                               threshold = 5, input_signatures = NULL,
                               sig_cols = mycol, plot = TRUE, resdir = NULL)
{
  "%ni%" <- Negate("%in%");
  requireNamespace("NMF", quietly = TRUE)
  mutSign_props <- c()
  for (s in unique(vcf[, "Sample"])) {
    print(s)
    vcf. <- vcf[which(vcf[, "Sample"] == s), ]
    if (plot == TRUE) {
      resdir_samp <- paste0(resdir,"Samples/");if (!file.exists(resdir_samp)) dir.create(resdir_samp)
      resdir.. <- file.path(resdir_samp, s)
      if (!file.exists(resdir..)) {
        dir.create(resdir..)
      }
        datas_for_plot <- generate_table(vcf)
        translocs = datas_for_plot[[1]]
        others = datas_for_plot[[2]]
        df <- t(as.data.frame(input_data[, s]))
        rownames(df) <- s
        pdf(file.path(resdir.., "Mean_proportion_38_sv_types.pdf"),
            width = 24, height = 6)
        plot.SV.sigs(df)
        dev.off()
        make_plot(translocs, others, s,resdir = resdir..)
    }
    res <- fcnnls(as.matrix(t(input_signatures)), input_data[,s], verbose = TRUE, pseudo = FALSE)
    sig.tot <- margin.table(res$x, 2)
    num.vec <- as.numeric(res$x)
    sig.fit <- data.frame(res$x)
    sig.fit$percent.fit <- sig.fit$res.x/sig.tot * 100
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
    if (plot == TRUE) {
      pdf(file.path(resdir.., paste(s, "_Signature_Contribution.pdf", sep = "")), width = 12, height = 10)
      pie(t(prop.sigs.), labels = rownames(prop.sigs.),
          explode = 0.1, main = paste("Mutational Signatures Contribution in: ", colnames(prop.sigs.), sep = ""), col = sig_cols[rownames(prop.sigs.)],
          border = sig_cols[rownames(prop.sigs.)])
      dev.off()
    }
    mutSign_props <- rbind(mutSign_props, signature_content)
  }
  mutSign_nums <- mutSign_props
  tot.muts <- data.frame(table(vcf[, "Sample"]))
  mutSign_nums$total_mutations <- tot.muts[match(rownames(mutSign_nums),
                                                 tot.muts$Var1), "Freq"]
  mutSign_nums <- round(mutSign_nums * mutSign_nums$total_mutations)
  mutSign_nums <- mutSign_nums[-dim(mutSign_nums)[2]]
  outputs <- list(sig_props = mutSign_props, sig_nums = mutSign_nums)
  return(outputs)
}
