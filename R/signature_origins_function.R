#' signature_origins
#'
#' Annotates each mutation in a VCF with the signature with which it is most likely associated. E.g. when Type is set to "DBS", each line in the VCF corresponding to a DBS mutation is annotated with the probility that each DBS signature caused it.
#' @param vcf The input VCF file to which signature origin annotations are to be added. 
#' @param Type Mutation type (SBS, DBS, ID or SV)
#' @param signature_contribution List of signatures exposure numbers and proportions matrices (output from deconvolution_fit function).
#' @param input_signatures Matrix of the input signatures with which the VCF is to be annotated 
#' @keywords Signatures
#' @export
#' @examples
#' vcf <- signature_origins(vcf=vcf, Type = "SBS", signature_contribution = signatures_exp, input_signatures = SBS_liver_signatures)


signature_origins <- function (input = NULL, Type = Type,  
                               signature_contribution = signatures_exp, input_signatures = NULL){
  signature_contribution <- signature_contribution$sig_nums
  
  if(Type == "SBS") mutcat.col <- "SBS_cat3"; if(Type == "DBS") mutcat.col <- "DBS_cat"; if(Type == "ID") mutcat.col <- "ID_cat"
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

