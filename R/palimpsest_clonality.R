#' cnaCCF_annot
#'
#' This function allows to annotate the input mutation data with clonality information (Cancer Cell Fraction- CCF) based on the Copy Number information provided
#' 
#' @param vcf vcf data frame containing somatic mutations
#' @param annot_data Annotation file with Gender information and Tumor Purity for each sample
#' @param cna_data Data frame with Copy Number Alteration Information
#' @param CCF_boundary percentage boundary threshold for clonal mutations
#'
#' @return
#' @export
#'
#' @examples
cnaCCF_annot <- function (vcf = NULL, annot_data = NULL, cna_data = NULL, CCF_boundary = 0.95)
{
  sample.col = "Sample"; CHROM.col = "CHROM"; POS.col = "POS"; VAF.col = "Tumor_Varprop"; VarCount.col = "Tumor_Varcount";
  TumorDepth.col = "Tumor_Depth"; NormalDepth.col = "Normal_Depth";cnv.chr = "CHROM"; cnv.start = "POS_START";cnv.end = "POS_END";
  cnv.LogR = "LogR";cnv.ntot = "ntot";cnv.nMaj = "Nmaj"; cnv.nMin = "Nmin"; cnv.sample = "Sample"
  vcf_all <- c()
  vcf[,VAF.col] <- vcf[,VarCount.col]/vcf[,TumorDepth.col] 
  vcf$Gender <- annot_data[match(vcf[, "Sample"], annot_data[, "Sample"]), "Gender"]
  vcf$Purity <- annot_data[match(vcf[, "Sample"], annot_data[, "Sample"]), "Purity"]
  for (samp_name in unique(vcf[, "Sample"])) {
    print(samp_name)
    cnv_samp <- cna_data[which(cna_data[, cnv.sample] == samp_name), ]
    #cnv_samp$chr <- paste("chr", cnv_samp[, cnv.chr], sep = "")
    vcf. <- vcf[which(vcf[, sample.col] == samp_name), ]
    rho <- unique(vcf.[, "Purity"])
    vcf. <- palimpsest_dfPosXSegm(vcf., CHROM.col, POS.col, cnv_samp, cnv.chr, cnv.start, cnv.end, c(cnv.LogR, cnv.ntot, cnv.nMaj, cnv.nMin), c(cnv.LogR, cnv.ntot, cnv.nMaj, cnv.nMin))
    vcf.$LogR <- log2(vcf.[, TumorDepth.col]/vcf.[, NormalDepth.col])
    if (unique(vcf.$Gender) == "M") {
      ind <- which(vcf.[, CHROM.col] == "chrX")
      vcf.$nmut <- vcf.[, VAF.col] * (rho * vcf.[, cnv.ntot] + (1 - rho) * 2)/rho
      vcf.[ind, "nmut"] <- vcf.[ind, VAF.col] * (rho * vcf.[ind, cnv.ntot] + (1 - rho) * 1)/rho
    }
    if (unique(vcf.$Gender) != "M") {
      vcf.$nmut <- vcf.[, VAF.col] * (rho * vcf.[, cnv.ntot] + (1 - rho) * 2)/rho
    }
    vcf.$mult <- round(vcf.$nmut)
    vcf.$CCF <- vcf.$nmut
    ind <- which(vcf.$mult > 1.5)
    vcf.[ind, "CCF"] <- vcf.[ind, "CCF"]/vcf.[ind, "mult"]
    vcf.$CCF.adj <- vcf.$CCF
    vcf.[which(vcf.$CCF > 1), "CCF.adj"] <- 1
    ccf.conf <- sapply(1:nrow(vcf.), function(i) binom.test(vcf.[i,VarCount.col], vcf.[i, TumorDepth.col])$conf.int)
    vcf.$CCF.min <- ccf.conf[1, ] * (rho * vcf.[, cnv.ntot] + (1 - rho) * 2)/rho
    vcf.[ind, "CCF.min"] <- vcf.[ind, "CCF.min"]/vcf.[ind,"mult"]
    vcf.$CCF.max <- ccf.conf[2, ] * (rho * vcf.[, cnv.ntot] + (1 - rho) * 2)/rho
    vcf.[ind, "CCF.max"] <- vcf.[ind, "CCF.max"]/vcf.[ind,"mult"]
    vcf.$Clonality <- c("subclonal", "clonal")[as.numeric(vcf.$CCF.max >= CCF_boundary) + 1]
    vcf_all <- rbind(vcf_all, vcf.)
  }
  return(vcf_all)
}
