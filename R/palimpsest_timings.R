#' chrTime_annot
#'
#' Function to annotate the mutation data with chromosomal duplication timing information
#' @param vcf vcf data frame with mutation information to annotate
#' @param cna_data Data frame with Copy Number Alteration Information
#' @param cyto Cytoband description table
#'
#' @return
#' @export
#'
#' @examples
chrTime_annot <- function (vcf = NULL, 
						   cna_data = NULL,
                           cyto = cyto)
{
  sample.col = "Sample"; CHROM.col = "CHROM";POS.col = "POS"; REF.col = "REF"; ALT.col = "ALT"; ntot.col = "ntot";
  nMin.col = "Nmin";nMaj.col = "Nmaj"; cnv.chr = "CHROM"; cnv.start = "POS_START"; cnv.end = "POS_END"
  cnv.LogR = "LogR";cnv.ntot = "ntot"; cnv.sample = "Sample"

  # Annotate mutation and CNA tables with chromosome arms
  cna_data$startarm <- NA
  cna_data$endarm <- NA
  for (a in unique(cyto$Arm)) {
    chr <- paste("chr", sub("p|q", "", a), sep = "")
    ind <- which(cyto$Arm == a)
    minpos <- min(cyto[ind, "Start"])
    maxpos <- max(cyto[ind, "End"])
    vcf[which(vcf[, CHROM.col] == chr & vcf[, POS.col] >= minpos & vcf[, POS.col] <= maxpos), "Arm"] <- a
    cna_data[which(cna_data[, cnv.chr] == chr & cna_data[, cnv.start] >= minpos & cna_data[, cnv.start] <= maxpos), "startarm"] <- a
    cna_data[which(cna_data[, cnv.chr] == chr & cna_data[, cnv.end] >= minpos & cna_data[, cnv.end] <= maxpos), "endarm"] <- a
  }
  
  # Identify aberrant chromosome segments and classify as deletions/gains/UPD
  cna_data$type <- "normal"
  cna_data$dup.lev <- NA
  cna_data$dup.ID <- NA
  cna_data[which(cna_data[, nMaj.col] <= 1 & cna_data[, nMin.col] == 0),"type"] <- "deletion"
  cna_data[which(cna_data[, nMaj.col] > 1 & cna_data[, nMin.col] == 0),"type"] <- "UPD"
  cna_data[which(cna_data[, nMaj.col] > 1 & cna_data[, nMin.col] > 0),"type"] <- "gain"
  ind <- which(cna_data[, nMaj.col] > 1)
  cna_data[ind, "dup.lev"] <- paste(cna_data[ind, "type"], round(cna_data[ind, ntot.col]), cna_data[ind, nMaj.col], cna_data[ind, nMin.col], sep = "_")
  cna_data[ind, "dup.ID"] <- paste(cna_data[ind, "Arm"], cna_data[ind, "dup.lev"], sep = "_")
	
  vcf$type <- "normal"
  vcf$dup.lev <- NA
  vcf$dup.ID <- NA
  vcf[which(vcf[, nMaj.col] <= 1 & vcf[, nMin.col] == 0),"type"] <- "deletion"
  vcf[which(vcf[, nMaj.col] > 1 & vcf[, nMin.col] == 0),"type"] <- "UPD"
  vcf[which(vcf[, nMaj.col] > 1 & vcf[, nMin.col] > 0),"type"] <- "gain"
  ind <- which(vcf[, nMaj.col] > 1)
  vcf[ind, "dup.lev"] <- paste(vcf[ind, "type"], round(vcf[ind, ntot.col]), vcf[ind, nMaj.col], vcf[ind, nMin.col], sep = "_")
  vcf[ind, "dup.ID"] <- paste(vcf[ind, "Arm"], vcf[ind, "dup.lev"], sep = "_")

  # Time chromosome duplications
  nmut.min <- 30 # minimum number of mutations to be able to time a segment
  colclust <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
                "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6",
                "#6A3D9A", "#FFFF99", "#B15928", "#bd18ea", "#2ef4ca",
                "#f4cced", "#f4cc03", "#05188a", "#e5a25a", "#06f106",
                "#85848f", "#000000", "#076f25", "#93cd7f", "#4d0776",
                "darkred", "darkgrey")
  point.mut.time <- data.frame(Sample = c(), ID = c(), CNA = c(),
                               Chrom = c(), Type = c(), Ntot = c(), Nmaj = c(), Nmin = c(),
                               Ndup = c(), Nsing = c(), Moltime = c(), col = c())
  vcf$uniqueID <- paste(vcf[, sample.col], vcf[, CHROM.col], vcf[, POS.col], vcf[, REF.col], vcf[, ALT.col], sep = " ")
  vcf$Timing <- NA
  for (s in unique(vcf[, sample.col])) {
    print(s)
    memvcf <- vcf[which(vcf[, sample.col] == s), ]
    vcf. <- vcf[which(vcf[, sample.col] == s & !is.na(vcf$dup.lev)),]
    if (!nrow(vcf.)) { next }
    tt <- table(vcf.[, c("Arm", "dup.lev")])
    mydups <- c()
    for (a in rownames(tt)) {
      mydups <- c(mydups, paste(a, names(which.max(tt[a,])), sep = "_"))
    }
    nmut <- sapply(mydups, function(z) sum(vcf.$dup.ID == z))
    mydups <- mydups[which(nmut >= nmut.min)]
    if (!length(mydups)) { next }
    mychrs <- sub("p|q", "", mydups)
    pq <- names(which(table(mychrs) > 1))
    if (length(pq)) {
      ind <- which(mychrs %in% pq)
      for (i in ind) {
        vcf.[which(vcf.$dup.ID == mydups[i]), "Arm"] <- sub("p|q", "", vcf.[which(vcf.$dup.ID == mydups[i]), "Arm"])
        vcf[which(vcf$dup.ID == mydups[i]), "Arm"] <- sub("p|q", "", vcf[which(vcf$dup.ID == mydups[i]), "Arm"])
        vcf.[which(vcf.$dup.ID == mydups[i]), "dup.ID"] <- sub("p|q", "", vcf.[which(vcf.$dup.ID == mydups[i]), "dup.ID"])
        vcf[which(vcf$dup.ID == mydups[i]), "dup.ID"] <- sub("p|q", "", vcf[which(vcf$dup.ID == mydups[i]), "dup.ID"])
        mydups[i] <- sub("p|q", "", mydups[i])
      }
    }
    mydups <- unique(mydups)
    vcf. <- vcf.[which(vcf.$dup.ID %in% mydups), ]
    ind <- match(unique(vcf.$dup.ID), vcf.$dup.ID)
    duptab <- data.frame(Sample = s, ID = vcf.[ind, "dup.ID"], CNA = paste(vcf.$dup.type, vcf.$Arm)[ind], Chrom = vcf.[ind,CHROM.col], Type = vcf.[ind, "type"], Ntot = round(vcf.[ind,ntot.col]), Nmaj = vcf.[ind, nMaj.col], Nmin = vcf.[ind, nMin.col], Ndup = NA, Nsing = NA, Moltime = NA)
    vcf.$Timing <- NA
    for (i in which(duptab$Type == "gain")) {
      ind <- which(vcf.$dup.ID == duptab[i, "ID"])
      tmp <- cit.peaks(vcf.[ind, "nmut"], percentHighestPeak = 0.01, bw.adj = 1.5)
      maxpeak <- tmp$`nb big peaks`
      vcf.[ind[which(vcf.[ind, "nmut"] >= tmp$`peaks x size`[maxpeak, "left born"] & vcf.[ind, "nmut"] <= tmp$`peaks x size`[maxpeak, "right born"])], "Timing"] <- "early"
      vcf.[ind[which(vcf.[ind, "mult"] >= 1 & vcf.[ind, "nmut"] < tmp$`peaks x size`[maxpeak, "left born"])], "Timing"] <- "late"
      duptab[i, "Ndup"] <- sum(vcf.[ind, "Timing"] == "early", na.rm = T)
      duptab[i, "Nsing"] <- sum(vcf.[ind, "Timing"] == "late", na.rm = T)
      if (duptab[i, "Nmaj"] == duptab[i, "Nmin"]) { npardup <- 2 } else { npardup <- 1 }
      duptab[i, "Moltime"] <- min(c(100, (duptab[i, "Ndup"]/npardup)/((duptab[i,"Ndup"]/npardup) + (duptab[i, "Nsing"] - duptab[i,"Ndup"] * (2 - npardup))/mean(c(3, duptab[i,"Ntot"]))) * 100))
    }
    for (i in which(duptab$Type == "UPD")) {
      ind <- which(vcf.$dup.ID == duptab[i, "ID"])
      tmp <- cit.peaks(vcf.[ind, "nmut"])
      maxpeak <- tmp$`nb big peaks`
      vcf.[ind[which(vcf.[ind, "nmut"] >= tmp$`peaks x size`[maxpeak, "left born"] & vcf.[ind, "nmut"] <= tmp$`peaks x size`[maxpeak, "right born"])], "Timing"] <- "early"
      vcf.[ind[which(vcf.[ind, "mult"] >= 1 & vcf.[ind, "nmut"] < tmp$`peaks x size`[maxpeak, "left born"])], "Timing"] <- "late"
      duptab[i, "Ndup"] <- sum(vcf.[ind, "Timing"] == "early", na.rm = T)
      duptab[i, "Nsing"] <- sum(vcf.[ind, "Timing"] == "late", na.rm = T)
      duptab[i, "Moltime"] <- min(c(100, duptab[i, "Ndup"]/(duptab[i, "Ndup"] + duptab[i, "Nsing"]/mean(c(2, duptab[i, "Ntot"]))) * 100))
    }
    duptab$Col <- colclust[match(sub("chr", "", duptab$Chrom), c(1:22, "X", "Y"))]
    point.mut.time <- rbind(point.mut.time, duptab)
    vcf[match(vcf.$uniqueID, vcf$uniqueID), "Timing"] <- vcf.$Timing
  }
  
  outputs <- list(vcf = vcf, point.mut.time = point.mut.time, cna_data=cna_data)
  return(outputs)
}
