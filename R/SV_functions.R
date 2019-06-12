#' preprocessInput_sv
#'
#' Annotating the mutation data with necessary fields for further analysis
#' @param input_data Table describing somatic structural rearrangements.
#' @param ensgene Gene table for annotations. A table of Ensembl genes is provided with the package.
#' @param resdir Results directory where graphical outputs should be exported.
#'
#' @export
#'

preprocessInput_sv = function (input_data = NULL, Sample.col = "Sample", CHROM_1.col = "CHROM_1",
                                 CHROM_2.col = "CHROM_2", POS_1.col = "POS_1", POS_2.col = "POS_2",
                                 type.col = "Type",resdir = resdir, genome_build = "hg19")
{
  packs <- c("NMF", "SomaticSignatures")
  if ("NMF" %in% ls() | "SomaticSignatures" %in% ls()) {
    detach("package:NMF", unload = TRUE)
    detach("package:SomaticSignatures", unload = TRUE)
  }
  
  if(genome_build == "hg19") ensgene <- ensgene_hg19
  if(genome_build == "hg38") ensgene <- ensgene_hg38
  
  chroms <- unique(input_data[, CHROM_1.col])
  if (1 %in% chroms == TRUE) {
    input_data[, CHROM_1.col] <- paste("chr", input_data[,
                                                         CHROM_1.col], sep = "")
    input_data[, CHROM_2.col] <- paste("chr", input_data[,
                                                         CHROM_2.col], sep = "")
  }
  input_data <- input_data[mixedorder(input_data[, CHROM_1.col]),
                           ]
  input_data$strand.mut <- "+"
  input_data$strand.gene <- NA
  print("Annotating mutation data:")
  vcf <- palimpsest_dfPosXSegm_2(input_data, dfPos.chrom.col = CHROM_1.col,
                                 dfPos.pos.col = POS_1.col, ensgene, dfSegm.chrom.col = "Chromosome.Name",
                                 dfSegm.start.col = "Gene.Start..bp.", dfSegm.end.col = "Gene.End..bp.",
                                 colsToAdd = "Associated.Gene.Name", namesColsToAdd = "Associated.Gene.Name")
  vcf$strand.gene <- c("-", NA, "+")[ensgene[match(vcf$Associated.Gene.Name,
                                                   ensgene$Associated.Gene.Name), "Strand"] + 2]
  vcf <- palimpsest_addSVcategoriesToVcf(vcf, type.col = type.col,
                                           sample.col = Sample.col, CHROM_1.col = CHROM_1.col,
                                           CHROM_2.col = CHROM_2.col, POS_1.col = POS_1.col, POS_2.col = POS_2.col,
                                           resdir = resdir, genome_build = genome_build)
  return(vcf)
}




palimpsest_dfPosXSegm_2 = function (dfPos = NULL, dfPos.chrom.col = "chrom", dfPos.pos.col = "pos",
                                    dfSegm = NULL, dfSegm.chrom.col = "chrom", dfSegm.start.col = "start",
                                    dfSegm.end.col = "end", colsToAdd = NULL, namesColsToAdd = NULL,
                                    multseg = c(NA, "first", "all")[1])
{
  for (col in namesColsToAdd) dfPos[, col] <- NA
  dfSegm. <- split(dfSegm, dfSegm[, dfSegm.chrom.col])
  dfPos. <- split(dfPos, dfPos[, dfPos.chrom.col])
  chr_listy <- mixedsort(intersect(names(dfSegm.), names(dfPos.)))
  total <- length(chr_listy)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  counter_start <- c()
  if (is.na(multseg)) {
    for (chr in chr_listy) {
      counter_start <- match(chr, chr_listy)
      Sys.sleep(0.1)
      setTxtProgressBar(pb, counter_start)
      ind <- unlist(sapply(dfPos.[[chr]][, dfPos.pos.col],
                           function(pos) {
                             tmp <- which(dfSegm.[[chr]][, dfSegm.start.col] <=
                                            pos & dfSegm.[[chr]][, dfSegm.end.col] >=
                                            pos)
                             if (length(tmp) == 1)
                               tmp
                             else {
                               NA
                             }
                           }))
      if(!all(is.na(ind))){ ## Quentin added, sinon ca plante si y'a aucun SV sur toute la série trouve dans tout un chromosome (ca m'est arrive dans GEPELIN avec le chrY)
        dfPos.[[chr]][, namesColsToAdd] <- dfSegm.[[chr]][ind, colsToAdd]}
    }
    close(pb)
  }else{
    if (multseg == "first") {
      for (chr in chr_listy) {
        counter_start <- match(chr, chr_listy)
        Sys.sleep(0.1)
        setTxtProgressBar(pb, counter_start)
        ind <- unlist(sapply(dfPos.[[chr]][, dfPos.pos.col],
                             function(pos) {
                               tmp <- which(dfSegm.[[chr]][, dfSegm.start.col] <=
                                              pos & dfSegm.[[chr]][, dfSegm.end.col] >=
                                              pos)[1]
                             }))
        dfPos.[[chr]][, namesColsToAdd] <- dfSegm.[[chr]][ind,
                                                          colsToAdd]
      }
      close(pb)
    }
    if (multseg == "all") {
      for (chr in chr_listy) {
        counter_start <- match(chr, chr_listy)
        Sys.sleep(0.1)
        setTxtProgressBar(pb, counter_start)
        ind <- sapply(dfPos.[[chr]][, dfPos.pos.col],
                      function(pos) {
                        tmp <- which(dfSegm.[[chr]][, dfSegm.start.col] <=
                                       pos & dfSegm.[[chr]][, dfSegm.end.col] >=
                                       pos)
                        tmp
                      })
        for (j in 1:length(colsToAdd)) {
          dfPos.[[chr]][, namesColsToAdd[j]] <- unlist(lapply(ind,
                                                              function(z) paste(dfSegm.[[chr]][z, colsToAdd[j]],
                                                                                collapse = ",")))
          dfPos.[[chr]][which(sapply(ind, length) ==
                                0), namesColsToAdd[j]] <- NA
        }
      }
      close(pb)
    }
  }
  dfPos <- unsplit(dfPos., dfPos[, dfPos.chrom.col])
  dfPos
}


palimpsest_addSVcategoriesToVcf = function (sv = sv, type.col = "Type", sample.col = "Sample",
                                              CHROM_1.col = "CHROM", CHROM_2.col = "CHROM", POS_1.col = "POS",
                                              POS_2.col = "POS", resdir = resdir, genome_build =  genome_build)
{
  if(genome_build == "hg19") cyto <- cytoband_hg19
  if(genome_build == "hg38") cyto <- cytoband_hg38
  
  sv$Size <- 0
  sv$size.cat <- NA
  ind <- which(sv[, type.col] %in% c("DEL", "DUP", "INV"))
  sv[ind, "Size"] <- abs(sv[ind, POS_2.col] - sv[ind, POS_1.col]) ### ICI
  sv[which(sv$Size <= 1000), "size.cat"] <- "a_0-1kb"
  sv[which(sv$Size > 1000 & sv$Size <= 10000), "size.cat"] <- "b_1-10kb"
  sv[which(sv$Size > 10000 & sv$Size <= 1e+05), "size.cat"] <- "c_10-100kb"
  sv[which(sv$Size > 1e+05 & sv$Size <= 1e+06), "size.cat"] <- "d_100kb-1Mb"
  sv[which(sv$Size > 1e+06 & sv$Size <= 1e+07), "size.cat"] <- "e_1-10Mb"
  sv[which(sv$Size > 1e+07), "size.cat"] <- "f_>10Mb"
  sv$Clustered <- 0
  sv$Cluster <- NA
  minbkp <- 10
  workdir <- file.path(resdir, "bkp_clustering")
  if (!file.exists(workdir))
    dir.create(workdir)
  setwd(workdir)
  
  tmp_1 = paste(sv[,CHROM_1.col], sv[, POS_1.col])
  tmp_2 = paste(sv[, CHROM_2.col], sv[,POS_2.col])
  
  for (s in unique(sv[, sample.col])) {
    print(s)
    chr <- pos <- c()
    indsamp <- which(sv[, sample.col] == s)
    for (i in indsamp) {
      # if(sv[i,"Type"]=="BND"){    # ADDED \\\\\\\\\\\\\\\\\\\\\\\\\\\\\ si 2lignes de BND / trandloc (4/reciproque) ## formater si ca change
      #   chr <- c(chr,sv[i,CHROM_1.col])
      #   pos <- c(pos,sv[i,POS_1.col])
      # }else{
      #   chr <- c(chr,rep(sv[i,CHROM_1.col],2))
      #   pos <- c(pos,as.numeric(sv[i,c(POS_1.col,POS_2.col)]))
      # }  # ADDED ///////////////////////////////////////////////
      # # OLD
      chr <- c(chr, sv[i, CHROM_1.col], sv[i, CHROM_2.col])
      pos <- c(pos, as.numeric(sv[i, c(POS_1.col, POS_2.col)]))
    }
    bkp <- data.frame(chr, pos, dist = NA, stringsAsFactors = F)
    tmp <- sub("X", 23, sub("Y", 24, bkp$chr))
    bkp <- bkp[order(tmp, bkp$pos), ]
    for (c in unique(bkp$chr)) {
      ind <- which(bkp$chr == c)
      if (length(ind) > 1) {
        bkp[ind, ] <- bkp[ind[order(bkp[ind, "pos"])],
                          ]
        bkp[ind[-1], "dist"] <- bkp[ind[-1], "pos"] -
          bkp[ind[-length(ind)], "pos"]
      }
    }
    
    # bed <- bkp[, c(1, 2, 2, 3)]
    # indi <- bed2index(bed, sort = TRUE)
    # tmp <- cluster.region(x = indi, distance = "1000000",
    #   verbose = T)
    # tmp[, "chr"] <- data.frame(do.call("rbind", strsplit(as.character(tmp$index),
    #   ":", fixed = TRUE)))[1]
    # tmp[, c("start", "end")] <- data.frame(do.call("rbind",
    #   strsplit(as.character(tmp$index), "-", fixed = TRUE)))[2]
    # bkp$cluster <- tmp[, "regionIndex"]
    
    bed <- bkp[, c(1, 2, 2)]
    write.table(bed,"tmp.bed",sep="\t",quote=F,col.names=F,row.names=F) # added
    system(paste("bedtools cluster -i tmp.bed -d 1000000","> res.bed")) # added
    tmp <- read.delim("res.bed",as.is=T,header=F) # added
    bkp$cluster <- tmp[,4] # added
    
    bkp$clustered <- 0
    tt <- table(bkp$cluster)
    for (clst in names(which(tt >= 10))) {
      bkp[which(bkp$cluster == clst), "clustered"] <- 1
    }
    # # attention ici créer les tmp1 et tmp2 a chaque itération (cf loop modifiee) prend un temps fou...
    # for (i in which(bkp$clustered == 1)) {
    #   ind <- which(sv[, sample.col] == s &
    #        (paste(sv[,CHROM_1.col], sv[, POS_1.col]) == paste(bkp$chr[i],bkp$pos[i]) |
    #         paste(sv[, CHROM_2.col], sv[,POS_2.col]) == paste(bkp$chr[i], bkp$pos[i])))
    #   sv[ind, "Clustered"] <- 1
    #   sv[ind, "Cluster"] <- paste(s, bkp[i, "cluster"], sep = ".")
    #   # ind <- ind[which(sv[ind,"Type"]=="BND")]
    #   # for(k in ind){
    #   #   sv[which(sv$POS_2 == sv[k,"POS_1"] & sv$POS_1 == sv[k,"POS_2"]), "Clustered"] <- 1
    #   #   # sv[match(sv[k,"MATEID"],sv$ID),"Cluster"] <- paste(s,bkp[i,"cluster"],sep=".")
    #   # }
    # }
    
    tmp_0 = sv[, sample.col] == s
    # tmp_1 = paste(sv[,CHROM_1.col], sv[, POS_1.col]) # deplace avant la boucle principale
    # tmp_2 = paste(sv[, CHROM_2.col], sv[,POS_2.col])
    
    tmp_3 = paste(bkp$chr, bkp$pos)
    
    for (i in which(bkp$clustered == 1)) {
      ind <- which(tmp_0 &
                     (tmp_1 == tmp_3[i] |
                        tmp_2 == tmp_3[i]))
      sv[ind, "Clustered"] <- 1
      sv[ind, "Cluster"] <- paste(s, bkp[i, "cluster"], sep = ".")
    }
    
    chr <- pos <- clustered <- c()
    indsamp <- which(sv[, sample.col] == s)
    for (i in indsamp) {
      chr <- c(chr, sv[i, CHROM_1.col], sv[i, CHROM_2.col])
      pos <- c(pos, as.numeric(sv[i, c(POS_1.col, POS_2.col)]))
      clustered <- c(clustered, rep(sv[i, "Clustered"],
                                    2))
    }
    bkp <- data.frame(chr, pos, clustered, stringsAsFactors = F)
    bkp$clustered[is.na(bkp$clustered)] <- 0
    pdf(paste("bkp_", s, ".pdf", sep = ""), width = 15,
        height = 3)
    for (chr in unique(bkp$chr)) {
      ind <- which(bkp$chr == chr)
      if (length(ind) > 1)
        plot(bkp$pos[ind], rep(1, length(ind)), col = bkp$clustered[ind] +
               1, main = chr, type = "h", ylim = c(0, 1.5),
             yaxt = "n", bty = "l", las = 1, ylab = "",
             xlim = c(min(cyto[which(cyto$Chromosome ==
                                       sub("chr", "", chr)), "Start"]), max(cyto[which(cyto$Chromosome ==
                                                                                         sub("chr", "", chr)), "End"])), xlab = "Position")
      else next
    }
    dev.off()
  }
  sv$Clustered[is.na(sv$Clustered)] <- 0
  sv$Category1 <- paste(sv[, type.col], substr(sv$size.cat,
                                               3, nchar(sv$size.cat)), "clust", sv$Clustered, sep = "_")
  sv$Category1 <- sub("BND_0-1kb", "BND", sv$Category1)
  return(sv)
}


