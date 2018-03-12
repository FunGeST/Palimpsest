#' palimpsest_dfPosXSegm
#'
#' Function to annotate positions in a data frame (dfPos) using segments in another data frame (dfSegm)
#' @param dfPos Data frame with positions to annotate
#' @param dfPos.chrom.col Chromosome column in dfPos
#' @param dfPos.pos.col Position column in dfPos
#' @param dfSegm Data frame with segments to use for annotating dfPos
#' @param dfSegm.chrom.col Chromosome column in dfSegm
#' @param dfSegm.start.col Start position column in dfSegm
#' @param dfSegm.end.col End position column in dfSegm
#' @param colsToAdd Names of columns in dfSegm that should be used to annotate dfPos
#' @param namesColsToAdd Column names to give in the columns added to dfPos
#' @param multseg required
#'
#' @return
#' @export
#' @import gtools
#' @examples
palimpsest_dfPosXSegm <- function(dfPos=NULL,
                               dfPos.chrom.col="chrom",
                               dfPos.pos.col="pos",
                               dfSegm=NULL,
                               dfSegm.chrom.col="chrom",
                               dfSegm.start.col="start",
                               dfSegm.end.col="end",
                               colsToAdd=NULL,
                               namesColsToAdd=NULL,
                               multseg=c(NA,"first","all")[1]
)
{
  for(col in namesColsToAdd)	dfPos[,col] <- NA
  dfSegm. <- split(dfSegm,dfSegm[,dfSegm.chrom.col])
  dfPos. <- split(dfPos,dfPos[,dfPos.chrom.col])
  chr_listy <- mixedsort(intersect(names(dfSegm.),names(dfPos.)))
  total <- length(chr_listy)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  counter_start <- c()
  if(is.na(multseg)){
    for(chr in chr_listy)
    {
      counter_start <- match(chr,chr_listy)
      Sys.sleep(0.1)
      # update progress bar
      setTxtProgressBar(pb, counter_start)

      ind <- unlist(sapply(dfPos.[[chr]][,dfPos.pos.col],function(pos){
        tmp <- which(dfSegm.[[chr]][,dfSegm.start.col] <= pos & dfSegm.[[chr]][,dfSegm.end.col] >= pos)
        if(length(tmp)==1)	tmp else{NA}
      }))
      dfPos.[[chr]][,namesColsToAdd] <- dfSegm.[[chr]][ind,colsToAdd]
    }
    close(pb)

  }else{
    if(multseg=="first"){
      for(chr in chr_listy)
      {
        counter_start <- match(chr,chr_listy)
        Sys.sleep(0.1)
        # update progress bar
        setTxtProgressBar(pb, counter_start)
        ind <- unlist(sapply(dfPos.[[chr]][,dfPos.pos.col],function(pos){
          tmp <- which(dfSegm.[[chr]][,dfSegm.start.col] <= pos & dfSegm.[[chr]][,dfSegm.end.col] >= pos)[1]
        }))
        dfPos.[[chr]][,namesColsToAdd] <- dfSegm.[[chr]][ind,colsToAdd]
      }
      close(pb)

    }
    if(multseg=="all"){
      for(chr in chr_listy)
      {
        counter_start <- match(chr,chr_listy)
        Sys.sleep(0.1)
        # update progress bar
        setTxtProgressBar(pb, counter_start)

        ind <- sapply(dfPos.[[chr]][,dfPos.pos.col],function(pos){
          tmp <- which(dfSegm.[[chr]][,dfSegm.start.col] <= pos & dfSegm.[[chr]][,dfSegm.end.col] >= pos)
          tmp
        })
        for(j in 1:length(colsToAdd)){
          dfPos.[[chr]][,namesColsToAdd[j]] <- unlist(lapply(ind,function(z)	paste(dfSegm.[[chr]][z,colsToAdd[j]],collapse=",")))
          dfPos.[[chr]][which(sapply(ind,length)==0),namesColsToAdd[j]] <- NA
        }
      }
      close(pb)
    }
  }
  dfPos <- unsplit(dfPos.,dfPos[,dfPos.chrom.col])
  dfPos
}

#' palimpsest_distCosine
#'
#' Function to calculate cosine distance
#' @param required
#'
#' @return
#' @export
#' @importFrom lsa cosine
#' @examples
palimpsest_distCosine <- function(m)
{
  nsamp <- nrow(m)
  res <- matrix(NA,nrow=nsamp,ncol=nsamp)
  rownames(res) <- colnames(res) <- rownames(m)
  for(i in 1:nsamp)
  {
    for(j in 1:nsamp)
    {
      res[i,j] <- cosine(m[i,],m[j,])
    }
  }
  as.dist(1-res)
}

#' palimpsest_addMutationContextToVcf
#'
#' Function to add mutation contexts to vcf
#' @param vcf vcf data frame containing the SNVs
#' @param Reference_Genome Reference genome (e.g. BSgenome.Hsapiens.UCSC.hg19)
#' @param chrom.col Name of column containing chromosomes (should be in chr1, chr2... format)
#' @param start.col Name of column containing start position
#' @param end.col Name of column containing end position
#' @param ref.col Name of column containing reference base
#' @param alt.col Name of column containing mutated base
#' @param strand.mut.col Name of column containing mutation strand
#' @param strand.gene.col Name of column containing gene strand
#' @param sample.col Name of column containing sample names
#' @param strand.ts.output.col Name of added column with transcribed or non transcribed strand
#' @param substype.output.col Name of added column with substitution type
#' @param ctx3.output.col Name of added column with 3 nucleotide context
#' @param ctx5.output.col Name of added column with 5 nucleotide context
#' @param mutcat3.output.col Name of added column with 3 nucleotide categories (includes substitution type and context)
#' @param mutcat5.output.col Name of added column with 5 nucleotide categories (includes substitution type and context)
#'
#' @return
#' @export
#' @import VariantAnnotation
#' @examples
palimpsest_addMutationContextToVcf <- function(vcf,
                                               Reference_Genome,
                                               chrom.col="chr",
                                               start.col = "start",
                                               end.col="end",
                                               ref.col="ref",
                                               alt.col="alt",
                                               strand.mut.col="strand.mut",
                                               strand.gene.col="strand.gene",
                                               sample.col="samplenames",
                                               strand.ts.output.col="strand.ts",
                                               substype.output.col="substype",
                                               ctx3.output.col="context3",
                                               ctx5.output.col="context5",
                                               mutcat3.output.col="mutcat3",
                                               mutcat5.output.col="mutcat5")
{
    requireNamespace("VariantAnnotation", quietly = TRUE)
    nostrand <- which(is.na(vcf[,strand.gene.col]) | vcf[,strand.gene.col]=="*")
    vcf[nostrand,strand.gene.col] <- "+"
    vr <- VRanges(seqnames=vcf[,chrom.col],ranges=IRanges(start=vcf[,start.col],end=vcf[,end.col]),ref=vcf[,ref.col],alt=vcf[,alt.col],sampleNames=vcf[,sample.col])
    vr@strand <- Rle(vcf[,strand.mut.col])
    vr3 <- palimpsest_addMutationContextToVR(vr,ref=Reference_Genome, k=3,unify = TRUE)
    vr5 <- palimpsest_addMutationContextToVR(vr,ref=Reference_Genome, k=5, unify = TRUE)
    vcf[,strand.mut.col] <- as.character(vr3@strand)
    vcf[,strand.ts.output.col] <- NA
    vcf[which(vcf[,strand.gene.col] == "+" & vcf[,strand.mut.col] =="-"), strand.ts.output.col] <- "ts"
    vcf[which(vcf[,strand.gene.col] == "+" & vcf[,strand.mut.col] =="+"), strand.ts.output.col] <- "nt"
    vcf[which(vcf[,strand.gene.col] == "-" & vcf[,strand.mut.col] =="+"), strand.ts.output.col] <- "ts"
    vcf[which(vcf[,strand.gene.col] == "-" & vcf[,strand.mut.col] =="-"), strand.ts.output.col] <- "nt"
    vcf[nostrand,strand.gene.col] <- NA;vcf[nostrand,strand.ts.output.col] <- NA
    vcf[,substype.output.col] <- as.character(vr3$alteration)
    vcf[,ctx3.output.col] <- as.character(vr3$context)
    vcf[,mutcat3.output.col] <- paste(vcf[,substype.output.col],vcf[,ctx3.output.col],sep="_")
    vcf[which(vcf[,mutcat3.output.col]=="NA_NA"),mutcat3.output.col] <- NA
    vcf[,ctx5.output.col] <- as.character(vr5$context)
    vcf[,mutcat5.output.col] <- paste(vcf[,substype.output.col],vcf[,ctx5.output.col],sep="_")
    vcf[which(vcf[,mutcat5.output.col]=="NA_NA"),mutcat5.output.col] <- NA
    return(vcf)
  }

#' palimpsest_addSVcategoriesToVcf
#'
#' Function to add SV mutation categories to vcf
#' @param sv vcf data frame containing the SVs
#' @param type.col SV type description column name in mutation_data ["DEL","DUP","INV"]
#' @param sample.col Sample column name in sv
#' @param CHROM_1.col Start chromosome column name in sv
#' @param CHROM_2.col End chromosome column name in sv
#' @param POS_1.col Start position column name in sv
#' @param POS_2.col End position column name in sv
#' @param resdir Result directory
#'
#' @return
#' @export
#'
#' @examples
palimpsest_addSVcategoriesToVcf <- function (sv = sv, type.col = "Type", sample.col = "Sample",
                                             CHROM_1.col = "CHROM", CHROM_2.col = "CHROM", POS_1.col = "POS",
                                             POS_2.col = "POS", resdir = resdir)
{
  sv$Size <- 0
  sv$size.cat <- NA
  ind <- which(sv[, type.col] %in% c("DEL", "DUP", "INV"))
  sv[ind, "Size"] <- sv[ind, POS_2.col] - sv[ind, POS_1.col]
  sv[which(sv$Size <= 1000), "size.cat"] <- "a_0-1kb"
  sv[which(sv$Size > 1000 & sv$Size <= 10000), "size.cat"] <- "b_1-10kb"
  sv[which(sv$Size > 10000 & sv$Size <= 1e+05), "size.cat"] <- "c_10-100kb"
  sv[which(sv$Size > 1e+05 & sv$Size <= 1e+06), "size.cat"] <- "d_100kb-1Mb"
  sv[which(sv$Size > 1e+06 & sv$Size <= 1e+07), "size.cat"] <- "e_1-10Mb"
  sv[which(sv$Size > 1e+07), "size.cat"] <- "f_>10Mb"
  sv$Clustered <- 0
  sv$Cluster <- NA
  minbkp <- 10
  for (s in unique(sv[, sample.col])) {
    print(s)
    chr <- pos <- c()
    indsamp <- which(sv[, sample.col] == s)
    for (i in indsamp) {
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
    bed <- bkp[, c(1, 2, 2, 3)]
    indi <- bed2index(bed, sort = TRUE)
    tmp <- cluster.region(x = indi, distance = "1000000",
                          verbose = T)
    tmp[, "chr"] <- data.frame(do.call("rbind", strsplit(as.character(tmp$index), ":", fixed = TRUE)))[1]
    tmp[, c("start", "end")] <- data.frame(do.call("rbind", strsplit(as.character(tmp$index), "-", fixed = TRUE)))[2]
    bkp$cluster <- tmp[, "regionIndex"]
    bkp$clustered <- 0
    tt <- table(bkp$cluster)
    for (clst in names(which(tt >= 10))) {
      bkp[which(bkp$cluster == clst), "clustered"] <- 1
    }
    for (i in which(bkp$clustered == 1)) {
      ind <- which(sv[, sample.col] == s & (paste(sv[,CHROM_1.col], sv[, POS_1.col]) == paste(bkp$chr[i],bkp$pos[i]) | paste(sv[, CHROM_2.col], sv[, POS_2.col]) == paste(bkp$chr[i], bkp$pos[i])))
      sv[ind, "Clustered"] <- 1
      sv[ind, "Cluster"] <- paste(s, bkp[i, "cluster"], sep = ".")
    }
    chr <- pos <- clustered <- c()
    indsamp <- which(sv[, sample.col] == s)
    for (i in indsamp) {
      chr <- c(chr, sv[i, CHROM_1.col], sv[i, CHROM_2.col])
      pos <- c(pos, as.numeric(sv[i, c(POS_1.col, POS_2.col)]))
      clustered <- c(clustered, rep(sv[i, "Clustered"], 2))
    }
    bkp <- data.frame(chr, pos, clustered, stringsAsFactors = F)
    bkp$clustered[is.na(bkp$clustered)] <- 0
  }
  sv$Clustered[is.na(sv$Clustered)] <- 0
  sv$Category <- paste(sv[, type.col], substr(sv$size.cat, 3, nchar(sv$size.cat)), "clust", sv$Clustered, sep = "_")
  sv$Category <- sub("BND_0-1kb", "BND", sv$Category)
  return(sv)
}

#' palimpsest_addMutationContextToVR
#'
#' Function to add mutation context to VRanges object
#' @param vr VRanges object containing the mutations
#' @param ref Reference genome (e.g. BSgenome.Hsapiens.UCSC.hg19)
#' @param k Context length
#' @param check.ref Should consistence with the reference genome be checked?
#' @param check.strand Should mutations on the minus strand be converted (only useful if the input file contains minus strand mutations)
#' @param unify Should substition types be unified to the common 6 mutation types (CA, CG, CT, TA, TC, TG)
#'
#' @return
#' @export
#' @import Biostrings
#' @examples
palimpsest_addMutationContextToVR <- function (vr, ref, k = 3, check.ref = TRUE, check.strand = FALSE, unify = TRUE)
{
  requireNamespace("Biostrings", quietly = TRUE);
  if (any(width(vr)) != 1) stop("SNVs must have width of 1.")
  if (k%%2 != 1) stop("'k' must be odd.")
  mid = (k + 1)/2
  gr = granges(vr)
  strand.mut = strand(gr)
  if (any(strand.mut == "*")) stop("The strand must be explicit, in order to read the correct strand.")
  ranges = resize(gr, k, fix = "center")
  context = getSeq(ref, ranges)
  ref_base = DNAStringSet(ref(vr))
  alt_base = DNAStringSet(alt(vr))
  if (check.ref) {
    ref0 = subseq(context, mid, mid)
    idx_invalid = (ref0 != ref_base)
    if (any(idx_invalid)) warning(sprintf("References do not match in %d cases",sum(idx_invalid)))
  }
  if (check.strand) {
    idx_minus = (strand.mut == "-")
    context[idx_minus] = reverseComplement(context[idx_minus])
    ref_base[idx_complement] = reverseComplement(ref_base[idx_complement])
    alt_base[idx_complement] = reverseComplement(alt_base[idx_complement])
    strand.mut[idx_minus] = "+"
  }
  if (unify) {
    idx_complement = as.character(ref_base) %in% c("A", "G")
    context[idx_complement] = reverseComplement(context[idx_complement])
    ref_base[idx_complement] = reverseComplement(ref_base[idx_complement])
    alt_base[idx_complement] = reverseComplement(alt_base[idx_complement])
    strand.mut[idx_complement] = "-"
  }
  alteration = as.character(xscat(ref_base, alt_base));alteration[idx_invalid] <- NA
  context = as.character(context);context[idx_invalid] <- NA
  vr$alteration = alteration
  vr$context = context
  vr@strand = Rle(strand.mut)
  return(vr)
}

#' palimpsest_makeMutypeMatFromVcf
#'
#' Function to create a matrix in mutation type x sample format with either counts or proportions
#' @param vcf vcf data frame containing the mutations
#' @param sample.col Name of column containing sample names
#' @param mutcat.col Name of column containing mutation categories
#' @param mutypes List of categories to be included in the output
#' @param proportion If TRUE, the output matrix will indicate mutation proportions instead of numbers
#'
#' @return
#' @export
#'
#' @examples
palimpsest_makeMutypeMatFromVcf <- function(vcf,
                                 sample.col="sample",
                                 mutcat.col="mutcat3",
                                 mutypes=c("CA","CG","CT","TA","TC","TG"),
                                 proportion=TRUE)
{
  tmp <- split(vcf,vcf[,sample.col])
  tmp <- lapply(tmp,function(d){
    sapply(mutypes,function(m){sum(d[,mutcat.col]==m,na.rm=T)})
  })
  tmp <- as.matrix(as.data.frame(tmp))
  if(proportion){for(j in 1:ncol(tmp))	tmp[,j] <- tmp[,j]/sum(tmp[,j])}
  tmp
}

#' bed2index
#'
#' convert bed to index object
#' @param x bed file
#' @param sort required
#'
#' @return
#' @export
#'
#' @examples
bed2index <- function(x, sort = TRUE) {
  index <- paste0(
    x[, 1],
    ':',
    x[, 2],
    '-',
    x[, 3]
  );
  if (sort) {
    index <- bedr.sort.region(index);
  }
  index;
}

#' cit.genomOrder
#'
#' @param d data.frame
#' @param chrom chromosome column name
#' @param pos sequence location column name
#' @param startPos the column in d giving the start base pair position
#' @param endPos the column in d giving the end base pair position
#' @param chromNum the name of added column whith chromosome from 1 to 24
#'
#' @return
#' @export
#'
#' @examples
cit.genomOrder <- function(d,
                            chrom     = "chrom",
                            pos       = "pos",
                            startPos  = "start",
                            endPos    = "end" ,
                            chromNum  = "chrNum"
)
{
  d <- factochar(d)
  if(!chromNum %in% names(d))  d <- cit.chromString2num(d,chrom=chrom)

  if (is.null(pos)) {
    pos <- "meanPos"
    if (options("verbose")$verbose)
      warning("missing value for parameter pos -> using default value (meanPos).")
  }

  if(!pos %in% colnames(d))
  {
    if(all(c(startPos,endPos) %in% colnames(d))){d$meanPos <- rowMeans(d[,c(startPos,endPos)],na.rm=T)}else{
      stop("Please set valid pos or startPos and endPos values.")
    }
    pos <- "meanPos"
  }
  ord <- order(d[,chromNum], as.numeric(d[,pos]), na.last = TRUE)
  d <- d[ord,]

  d
}


#' cit.pangenomCoord
#'
#' @param x data.frame which columns include genomic position information
#' @param chrom the column in x giving the chromosome position (in 1:22,X,Y)
#' @param pos the column in x giving the base pair position   (optional)
#' @param startPos the column in x giving the start base pair position
#' @param endPos the column in x giving the end base pair position
#' @param cyto cytoband object (optional)
#' @param absPos name of added column giving absolute pangenomic order
#' @param chromNum the column in x giving the chromosome position with X = 23 & Y =24
#'
#' @return
#' @export
#'
#' @examples
cit.pangenomCoord <- function( x,
                                chrom="chrom",
                                pos="meanPos",
                                startPos="start",
                                endPos="end",
                                cyto=NULL,
                                absPos="absPos",
                                chromNum="chrNum"
)
{
  if(is.null(cyto))  stop("Cytoband object required")

  x <- cit.genomOrder(d=x,
                       chrom     = chrom,
                       pos       = pos,
                       startPos  = startPos,
                       endPos    = endPos ,
                       chromNum  = chromNum)

  if(absPos %in% colnames(x))  return(x)

  if(!endPos %in% colnames(x)){
    endPos <- "endPos"
    if(pos %in% colnames(x)){
      x[,endPos] <- x[,pos]
      if(options("verbose")$verbose) warning("uncorrect value for parameter endPos -> using pos value.")
    }else{
      stop("missing value for pos or endPos.")
    }
  }

  if(is.null(absPos)){
    absPos <- "absPos"
    if(options("verbose")$verbose) warning("missing value for parameter absPos -> using default value.")
  }

  if(!absPos %in% names(x)){
    x$absPos <- as.numeric(x[,pos])
    if(pos %in% colnames(x)){
      x[,absPos] <- as.numeric(x[,pos])
      if(options("verbose")$verbose) warning("uncorrect value for parameter absPos -> using pos value.")
    }else{
      stop("missing value for absPos or pos.")
    }
  }

  if(!is.null(cyto)){
    if(!"ChrNumeric" %in% names(cyto) ) cyto <-  cit.chromString2num(cyto,chrom="Chromosome",chromNum="ChrNumeric")

    allchrnum <- sort(unique(c(cyto$ChrNumeric,x[,chromNum])))

    maxByChr <- sapply(allchrnum,function(z){
      m <- 0
      w <- which(cyto$"ChrNumeric"==z)
      if(length(w)>0) m <- max(cyto[w,"End"],na.rm=TRUE)
      w <- which(x[,chromNum]==z)
      if(length(w)>0){
        if( which.max(c(m,max(x[w,endPos])))==2 )
          warning("In cit.pangenomCoord, some positions in your data are superior to max cytoband in cyto.\ncyto object used does not seem to be the right cytoband object version.\n")
        m <- max(c(m,max(x[w,endPos],na.rm=TRUE)))
      }
      m})

    maxByChr <- cumsum(as.numeric(maxByChr))

  }else{
    maxByChr <- cumsum(as.numeric(sapply(split(x[,endPos],x[,chromNum]),max)))
  }


  for(i in setdiff(unique(x[,chromNum]),1)){
    w <- which(x[,chromNum] == i)
    x[w,absPos] <-  x[w,absPos] + maxByChr[i-1]
  }
  w <- which(is.na(x[,chromNum]))
  if(length(w)>0) x[w,absPos] <- NA
  x
}

#' cit.pangenomPlot
#'
#' pangenomic plot
#' @param d data.frame which columns include genomic position information and the wanted y axis
#' @param ycol the column in d to be used as y axis
#' @param chrom the column in  d giving the chromosome position
#' @param pos the column in d giving the base pair position   (optional)
#' @param startPos the column in d giving the start base pair position
#' @param endPos the column in  d giving the end base pair position
#' @param colorscol (optional) column in d giving the color for each point from ycol or a unique color name or number (default black)
#' @param cyto cytoband object (optional)
#' @param absPos  the column in d giving the absolute base pair position (will be calculated if not present in d)
#' @param rectangleOption boolean : if {TRUE} rectangles are plotted instead of points (default) - NB : for points, several 'looks' can be obtained (ex using parameter type , pch,...)
#' @param plotCentro boolean : should centromeric delimitation be plotted ?
#' @param plotCentro.color color of the centromeric delimitation
#' @param plotChrom boolean : should chromosome delimitation be plotted ?
#' @param plotnew boolean : should the a new graph be plotted (default={TRUE}), or should something be added to the current graph (->{plotnew=FALSE})
#' @param plotyaxis boolean : should an Y axis be plotted ? (considered only if {plotnew = TRUE})
#' @param yaxmark  y axis 'at' parameter (considered only if {plotyaxis = TRUE})
#' @param yaxlab y axis 'labels' parameter (considered only if {plotyaxis = TRUE})
#' @param chromToPlot a vector of the chromosomes to be plotted (default c(1:22,"X","Y"))

#'
#' @return
#' @export
#'
#' @examples
cit.pangenomPlot <- function(d=NULL,
                              ycol=NULL,
                              chrom="chrom",
                              pos="meanPos",
                              startPos="start",
                              endPos="end",
                              colorscol=NULL,
                              cyto=NULL,
                              absPos="absPos",
                              rectangleOption=FALSE,
                              plotCentro = TRUE,
                              plotCentro.color = "lightgrey",
                              plotChrom =TRUE,
                              plotnew=TRUE,
                              plotyaxis=TRUE,
                              plotxaxis=TRUE,
                              yaxmark=NULL,
                              yaxlab=NULL,
                              chromToPlot=c(1:22,"X","Y"),
                              xlim=NULL,
                              CexAxis=0.8,
                              ...)
{
  if(! ycol %in% names(d))  stop("ycol ", ycol, " not found in d.\n")
  if(is.null(cyto))  stop("Cytoband object required")

  d <- cit.pangenomCoord(d,chrom=chrom,pos=pos,startPos=startPos,endPos=endPos,cyto=cyto)
  y <- d[,ycol]

  if(is.null(colorscol) ){colorY <-"black"}else{colorY <- d[,colorscol]}

  if(!is.null(chromToPlot))	cyto <- cyto[which(cyto$Chromosome %in% chromToPlot),] # added to rm X and Y chromosome from graphics
  cyto <- cit.pangenomCoord (cyto, chrom="ChrNumeric", pos="End", startPos="Start", endPos="End", cyto=cyto, absPos="absPos")
  if(is.null(xlim)) {xlim <- c(0,max(cyto$absPos,na.rm=T))}
  delta <- diff(xlim)*0.02;xlim <- xlim+c(-delta,delta)
  delta <- diff(xlim)*0.04/0.92;xlim <- xlim+c(delta,-delta)

  if(plotnew){
    plot( d[,absPos],y,axes=FALSE,col=colorY,xlim=xlim, ...)

    box()

    if(plotyaxis){
      okk <- FALSE
      if(is.null(yaxmark)){
        yaxmark <- pretty(y,n=10)
        okk <- TRUE
      }
      if(is.null(yaxlab)) yaxlab <- yaxmark
      axis(2, #1=below, 2=left, 3=above and 4=right
           at = yaxmark,
           labels = yaxlab, cex.axis=CexAxis, las=2)
    }

  }else{
    points( d[,absPos],y, col=colorY, ...)
  }

  wat <- NULL
  if(!is.null(cyto)){
    w <- which(cyto$"Centro"==1)
    wat <- sapply(split(cyto[w,"absPos"],cyto[w,"ChrNumeric"]),mean)
    if(length(w) == 0) wat <- sapply(split(cyto[,"absPos"],cyto[,"ChrNumeric"]),mean)

    if(plotCentro & length(w) != 0)abline(v=wat,lty=3,col=plotCentro.color)

    if(plotChrom ){

      if(plotxaxis){
        axis(1, #1=below, 2=left, 3=above and 4=right
             at = wat,
             labels = c(1:22,"X","Y")[c(1:22,"X","Y") %in% unique(cyto$Chromosome)],#c(1:22,"X","Y"),
             cex.axis=CexAxis, las=2)
      }

      wat <-  sapply(split(cyto[,"absPos"],cyto$"ChrNumeric"),max)
      abline(v=wat[-length(wat)],lty=1)
    }
  }
  d
}


#' factochar
#'
#' Function to convert factors factors in a data frame to characters
#' @param d Data frame to be converted
#'
#' @return
#' @export
#'
#' @examples
factochar <- function(d) {
  for(i in 1:ncol(d)) if(is.factor(d[,i])) d[,i] <- as.character(d[,i])
  d
}

#' cit.chromString2num
#'
#' convert chromosomes X, Y and MT into numerics 23, 24 and 25
#' @param d data.frame which columns include genomic position information
#' @param chrom the column in \code{d} giving the chromosome position
#'
#' @return
#' @export
#'
#' @examples
cit.chromString2num <- function (d,
                                  chrom = "chrom",
                                  chromNum  = "chrNum"
)
{
  d[, chromNum] <- as.character(sub("chr","",d[, chrom]))

  maxchr <- suppressWarnings( max(as.numeric(d[,chromNum]),na.rm=TRUE) )

  w <- which(d[, chromNum] == "X")
  if (length(w) > 0)d[w, chromNum] <- maxchr+1

  w <- which(d[, chromNum] == "Y")
  if (length(w) > 0)d[w, chromNum] <- maxchr+2

  w <- which(d[, chromNum] %in% c("M","MT"))
  if (length(w) > 0)d[w, chromNum] <- maxchr+3

  w <- which(!(d[, chromNum] %in% as.character(1:25)))
  if (length(w) > 0) d[w, chromNum] <- NA

  d[, chromNum] <- as.numeric(d[, chromNum])
  d
}

#' cit.density
#'
#' wrapping of the function 'density' adding down and top points to the result
#' @param x a numeric vector
#' @param doplot boolean
#' @param pc
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
cit.density <- function(x,doplot=FALSE,pc=.05,...){
  dx <- density(x,na.rm=TRUE,...)
  ymax <- diff(range(dx$y))
  n <- length(dx$y)
  croissance <- as.numeric((dx$y[-1] - dx$y[-n])>0)
  wB <- 1+which(diff(croissance)== 1)
  wH <- 1+which(diff(croissance)== -1)
  pointsBas   <- dx$x[wB]
  pointsHauts <- dx$x[intersect(wH,which(dx$y>pc*ymax))]

  if(length(pointsHauts)>0)pointsBas   <- sapply(split(pointsBas,cit.discretize(pointsBas,pointsHauts)),median)

  if(doplot){
    plot(dx,...)
    abline(v=pointsBas,lty=3,col="blue")
    abline(v=pointsHauts,lty=3,col="red")
  }
  L <- c(dx,list("down"=pointsBas,"top"=pointsHauts))
  attr(L,"class") <- "density"
  L
}

#' cit.peaks
#'
#' wrapping of the function 'cit.density' to control the maximum number of peaks
#' @param x a numeric vector
#' @param percentHighestPeak
#' @param maxNbPeaks
#' @param minDeltaBetweenPeaks
#' @param deltaApproach
#' @param doplot boolean
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
cit.peaks <- function(x,
                       percentHighestPeak=.2,
                       maxNbPeaks=NULL,
                       minDeltaBetweenPeaks=.03,
                       deltaApproach=1,
                       doplot=FALSE,
                       ...){
  if(is.na( minDeltaBetweenPeaks ))minDeltaBetweenPeaks <- NULL
  v <- cit.density(x,pc =percentHighestPeak)
  m <- v$y[which(v$x %in% v$top)]
  w <- which(m/max(m) > percentHighestPeak)
  xw <- v$top[w]
  which.eq <- function(z){ order(abs(v$x-z))[1]}
  if(deltaApproach==1 & !is.null( minDeltaBetweenPeaks ) & length(xw)>1){
    for(i in 1:(length(xw)-1)){
      if(xw[i+1]-xw[i] <  minDeltaBetweenPeaks) {
        xw[i:(i+1)] <- xw[c(i:(i+1))][which.max(m[i:(i+1)])]
      }
    }

    xw <- unique(xw)

  }
  if(!is.null(maxNbPeaks)){
    if(length(xw)> maxNbPeaks)
      for(i in 1:(length(xw)-maxNbPeaks)){
        oneamong <- xw[c(0,1)+which.min(diff(xw))]
        w1 <-  which.eq(oneamong[1])
        w2 <-  which.eq(oneamong[2])
        out <- oneamong[which.min(v$y[c(w1,w2)])]
        xw <- setdiff(xw,out)

      }
  }
  umin <- NULL
  if(length(xw)>1){
    for(i in 1:(length(xw)-1)){
      possiblemin <- v$down[which(v$down >= xw[i] & v$down <= xw[i+1])]
      if(length(possiblemin)>1){
        w<- sapply(possiblemin,which.eq)
        possiblemin <- possiblemin[which.min(v$y[w])]
      }
      umin <- c(umin,possiblemin)
    }

  }
  L <- list("x abciss big peaks"=xw,
            "x abciss inter-peaks"=umin,
            "nb big peaks"=length(xw))

  peaks <- L[[1]]
  interpeaks <- L[[2]]
  temp <- as.data.frame(t(sapply(peaks,function(pic) {
    wleft <- which(interpeaks < pic)
    if(length(wleft)>0){
      wleft <- max(interpeaks[wleft])
    }else{
      wleft <- min(x,na.rm=TRUE)
    }
    wright <- which(interpeaks > pic)
    if(length(wright)>0){
      wright <- min(interpeaks[wright])
    }else{
      wright <- max(x,na.rm=TRUE)
    }
    c(pic,wleft,wright,length(which(x>=wleft & x<=wright)))
  }) ))
  names(temp) <- c("peak","left born","right born","size")
  L$"peaks x size" <- temp


  if(deltaApproach==2 &!is.null( minDeltaBetweenPeaks ) & L$'nb big peaks'>1){
    Lini <- L
    names(Lini) <- paste("initial",names(L))

    while(any( diff(L$"x abciss big peaks") < minDeltaBetweenPeaks) & L$"nb big peaks" >1){

      w <- which.min(abs(diff(L$"x abciss big peaks")) )


      lb <- min(L$"peaks x size"[w:(w+1),"left born"])
      rb <- max(L$"peaks x size"[w:(w+1),"right born"])
      si <- sum(L$"peaks x size"[w:(w+1),"size"])
      pe <- median(x[which(x>=lb & x<=rb)])

      L$"peaks x size"[w,] <- c(pe,lb,rb,si)
      L$"peaks x size" <- L$"peaks x size"[-(w+1),]

      L$"nb big peaks" <- L$"nb big peaks" - 1
      L$"x abciss big peaks"[w] <- pe
      L$"x abciss big peaks" <- L$"x abciss big peaks"[-(w+1)]
      L$"x abciss inter-peaks" <- L$"x abciss inter-peaks"[-(w+1)]
    }

    L <- c(Lini,L)
  }

  if(doplot){
    plot(v,...)
    abline(v=L$"x abciss big peaks",col="red")
    if(length(L$"x abciss inter-peaks")>0)abline(v=L$"x abciss inter-peaks",col="green",lty=3)
  }

  L
}

#' cit.discretize
#'
#' Discretize a continous variable by specified \code{lim} cut-off(s)
#' @param x a vector
#' @param lim cut-off(s)
#' @param quant if TRUE (default FALSE) the \code{x} is discretize by quantile and \code{lim} is considered as cut-off(s) for quantile, ie 0<lim<1
#' @param addlevels add character levels indicating the cut-offs (i.e. for un cut-off iqq levels=c("<iqq",">=iqq") )
#'
#' @return
#' @export
#'
#' @examples
cit.discretize <- function(x, lim, quant = FALSE, addlevels = FALSE){

  lim <- sort(lim)

  if( quant & ( any(lim>1) | any(lim<0) ) )
    stop("lim must be [0;1] as quant=T\n")

  res <- rep(NA, length(x))

  if(quant) lim  <- quantile(x, probs=lim, na.rm=TRUE)
  n <- length(lim)
  for(i in n:1) res[which(x<lim[i])] <- i
  res[which(x>=lim[n])] <- n+1

  if (addlevels) {
    res <- as.factor(res)

    if(quant)
      lim <- gsub(" ", "", prettyNum(lim, format="g", digits=1))

    if(length(lim)==1)
      lev <- c(paste("<", lim[1], sep = ""), paste(">=", lim[1], sep = ""))
    else
      lev <- c(paste("<", lim[1], sep = ""), paste("[",lim[-length(lim)],";",lim[-1],"[",sep=""), paste(">=", lim[length(lim)], sep = ""))
    levels(res) <- lev[as.numeric(levels(res))]
  }
  res
}


#' factoall
#'
#' Function to convert factors in a data frame to the appropriate format (character, numeric...)
#' @param d Data frame to be converted
#' @param ncmax Maximum number of characters for numeric fields
#'
#' @return
#' @export
#'
#' @examples
factoall <- function (d, ncmax = 10)
{
  n <- ncol(d)
  for (i in 1:n) {
    if (is.factor(d[, i])) {
      d[, i] <- as.character(d[, i])
      na <- which(is.na(d[, i]))
      num <- suppressWarnings(as.numeric(d[, i]))
      nanum <- which(is.na(num))
      if (length(nanum) == length(na)) {
        d[, i] <- num
      }
    }
  }
  d
}

