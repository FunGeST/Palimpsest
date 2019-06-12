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
#' @export
#' @import gtools

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
      if (all(is.na(ind)) == TRUE)
        dfPos.[[chr]][, namesColsToAdd] <-
        NA
      if (all(is.na(ind)) != TRUE) {
        dfPos.[[chr]][, namesColsToAdd] <- dfSegm.[[chr]][ind,
                                                          colsToAdd]
      }
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
#' @param required m 
#'
#' @export
#' @importFrom lsa cosine

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
#' @export
#' @import VariantAnnotation

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
    vr@strand <- Rle(strand(vcf[, strand.mut.col]))
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

#' palimpsest_makeMutypeMatFromVcf
#'
#' Function to create a matrix in mutation type x sample format with either counts or proportions
#' @param vcf vcf data frame containing the mutations
#' @param sample.col Name of column containing sample names
#' @param mutcat.col Name of column containing mutation categories
#' @param mutypes List of categories to be included in the output
#' @param proportion If TRUE, the output matrix will indicate mutation proportions instead of numbers
#'
#' @export

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
#' @export

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
#' @export

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
#' @export

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
#' @export

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
#' @export

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
#' @export

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
#' @param pc ?
#' @param ...
#'
#' @export

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
#' @param percentHighestPeak percentHighestPeak
#' @param maxNbPeaks maxNbPeaks
#' @param minDeltaBetweenPeaks minDeltaBetweenPeaks
#' @param deltaApproach deltaApproach
#' @param doplot boolean
#' @param ...
#' 
#' @export

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
#' @export

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
#' @export

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


#' "%!in%" 
#'
#' The opposite of "%in%"
#' @keywords Signatures
#' @export
#' @examples
#' days <- c("monday","tuesday","wednesday","thursday","friday","saturday","sunday")
#' weekdays <- days[days %!in% c("saturday", "sunday")]

'%!in%' <- function(x,y)!('%in%'(x,y))



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



load2object <- function (filename) 
{
  if (file.exists(filename)) 
    return(eval(parse(text = load(filename))))
  cat(paste("error! the file ", filename, 
            " was not found :("))
  NULL
}




#' palimpsest_distCosine
#'
#' Function to calculate cosine distance
#' @param m input
#'
#' @export
#' @importFrom lsa cosine

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



compare_results <- function(reference_sigs = NA, extraction_1 = NA, extraction_2 = NULL, extraction_1_name ="Palimp", extraction_2_name = NA, lower_threshold = 0.6, upper_threshold = 0.9){

  #ifelse(!is.na(extraction_2),print("good to go"),stop("function currently needs 2 denovo sets to work"))
  if(!missing(extraction_2)){
    rownames(extraction_1) <- c(paste0(extraction_1_name,"_",c(1:nrow(extraction_1))))
    rownames(extraction_2) <- c(paste0(extraction_2_name,"_",c(1:nrow(extraction_2))))
    
    refs <- rownames(reference_sigs)
    denovs1 <- rownames(extraction_1)
    denovs2 <- rownames(extraction_2)
    
    longest <- max(c(length(denovs1),length(denovs2)))
    
    res <- as.data.frame(matrix(nrow = nrow(reference_sigs)+longest,ncol = 0)) %>% 
      mutate(References = NA, denovo1 = NA, ref_denovo1_cos = NA, denovo2 = NA, ref_denovo2_cos = NA, denovo1_denovo2_cos = NA, keep_version = NA)
  
    res$References[1:nrow(reference_sigs)] <- c(rownames(reference_sigs))
                                                
                                                
    # perfrom cosine similarity analysis
    mutmat1 <- t(rbind(extraction_1, reference_sigs))
    m1 <- as.matrix(palimpsest_distCosine(t(mutmat1)))
    m1 <- 1-m1
    mutmat2 <- t(rbind(extraction_2, reference_sigs))
    m2 <- as.matrix(palimpsest_distCosine(t(mutmat2)))
    m2 <- 1-m2
    mutmat3 <- t(rbind(extraction_1, extraction_2))
    m3 <- as.matrix(palimpsest_distCosine(t(mutmat3)))
    m3 <- 1-m3
    
    # make results table
    for(i in 1:nrow(reference_sigs)){
      sig <- res$References[i]
      if(max(m1[sig,colnames(m1) %in% denovs1]) > lower_threshold){
        res$ref_denovo1_cos[i] <- max(m1[sig,colnames(m1) %in% denovs1])
        res$denovo1[i] <- rownames(m1)[m1[sig,]==res$ref_denovo1_cos[i]]
      }
      if( max(m2[sig,colnames(m2) %in% denovs2]) > lower_threshold){
        res$ref_denovo2_cos[i] <- max(m2[sig,colnames(m2) %in% denovs2])
        res$denovo2[i] <- rownames(m2)[m2[sig,]==res$ref_denovo2_cos[i]]
      }
    }
    
    manque1 <- denovs1[denovs1 %!in% res$denovo1]
    manque2 <- denovs2[denovs2 %!in% res$denovo2]
    
    if(length(manque1) > 0){
      for(i in 1:length(manque1)){
          res$denovo1[length(refs)+i] <- manque1[i]
          sig <- manque1[i]
          res$denovo1_denovo2_cos[length(refs)+i] <- max(m3[sig,colnames(m3) %in% denovs2])
          res$denovo2[length(refs)+i] <- rownames(m3)[m3[sig,]==res$denovo1_denovo2_cos[length(refs)+i]]
      }
    }
    for(i in 1:(length(refs)+length(manque1))){
      if(!is.na(res$denovo1[i]) & !is.na(res$denovo2[i])) res$denovo1_denovo2_cos[i]<- m3[res$denovo1[i],res$denovo2[i]]
    }
    
    manque2x <- manque2[manque2 %!in% res$denovo2]
    if(length(manque1) > 0){
      for(i in 1:length(manque2x)){
        if(length(manque2x)==0) break
        res$denovo2[length(refs)+i+length(manque1)] <- manque2x[i]
      }
    }
    
    for(i in (length(refs)+1):nrow(res)){
      sig1 <- res$denovo1[i]; sig2 <- res$denovo2[i]
      if(!is.na(sig1)){
        if(max(m1[sig1,colnames(m1) %in% refs]) > lower_threshold){
          res$ref_denovo1_cos[i] <- max(m1[sig1,colnames(m1) %in% refs])
          res$References[i] <- colnames(m1)[m1[sig1,]==res$ref_denovo1_cos[i]]
          if(!is.na(sig2)) res$ref_denovo2_cos[i] <- m2[sig2,res$References[i]]
        }
      }
      if(!is.na(sig2)){
        if(max(m2[sig2,colnames(m2) %in% refs]) > lower_threshold){
          res$ref_denovo2_cos[i] <- max(m2[sig2,colnames(m2) %in% refs])
          res$References[i] <- colnames(m2)[m2[sig2,]==res$ref_denovo2_cos[i]]
          if(!is.na(sig1)) res$ref_denovo1_cos[i] <- m1[sig1,res$References[i]]
        }
      }
    }
    
    # add suggestion annotations
    res_f = res %>% 
      filter(rowSums(is.na(.)) < ncol(.)) %>% 
      mutate(keep = ifelse((ref_denovo1_cos > lower_threshold & ref_denovo2_cos > lower_threshold)|ref_denovo1_cos > upper_threshold | ref_denovo2_cos > upper_threshold,"Yes",NA)) #mark those to keep
    
    res_f$keep[res_f$denovo1_denovo2_cos > upper_threshold] <- "Yes"
    res_f$keep[is.na(res_f$keep)] <- "No"
    
    res_f = res_f %>%
      mutate(keep_version = ifelse(keep == "Yes" & References %in% refs & ref_denovo1_cos > upper_threshold | ref_denovo2_cos > upper_threshold , "Keep_Reference",keep_version)) #highlight for which ref is defo best
      
    res_f = res_f %>% 
      mutate(keep_version = ifelse(keep == "Yes" & is.na(keep_version),"Check_denovo_results",keep_version)) #highlight for which user should check
    res_f$keep_version[res_f$keep == "No"] <- "Discard"
    duped <- unique(res_f$References[duplicated(res_f$References) & !is.na(res_f$References)])
    res_f$keep_version[res_f$References %in% duped] <- "Check_denovo_results"
    
    colnames(res_f) <- c("Ref_Signature",paste0(extraction_1_name,"_Equivalent"),paste0("Ref_",extraction_1_name,"_cosine_score"),paste0(extraction_2_name,"_Equivalent"),
                         paste0("Ref_",extraction_2_name,"_cosine_score"), paste0(extraction_1_name,"_",extraction_2_name,"_cosine_score"),"Suggested_Action","Keep")
    
    res_f <- res_f[,c(1:7)]
    res_f
  }
  else{
    rownames(extraction_1) <- c(paste0(extraction_1_name,"_",c(1:nrow(extraction_1))))

  refs <- rownames(reference_sigs)
  denovs1 <- rownames(extraction_1)

  
  res <- as.data.frame(matrix(nrow = nrow(reference_sigs)+length(denovs1),ncol = 0)) %>% 
    mutate(References = NA, denovo1 = NA, ref_denovo1_cos = NA, keep_version = NA)
  
  res$References[1:nrow(reference_sigs)] <- c(rownames(reference_sigs))
  
  # perfrom cosine similarity analysis
  mutmat1 <- t(rbind(extraction_1, reference_sigs))
  m1 <- as.matrix(palimpsest_distCosine(t(mutmat1)))
  m1 <- 1-m1

  
  # make results table
  for(i in 1:nrow(reference_sigs)){
    sig <- res$References[i]
    if(max(m1[sig,colnames(m1) %in% denovs1]) > lower_threshold){
      res$ref_denovo1_cos[i] <- max(m1[sig,colnames(m1) %in% denovs1])
      res$denovo1[i] <- rownames(m1)[m1[sig,]==res$ref_denovo1_cos[i]]
    }
  }
  
  manque1 <- denovs1[denovs1 %!in% res$denovo1]
  
  if(length(manque1) > 0){
    for(i in 1:length(manque1)){
      res$denovo1[length(refs)+i] <- manque1[i]
    }
  }

  
  
  # add suggestion annotations
  res_f = res %>% 
    filter(rowSums(is.na(.)) < ncol(.)) %>% 
    mutate(keep = ifelse((ref_denovo1_cos > lower_threshold)|ref_denovo1_cos > upper_threshold,"Yes",NA)) #mark those to keep
  
  res_f$keep[is.na(res_f$keep)] <- "No"
  
  res_f = res_f %>%
    mutate(keep_version = ifelse(keep == "Yes" & References %in% refs & ref_denovo1_cos > upper_threshold, "Keep_Reference",keep_version)) #highlight for which ref is defo best
  
  res_f = res_f %>% 
    mutate(keep_version = ifelse(keep == "Yes" & is.na(keep_version),"Check_denovo_results",keep_version)) #highlight for which user should check
  res_f$keep_version[res_f$keep == "No"] <- "Discard"
  duped <- unique(res_f$References[duplicated(res_f$References) & !is.na(res_f$References)])
  res_f$keep_version[res_f$References %in% duped] <- "Check_denovo_results"
  
  colnames(res_f) <- c("Ref_Signature",paste0(extraction_1_name,"_Equivalent"),paste0("Ref_",extraction_1_name,"_cosine_score"),
                       "Suggested_Action","Keep")
  
  res_f <- res_f[,c(1:4)]
  res_f
    
}
}


