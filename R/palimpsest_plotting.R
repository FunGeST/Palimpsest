#' plotSubstypeNumberPerSample
#'
#' @param vcf vcf data frame containing the mutations
#' @param sample.col Sample column name in vcf
#' @param substype.col Substitution types column name in vcf
#' @param plot.file File to which the plot should be saved (if left NULL, plot will be sent directly to the active window)
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plotSubstypeNumberPerSample <- function(vcf,
                                        sample.col="sample",
                                        substype.col="substype",
                                        plot.file=NULL, ...)
{
  mt <- c("CA","CG","CT","TA","TC","TG")
  tt = table(vcf[,substype.col],vcf[,sample.col])[mt,]
  nmut <- apply(tt,2,sum)
  if(!is.null(plot.file)){pdf(plot.file,width=10,height=5)}
  barplot(tt[,order(nmut,decreasing=T)],las=2,...,col=c("skyblue3","black","red","grey","green","pink"),legend.text=mt,args.legend=list(x="topright",bty ="n"),ylab="# of mutations")
  if(!is.null(plot.file)){dev.off()}
}


#' plot6mutationSpectrumFromVcf
#'
#' Function to plot the proportions of each of the 6 substitution types from a vcf file (containing one or a series of samples)
#' @param vcf vcf data frame containing the mutations
#' @param sample.col Sample column name in vcf
#' @param substype.col Substitution types column name in vcf
#' @param averageProp If TRUE, the proportions returned will be the average of sample-wise proportions, otherwise the general proportion over the vcf is used (induces slight changes)
#' @param plot.file File to which the plot should be saved (if left NULL, plot will be sent directly to the active window)
#'
#' @return
#' @export
#'
#' @examples
plot6mutationSpectrumFromVcf <- function(vcf,
                                         sample.col="sample",
                                         substype.col="substype",
                                         averageProp=FALSE,
                                         plot.file=NULL)
{
  mt <- c("CA","CG","CT","TA","TC","TG")
  nsamp <- length(unique(vcf[,sample.col]))
  if(averageProp & nsamp > 1){
    tmp <- makeMutypeMatFromVcf(vcf,sample.col="CHCID",mutcat.col=substype.col,mutypes=mt)
    freq <- apply(tmp,1,mean)
  }else{freq <- sapply(mt,function(z){mean(vcf[,substype.col]==z,na.rm=T)})}
  if(!is.null(plot.file)){pdf(plot.file,width=5)}
  col6 <- c("skyblue3","black","red","grey","green","pink")
  bp <- barplot(freq,col=col6,border=col6,las=2,width=1,space=1,las=1,ylab="Percentage of mutations")
  if(!is.null(plot.file)){dev.off()}
}

#' plot96mutationSpectrumFromVcf
#'
#' Function to plot the proportions of each of the 96 mutation types from a vcf file (containing one or a series of samples)
#' @param vcf vcf data frame containing the mutations
#' @param sample.col Sample column name in vcf
#' @param mutcat3.col Mutation categories column name in vcf
#' @param ymax Upper value for the y-axis (if left NULL it will be computed automatically)
#' @param averageProp If TRUE, the proportions returned will be the average of sample-wise proportions, otherwise the general proportion over the vcf is used (induces slight changes)
#' @param plot.file File to which the plot should be saved (if left NULL, plot will be sent directly to the active window)
#'
#' @return
#' @export
#'
#' @examples
plot96mutationSpectrumFromVcf <- function(vcf,
                                          sample.col="sample",
                                          mutcat3.col="mutcat3",
                                          ymax=NULL,
                                          averageProp=FALSE,
                                          plot.file=NULL
)
{
  bases <- c("A","C","G","T")
  ctxt16 <- paste(rep(bases,each=4),rep(bases,4),sep=".")
  mt <- c("CA","CG","CT","TA","TC","TG")
  types96 <- paste(rep(mt,each=16),rep(ctxt16,6),sep="_")
  types96 <- sapply(types96,function(z){sub("\\.",substr(z,1,1),z)})
  context <- substr(types96,4,6)
  nsamp <- length(unique(vcf[,sample.col]))
  if(averageProp & nsamp > 1){
    tmp <- makeMutypeMatFromVcf(vcf,sample.col="CHCID",mutcat.col="mutcat3",mutypes=types96)
    freq <- apply(tmp,1,mean)
  }else{freq <- sapply(types96,function(z){mean(vcf[,mutcat3.col]==z,na.rm=T)})}
  if(!is.null(plot.file)){pdf(plot.file,width=24,height=5)}
  col96 <- c(rep("skyblue3",16),rep("black",16),rep("red",16),rep("grey",16),rep("green",16),rep("pink",16))
  labs <-c(rep("C>A",16),rep("C>G",16),rep("C>T",16),rep("T>A",16),rep("T>C",16),rep("T>G",16))
  if(is.null(ymax)){ymax <- ceiling(max(freq)*100)/100}
  bp <- barplot(freq,col=col96,border=col96,las=2,width=1,space=1,yaxt="n",xaxt="n",ylim=c(0,ymax*1.2))
  title(ylab="Percentage of mutations",mgp=c(1,1,0), cex.lab=1.6)
  axis(1,at=bp,labels=context,pos=0,las=2,cex.axis=1.5,tick=T,cex.axis=1)
  axis(2,at=round(seq(0,ceiling(ymax*100),length.out=3),digits=1)/100,pos=0,las=1,cex.axis=1.5)
  for(i in seq(1,81,by=16)){
    rect(bp[i],par()$usr[4],bp[i+15],par()$usr[4]-0.05*diff(par()$usr[3:4]),col=col96[i],border=col96[i])
    text((bp[i]+bp[i+15])/2,par()$usr[4]+0.09*diff(par()$usr[3:4]),labels=labs[i], xpd= TRUE, cex = 2)
  }
  if(!is.null(plot.file)){dev.off()}
}

#' plotStrandBias6types
#'
#' Function to calculate and plot the transcriptional strand bias over the 6 substitution types
#' @param vcf vcf data frame containing the mutations
#' @param substype.col Substitution types data frame containing the mutations
#' @param strand.ts.col Name of column containing indicating whether mutations are on the transcribed or untranscribed strand
#' @param plot.file File to which the plot should be saved (if left NULL, plot will be sent directly to the active window)
#'
#' @return
#' @export
#'
#' @examples
plotStrandBias6types <- function(vcf,
                                 substype.col="substype",
                                 strand.ts.col="strand.ts",
                                 plot.file=NULL)
{
  mt <- c("CA","CG","CT","TA","TC","TG")
  nt <- sapply(mt,function(m){sum(vcf[which(vcf[,substype.col]==m),strand.ts.col]=="nt",na.rm=T)})
  ts <- sapply(mt,function(m){sum(vcf[which(vcf[,substype.col]==m),strand.ts.col]=="ts",na.rm=T)})
  tt <- rbind(ts,nt)
  pval <- as.numeric(apply(tt,2,function(z){try(prop.test(z[1],sum(z),p=0.5)$p.value,silent=TRUE)}))
  pval.lab <- rep("ns",length(pval));pval.lab[which(pval < 0.05)] <- "*";pval.lab[which(pval < 0.01)] <- "**";pval.lab[which(pval < 0.001)] <- "***"
  ymax <- max(tt)*1.2
  if(!is.null(plot.file)){pdf(plot.file)}
  bp <- barplot(tt,beside= TRUE,col= c("royalblue","red"),border=c("royalblue","red"),legend.text= c("Transcribed","Untranscribed"),args.legend=list(x="topright",bty = "n",inset=c(-0.04,-.08)),las=1,ylab="Number of mutations",xlab="Mutation type",ylim=c(0,ymax))
  text(x=apply(bp,2,mean),y=apply(tt,2,max)+0.03*diff(par()$usr[3:4]),pval.lab)
  abline(h=0)
  if(!is.null(plot.file)){dev.off()}
  return(rbind(tt,pval))
}

#' plotStrandBias96types
#'
#' Function to calculate and plot the transcriptional strand bias over the 96 mutation types
#' @param vcf data frame containing the mutations
#' @param mutcat3.col Mutation categories column name in vcf
#' @param strand.ts.col Name of column containing indicating whether mutations are on the transcribed or untranscribed strand
#' @param plot.file File to which the plot should be saved (if left NULL, plot will be sent directly to the active window)
#'
#' @return
#' @export
#'
#' @examples
plotStrandBias96types <- function(vcf,
                                  mutcat3.col="mutcat3",
                                  strand.ts.col="strand.ts",
                                  plot.file=NULL)
{
  bases <- c("A","C","G","T")
  ctxt16 <- paste(rep(bases,each=4),rep(bases,4),sep=".")
  mt <- c("CA","CG","CT","TA","TC","TG")
  types96 <- paste(rep(mt,each=16),rep(ctxt16,6),sep="_")
  types96 <- sapply(types96,function(z){sub("\\.",substr(z,1,1),z)})
  context <- substr(types96,4,6)

  nt <- sapply(types96,function(m){sum(vcf[which(vcf[,mutcat3.col]==m),strand.ts.col]=="nt",na.rm=T)})
  ts <- sapply(types96,function(m){sum(vcf[which(vcf[,mutcat3.col]==m),strand.ts.col]=="ts",na.rm=T)})
  tt <- rbind(ts,nt)
  pval <- as.numeric(apply(tt,2,function(z){try(prop.test(z[1],sum(z),p=0.5)$p.value,silent=TRUE)}))
  pval.lab <- rep("ns",length(pval));pval.lab[which(pval < 0.05)] <- "*";pval.lab[which(pval < 0.01)] <- "**";pval.lab[which(pval < 0.001)] <- "***"
  ymax <- max(tt)*1.2
  if(!is.null(plot.file)){pdf(plot.file,width=24,height=5)}
  col192 <- c(rep("skyblue3",32),rep("black",32),rep("red",32),rep("grey",32),rep("green",32),rep("pink",32))
  labs <-c(rep("C>A",32),rep("C>G",32),rep("C>T",32),rep("T>A",32),rep("T>C",32),rep("T>G",32))
  bp <- barplot(tt,beside= TRUE,col= c("royalblue","red"),border=c("royalblue","red"),legend.text= c("Transcribed","Untranscribed"),args.legend=list(x="topright",bty = "n",inset=c(0,0.1)),las=1,ylab="Number of mutations",xlab="Mutation type",ylim=c(0,ymax),xaxt="n",cex.lab=1.4)
  axis(1,at=apply(bp,2,mean),labels=context,pos=0,las=2,cex.axis=1.5,tick=T,cex.axis=1)
  for(i in seq(1,161,by=32)){
    rect(bp[i],par()$usr[4],bp[i+31],par()$usr[4]-0.05*diff(par()$usr[3:4]),col=col192[i],border=col192[i])
    text((bp[i]+bp[i+31])/2,par()$usr[4]+0.09*diff(par()$usr[3:4]),labels=labs[i], xpd= TRUE, cex = 2)
  }
  text(x=apply(bp,2,mean),y=apply(tt,2,max)+0.03*diff(par()$usr[3:4]),pval.lab,srt=90)
  if(!is.null(plot.file)){dev.off()}
  return(rbind(tt,pval))
}

#' rainfallz_plot
#'
#' Function to plot rainfall plots
#' @param vcf vcf data frame containing the mutations
#' @param mutype.col Substitution information column name in vcf
#' @param cyto required
#' @param resdir Result directory
#'
#' @return
#' @export
#'
#' @examples
rainfallz_plot <- function (vcf,mutype.col = "mutype",cyto = NULL, resdir = resdir)
{
  chrom.col = "CHROM"; pos.col = "POS"; sample.col = "Sample"; mutdist.col = "mutdist";colmut.col = "colmut"
  for (samp in unique(vcf[, sample.col])) {
    vcf. <- vcf[which(vcf[, sample.col] == samp), ]
    tmp <- cit.genomOrder(vcf., chrom = chrom.col, pos = pos.col)
    for (chr in unique(tmp[, chrom.col])) {
      ind <- which(tmp[, chrom.col] == chr)
      tmp[ind, mutdist.col] <- c(NA, diff(tmp[ind, pos.col]))
    }
    tmp$ycol <- log10(tmp[, mutdist.col])
    substype <- c("CA", "CG", "CT", "TA", "TC", "TG")
    mycol <- c("skyblue3", "black", "red", "grey", "green","pink")
    names(mycol) <- substype
    z <- as.factor(tmp[, mutype.col])
    levels(z) <- mycol[levels(z)]
    tmp[, colmut.col] <- as.character(z)
    pdf(file.path(resdir, paste(samp, "rainfall_plot.pdf",sep = "_")), width = 16, height = 6)
    cit.pangenomPlot(d = tmp, ycol = "ycol", chrom = chrom.col,
                      pos = pos.col, cyto = cyto, colorscol = colmut.col,
                      pch = 16, xlab = "Genomic Position", ylab = "-log10(mutation distance)")
    legend("top", horiz = T, legend = gsub("^(.{1})(.*)$","\\1>\\2", substype), pch = 19, col = mycol, bty = "n",inset = c(0, -0.16), xpd = T, cex = 1.2)
    dev.off()
  }
}

#' chrTime_plot
#'
#' Function to represent the timing of chromosomal duplications in point mutation time
#' @param vcf vcf data frame containing the mutations with multiplicity annotation, as obtained using chrTime_annot function
#' @param point.mut.time Timing of chromosome duplications in point mutation time, as obtained using chrTime_annot function
#' @param resdir Result directory
#'
#' @return
#' @export
#'
#' @examples
chrTime_plot <- function(vcf = NULL,
						 point.mut.time = NULL,
						 resdir=resdir)
  {
  sample.col <- "Sample"
  point.mut.time <- factoall(point.mut.time)
  for(samp in unique(point.mut.time[,sample.col]))
  {
    	resdir. <- file.path(resdir,samp);if(!file.exists(resdir.)){dir.create(resdir.)}
    	point.mut.time. <- point.mut.time[which(point.mut.time[,sample.col]==samp),]

    	pdf(file.path(resdir.,"Duplication_timing.pdf"),width=8,height=4)
    	par(mai=c(0.5,0.5,0.5,0.4))
    	plot(-10,xlim=c(0,105),ylim=c(-1,1),axes=F,ylab="",xlab="Point mutation time")
    	arrows(0,0,105,0)
    	segments(c(0,100),-0.15,c(0,100),0.15)
    	text(seq(0,100,by=20),-0.3,labels=seq(0,100,by=20),adj=c(0.5,1))
    	for(i in which(point.mut.time.[,sample.col]==samp)){
      		segments(point.mut.time.[i,"Moltime"],-0.15, point.mut.time.[i,"Moltime"],0.15,col=point.mut.time.[i,"Col"])
      		text(point.mut.time.[i,"Moltime"],0.5,labels=sub("gain","+",point.mut.time.[i,"CNA"]),adj=c(0.5,1),srt=45,col=point.mut.time.[i,"Col"])
   		}
    	text(50,-1,"Point mutation time (%)")
    	dev.off()

    	vcf. <- vcf[which(vcf$Sample==samp),]
    	VAF.color <- c("darkgrey","red","royalblue")[match(vcf.[,"Timing"],c(NA,"early","late"))]
    	png(file.path(resdir., "Somatic_mutations_logR_VAF.png"), width = 2400, height = 1600, res = 200)
    	palimpsest_cnvProfile(vcf = vcf.,VAF.color = VAF.color)
    	dev.off()
  }
}

#' palimpsest_cnvProfile
#'
#' Funciton to plot copy number variation plots
#' @param vcf vcf data frame containing the mutations
#' @param VAF.color vector of colors for the VAF plot (optionnal).
#'
#' @return
#' @export
#'
#' @examples
palimpsest_cnvProfile <- function(vcf=NULL,
								  VAF.color=NULL
)
{
  CHROM.col <- "CHROM";POS.col <- "POS";VAF.col <- "Tumor_Varprop";LogR.col = "LogR";ntot.col="ntot";nMin.col="Nmin"
  layout(matrix(1:3,3))
  d <- cit.pangenomPlot(vcf,LogR.col,CHROM.col,POS.col,cyto=cyto,pch=19,ylab="LogR",xlab="Genomic position",cex=0.2)
  if(is.null(VAF.color)){vcf$tmpcol <- rep("black",nrow(vcf))}else{vcf$tmpcol <- VAF.color}
  d <- cit.pangenomPlot(vcf,VAF.col,CHROM.col,POS.col,cyto=cyto,col="tmpcol",pch=19,ylab="VAF",xlab="Genomic position",cex=0.2)
  vcf$tmpcol <- "darkred"
  d <- cit.pangenomPlot(vcf,ntot.col,CHROM.col,POS.col,cyto=cyto,col="tmpcol",pch=".",cex=5,ylab="Absolute CN",xlab="Genomic position",ylim=c(0,min(10,max(vcf[,ntot.col],na.rm=T))))
  vcf$tmpcol <- "royalblue"
  d <- cit.pangenomPlot(vcf,nMin.col,CHROM.col,POS.col,cyto=cyto,col="tmpcol",pch=".",cex=5,ylab="Absolute CN",xlab="Genomic position",plotnew=F)
  legend("top",legend = c("Total Copy Number","Minor Allele Copy Number"),horiz = T,bty="n",fill = c("darkred","royalblue"),xpd=T,inset = c(0,-0.20))
}

#' palimpsest_Fireworks
#'
#' Function to plot fireworks plotfor variant allele fractions
#' @param vcf vcf data frame containing the mutations
#' @param CHROM.col Chromosome column name in vcf.
#' @param POS.col Position column name in vcf.
#' @param VAF.col Variant Allele Fraction information column name in vcf.
#' @param LogR.col LogR information column name in vcf.
#' @param nPurity.col Tumor purity column name in vcf.
#'
#' @return
#' @export
#'
#' @examples
palimpsest_Fireworks <- function(vcf=NULL,
                                  CHROM.col="CHROM",
                                  POS.col="POS",
                                  VAF.col="Tumor_Varprop",
                                  LogR.col="LogR",
                                  nPurity.col="nPurity")
{
  rho <- unique(vcf[,nPurity.col])
  colclust <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#bd18ea", "#2ef4ca", "#f4cced", "#f4cc03", "#05188a", "#e5a25a", "#06f106", "#85848f", "#000000", "#076f25", "#93cd7f", "#4d0776", "darkred", "darkgrey")
  vcf$chromcol <- colclust[match(vcf[,CHROM.col],paste("chr",c(1:22,"X","Y","M"),sep=""))]
  zones=matrix(1:2,2)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  par(mar=c(0,4.5,1,4.5))
  hist(vcf[,VAF.col],xaxt='n',xlab = "",main="",breaks =100 ,col="royalblue",border="royalblue",las=2,xlim=c(0,1))
  par(mar=c(6.4,4.5,1,4.5))
  plot(vcf[,VAF.col],vcf[,LogR.col],xlim=c(0,1),pch=".",cex=2,col=vcf$chromcol,ylim=c(-2,3), xlab="VAF",ylab = "logR")#,ylim=c(max(c(quantile(vcf$LogR,probs=0.01),-0.5)),min(c(quantile(vcf$LogR,probs=0.99),2))),las=2)
  abline(v=rho/2,lty=2)
  legend("right",legend = paste("chr",c(1:22,"X","Y","M"),sep=""),bty="n",xpd=T,pch=16,col = colclust,inset = c(-0.12,1),cex = 0.8)

}

#' palimpsest_CCFplot
#'
#' Function to plot Variant Allele Fraction and CCF distribution in a Sample
#' @param vcf vcf data frame containing the mutations
#' @param CHROM.col Chromosome column name in vcf.
#' @param VAF.col Variant Allele Fraction information column name in vcf.
#' @param clonality.col Clonaility information column name in vcf.
#'
#' @return
#' @export
#'
#' @examples
palimpsest_CCFplot <- function(vcf=NULL,
                            CHROM.col="CHROM",
                            VAF.col="Tumor_Varprop",
                            clonality.col="CCF_Clonality")
  {
  par(mfrow=c(2,1), mai = c(1, 1, 0.4, 1))
  hist.var <- hist(vcf[,VAF.col],br=seq(0,1,by=0.02),plot=F)
  ymax.var <- max(hist.var$counts)
  hist(vcf[,VAF.col][which(vcf[,clonality.col]=="clonal")],br=seq(0,1,by=0.02),xlim=c(0,1),col="tomato",ylim=c(0,ymax.var), xlab = "Tumor VAF",main="Distribution of clonal and subclonal mutations")
  hist(vcf[,VAF.col][which(vcf[,clonality.col]=="subclonal")],br=seq(0,1,by=0.02),xlim=c(0,1),add=T,col="steelblue")
  legend("topright",legend = c("Clonal","Subclonal"),bty="n",fill = c("tomato","steelblue"),border = c("tomato","steelblue"),xpd=T,inset = c(-0.05,-0.14))
  hist.ccf <-  hist(vcf$CCF,br=seq(0,2,by=0.02), plot = F)
  ymax.ccf <- max(hist.ccf$counts)
  hist(vcf$CCF[which(vcf[,clonality.col]=="clonal")],br=seq(0,2,by=0.02),xlim=c(0,max(vcf$CCF,na.rm = T)),ylim=c(0,ymax.ccf),col="tomato", main="", xlab = "Cancer Cell Fraction (CCF)")
  hist(vcf$CCF[which(vcf[,clonality.col]=="subclonal")],br=seq(0,2,by=0.02),xlim=c(0,max(vcf$CCF,na.rm = T)),add=T,col="steelblue")

}

#' palimpsest_DissectSigs
#'
#' Function to plot temporal dissection of mutational signatures within a sample
#' @param vcf vcf data frame containing the mutations
#' @param signatures_exp_clonal Mutational signature exposures extracted from clonal mutations only using deconvolution_fit function.
#' @param signatures_exp_subclonal Mutational signature exposures extracted from subclonal mutations only using deconvolution_fit function.
#' @param sig_cols Character vector indicating the colors representing each signature in graphical outputs. Must match to the total number of provided signatures
#' @param resdir Result directory
#'
#' @return
#' @export
#'
#' @examples
palimpsest_DissectSigs <- function(vcf=NULL,
                                   signatures_exp_clonal=NULL,
                                   signatures_exp_subclonal=NULL,
                                   sig_cols=NULL,
                                   resdir=NULL)
{
  substype <- c("CA","CG","CT","TA","TC","TG");substype_cols <- c("skyblue3", "black", "red", "grey", "green", "pink")
  for (samp in unique(vcf[,"Sample"])) {
    vcf. <- vcf[which(vcf$Sample==samp & vcf$Type=="SNV"),]
    resdir. <- file.path(resdir,samp);if(!file.exists(resdir.)){dir.create(resdir.)}
    mat <- matrix(0,6,2);rownames(mat) <- substype;colnames(mat)<-c("clonal","subclonal");
    tt <- table(vcf.[,"substype"],vcf.[,"Clonality"])
    mat[rownames(tt),colnames(tt)] <- tt
    ptt <- prop.table(mat,2)
    nearly <- sum(vcf.[,"Clonality"]=="clonal",na.rm=T);nlate <- sum(vcf.[,"Clonality"]=="subclonal",na.rm=T)
    pdf(file.path(resdir.,"Clonal_vs_Subclonal_6_substitution_types.pdf"),width=4,height=4)
    barplot(ptt,col=substype_cols,las=1,ylab="Proportion of mutations",names.arg=c(paste("clonal\nn=",nearly,sep=""),paste("subclonal\nn=",nlate,sep="")),xlim=c(0,4))
    title(paste("p=",chisq.test(tt)$p.value))
    legend("topright",legend=c("C>A","C>G","C>T","T>A","T>C","T>G"),fill=substype_cols,bty="n")
    dev.off()
    pdf(file.path(resdir.,"Clonal_vs_Subclonal_96_substitution_types.pdf"),width=24,height=10)
    layout(matrix(1:2,2))
    plot96mutationSpectrumFromVcf(vcf.[which(vcf.[,"Clonality"]=="clonal"),],sample.col="Sample")
    plot96mutationSpectrumFromVcf(vcf.[which(vcf.[,"Clonality"]=="subclonal"),],sample.col="Sample")
    dev.off()
    pdf(file.path(resdir.,"Clonal_vs_Subclonal_signatures.pdf"),width=5,height=4)
    mat <- t(as.matrix(rbind(signatures_exp_clonal$sig_nums[samp,],signatures_exp_subclonal$sig_nums[samp,])))
    ind <- which(apply(mat,1,sum) > 0)
    mat <- matrix(mat[ind,],ncol=ncol(mat),byrow=F,dimnames=list(rownames(mat)[ind],colnames(mat)))
    ptt <- prop.table(mat,2)
    barplot(ptt,col=sig_cols[rownames(mat)],las=1,ylab="Proportion of mutations",names.arg=c(paste("clonal\nn=",nearly,sep=""),paste("subclonal\nn=",nlate,sep="")),xlim=c(0,4.5))
    title(paste("p=",chisq.test(mat)$p.value))
    legend("topright",legend=rownames(mat),fill=sig_cols[rownames(mat)],bty="n")
    dev.off()
  }
}

#' cnaCCF_plots
#'
#' @param vcf vcf data frame containing the mutations
#' @param resdir Result directory
#'
#' @return
#' @export
#' @examples
cnaCCF_plots <- function (vcf = NULL, resdir = resdir)
{
  sample.col = "Sample"; CHROM.col = "CHROM";
  POS.col = "POS"; VAF.col = "Tumor_Varprop"; LogR.col = "LogR";
  clonality.col = "Clonality"; ntot.col = "ntot"; nPurity.col = "Purity";
  for (samp in unique(vcf[, sample.col])) {
    vcf. <- vcf[which(vcf[, sample.col] == samp), ]
    resdir. <- file.path(resdir, samp)
    if (!file.exists(resdir.)) {
      dir.create(resdir.)
    }
    pdf(file.path(resdir., "Fireworks_Plot.pdf"), width = 8,
        height = 8)
    palimpsest_Fireworks(vcf = vcf., CHROM.col = CHROM.col,
                         POS.col = POS.col, VAF.col = VAF.col, LogR.col = LogR.col,
                         nPurity.col = nPurity.col)
    dev.off()
    pdf(file.path(resdir., "VAF_CCF_plot.pdf"), width = 8,
        height = 6)
    palimpsest_CCFplot(vcf = vcf., CHROM.col = CHROM.col,
                       VAF.col = VAF.col, clonality.col = clonality.col)
    dev.off()
  }
}

#' deconvolution_exposure
#'
#' Function to plot the exposure of mutational signatures in samples across the entire series
#' @param mutSign_nums Matrix in sample x mutational signature exposure format in numbers
#' @param mutSign_props Matrix in sample x mutational signature exposure format in proportions
#' @param sig_cols Character vector indicating the colors representing each signature in graphical outputs. Must match to the total number of provided signatures
#'
#' @return
#' @export
#' @import ggplot2
#' @importFrom ggplot2 theme_bw
#' @import gplots
#' @import reshape2
#' @importFrom reshape2 melt
#' @examples
deconvolution_exposure <- function(mutSign_nums,mutSign_props,sig_cols=colclust) {
  scale<-1
  .theme_ss <- theme_bw(base_size=14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=8*scale, family="mono"),
          axis.text.y = element_text(hjust = 0.5,size=12*scale, family="mono"),
          axis.text = element_text(size = 12*scale, family = "mono"))
  ordering <- order(colSums(t(mutSign_nums)),decreasing=T)
  mutSign_nums <- t(mutSign_nums)
  mutSign_props <- t(mutSign_props)
  mutSign_nums <- mutSign_nums[,ordering]
  mutSign_props <- mutSign_props[,ordering]
  sample.ordering <- colnames(mutSign_nums)
  x1 <- melt(mutSign_nums)
  x2 <- melt(mutSign_props)
  colnames(x1) <- c("Signature","Sample","Activity")
  colnames(x2) <- c("Signature","Sample","Activity")
  x1[,"class0"] <- c("Counts")
  x2[,"class0"] <- c("Proportions")
  df2 <- rbind(x1,x2)
  df2$class0 <- factor(df2$class0,c("Counts","Proportions"))
  df2$Sample <- factor(df2$Sample,sample.ordering)
  p = ggplot(df2,aes(x=factor(Sample),y=Activity,fill=Signature))
  p = p+geom_bar(stat="identity",position='stack')
  p = p + scale_fill_manual(values=sig_cols)
  p = p + ggtitle("Mutational Signature Exposures")
  p = p + facet_grid(class0 ~ ., scale = "free_y")
  p = p + theme(plot.title=element_text(lineheight=1.0,face="bold",size=15*scale))
  p = p + xlab("Samples") + ylab("Mutational Signature Content")
  p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=15*scale))
  p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=15*scale))
  p = p + theme(axis.text.x = element_text(angle=90,vjust=0.5,size=12*scale,face="bold",colour="black"))
  p = p + theme(axis.text.y = element_text(size=10*scale,face="bold",colour="black"))
  p = p + theme(legend.title=element_blank())
  p = p + .theme_ss
  p = p + theme(legend.position="top")
  return(p)
}


#' plot.SNV.sigs
#'
#' Function to plot single nucleotide variants signatures (96 nucleotide type)
#' @param spec Matrix of mutational signatures x  mutation types
#'
#' @return
#' @export
#'
#' @examples
plot.SNV.sigs <- function(spec){
  bases <- c("A", "C", "G", "T")
  ctxt16 <- paste(rep(bases, each = 4), rep(bases, 4), sep = ".")
  mt <- c("CA", "CG", "CT", "TA", "TC", "TG")
  types96 <- paste(rep(mt, each = 16), rep(ctxt16, 6), sep = "_")
  types96 <- sapply(types96, function(z) {
    sub("\\.", substr(z, 1, 1), z)
  })
  context <- substr(types96, 4, 6)
  col96 <- c(rep("skyblue3",16),rep("black",16),rep("red",16),rep("grey",16),rep("green",16),rep("pink",16))
  labs <-c(rep("C>A",16),rep("C>G",16),rep("C>T",16),rep("T>A",16),rep("T>C",16),rep("T>G",16))
  colclust <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928","gray36", "#bd18ea", "#2ef4ca", "#f4cced", "#f4cc03", "#05188a","mediumvioletred", "#e5a25a", "#06f106", "#85848f", "#000000", "blue","#076f25", "#93cd7f", "#4d0776", "darkred", "yellow","cadetblue1")
  for(i in 1:nrow(spec)){
    ymax <- ceiling(max(spec[i,])*100)/100
    bp <- barplot(as.numeric(spec[i,]),col=col96,border=col96,las=2,width=1,space=1,yaxt="n",xaxt="n",ylim=c(0,ymax*1.2),xlab = rownames(spec)[i],cex.lab=1.4,ylab="Percentage of mutations")
    axis(1,at=bp,labels=context,pos=0,las=2,cex.axis=1.5,tick=T,cex.axis=1)
    axis(2,at=round(seq(0,ceiling(ymax*100),length.out=3),digits=1)/100,pos=0,las=1,cex.axis=1.5)
    for(i in seq(1,81,by=16)){
      rect(bp[i],par()$usr[4],bp[i+15],par()$usr[4]-0.05*diff(par()$usr[3:4]),col=col96[i],border=col96[i])
      text((bp[i]+bp[i+15])/2,par()$usr[4]+0.09*diff(par()$usr[3:4]),labels=labs[i], xpd= TRUE, cex = 2)
    }
  }
}

#' plot.SV.sigs
#'
#' Function to plot structural variants signatures (38 sv categories type)
#' @param sigs Matrix of mutational signatures x  mutation types
#'
#' @return
#' @export
#'
#' @examples
plot.SV.sigs <- function(sigs){
  sigs <- t(sigs)
  col15 <- rep(c(rep("red",6),rep("olivedrab1", 6), rep("lightskyblue", 6),"grey"),2)
  myord <- c("DEL_0-1kb_clust_1","DEL_1-10kb_clust_1","DEL_10-100kb_clust_1","DEL_100kb-1Mb_clust_1","DEL_1-10Mb_clust_1","DEL_>10Mb_clust_1","DUP_0-1kb_clust_1","DUP_1-10kb_clust_1","DUP_10-100kb_clust_1","DUP_100kb-1Mb_clust_1","DUP_1-10Mb_clust_1","DUP_>10Mb_clust_1","INV_0-1kb_clust_1","INV_1-10kb_clust_1","INV_10-100kb_clust_1","INV_100kb-1Mb_clust_1","INV_1-10Mb_clust_1","INV_>10Mb_clust_1","BND_clust_1","DEL_0-1kb_clust_0","DEL_1-10kb_clust_0","DEL_10-100kb_clust_0","DEL_100kb-1Mb_clust_0","DEL_1-10Mb_clust_0","DEL_>10Mb_clust_0","DUP_0-1kb_clust_0","DUP_1-10kb_clust_0","DUP_10-100kb_clust_0","DUP_100kb-1Mb_clust_0","DUP_1-10Mb_clust_0","DUP_>10Mb_clust_0","INV_0-1kb_clust_0","INV_1-10kb_clust_0","INV_10-100kb_clust_0","INV_100kb-1Mb_clust_0","INV_1-10Mb_clust_0","INV_>10Mb_clust_0","BND_clust_0")
  labs=rep(c(rep("del",6),rep("tds",6),rep("inv",6),"trans"),2)
  sigs <- matrix(sigs[myord,],ncol=ncol(sigs),dimnames=dimnames(sigs),byrow=F)
  if(max(range(sigs)) > 1) {axis_name <- "Number of mutations"}
  if(max(range(sigs)) <= 1) {axis_name <- "Percentage of mutations"}
  for(siggy in colnames(sigs))
  {
    par(mai=c(1.2,1,2,1))
    mylim <- max(sigs[,siggy])*1.2
    bp <- barplot(sigs[,siggy], col = col15, border = col15, las = 2, space = c(rep(1,19),5,rep(1,18)),ylim=c(0,mylim), xaxt = "n", ylab = axis_name, main = siggy)
    abline(v=mean(bp[19:20]))
    axis(1, at = bp, labels = sub("clust","",sapply(myord,function(z)	unlist(strsplit(z,"_"))[2])), pos = 0, las = 2, cex.axis = 1.5, tick = T, cex.axis = 1)
    for (i in seq(1, 13, by = 6)) {
      rect(bp[i], par()$usr[4], bp[i]+10, par()$usr[4] - 0.05 * diff(par()$usr[3:4]), col = col15[i], border = col15[i])
      text((bp[i] + bp[i]+8)/2, par()$usr[4] + 0.09 * diff(par()$usr[3:4]), labels = labs[i], xpd = TRUE, cex = 1)
    }
    for (i in seq(20, 32, by = 6)) {
      rect(bp[i], par()$usr[4], bp[i]+10, par()$usr[4] - 0.05 * diff(par()$usr[3:4]), col = col15[i], border = col15[i])
      text((bp[i] + bp[i]+8)/2, par()$usr[4] + 0.09 * diff(par()$usr[3:4]), labels = labs[i-16], xpd = TRUE, cex = 1)
    }
    for (i in c(19,38)) {
      rect(bp[i]-1, par()$usr[4], bp[i]+1, par()$usr[4] - 0.05 * diff(par()$usr[3:4]), col = "grey", border = "grey")
      text((bp[i]-1 + bp[i]+1)/2, par()$usr[4] + 0.09 * diff(par()$usr[3:4]), labels = "trans", xpd = TRUE, cex = 1)
    }
  }
}

#' palimpsest_plotTumorHistories
#'
#' Function to represent the natural histories of tumors
#' @param vcf Annotated vcf file for mutations
#' @param sv.vcf Annotated vcf file for SVs
#' @param cna_data Data frame with Copy Number Alteration Information
#' @param point.mut.time Duplication timing estimated from point mutations
#' @param clonsig Matrix giving the proportion of clonal mutations attributed to each signature
#' @param subsig Matrix giving the proportion of subclonal mutations attributed to each signature
#' @param msigcol Color codes for mutational signatures
#' @param msigcol.sv Color codes for SV signatures
#' @param spacechar Space for characters in the legends
#' @param ystep Step between annotation lines
#' @param resdir
#' @return
#' @export
#' @importFrom Rgraphviz pieGlyph
#' @importFrom plotrix draw.circle
#' @import plotrix
#' @examples
palimpsest_plotTumorHistories <- function(vcf=NULL,
							 sv.vcf=NULL,
							 cna_data=NULL,
							 point.mut.time=NULL,
							 clonsig=NULL,
							 subsig=NULL,
							 msigcol=NULL,
							 msigcol.sv=NULL,
							 spacechar=1,
							 ystep=2,
							 resdir=NULL
							 )
{
	nMin.col = "Nmin";nMaj.col = "Nmaj"

	for(Sample_to_plot in unique(vcf$Sample)){
		xposstart <- 0

	  	vcf. <- vcf[which(vcf$Sample==Sample_to_plot),]

		pdf(file.path(resdir,paste0(Sample_to_plot,".pdf")),width=12,height=5)

		# Plot frame
		posClon <- 65;posSub <- 100
		par(mai=rep(0,4))
		plot(-100,xlim=c(0,110),ylim=c(-22,27),axes=F,bty="n",xlab="",ylab="")
		rect(0,-2,posClon,2,col="aliceblue",border="aliceblue")
		rect(posClon,-1,posSub,1,col="steelblue1",border="steelblue1")
		draw.circle(posClon,0,3,border="royalblue",col="royalblue")
		draw.circle(posSub,0,3,border="darkblue",col="darkblue")
    	text(35,0,paste(sum(vcf.$Clonality=="clonal",na.rm=T),"mutations"))
	    text(82,0,paste(sum(vcf.$Clonality=="subclonal",na.rm=T),"mutations"))

		# Write annot summary
	    legend("topleft",legend = Sample_to_plot,bty="n")

		# Write mutations
		yposclon <- ypossub <- 6
		#clonal
		xpos <- xposstart
		ind <- which(vcf$Sample== Sample_to_plot & !is.na(vcf$Driver) & vcf$Clonality=="clonal")
		ind <- ind[match(unique(vcf[ind,"Driver"]),vcf[ind,"Driver"])]
		if(length(ind)){
			for(i in ind){
				text(xpos,yposclon,vcf[i,"Driver"],pos=4,col=msigcol[vcf[i,"Sig.max"]])  #vcf[i,"MutSig"]
				xpos <- xpos+strwidth(vcf[i,"Driver"])*spacechar
				if(i!=max(ind)){
					text(xpos,yposclon,",",pos=4,col="black")
					xpos <- xpos+1*spacechar
					if(xpos > 50){xpos <- xposstart;yposclon <- yposclon-ystep}
				}
			}
		}
		#subclonal
		xpos <- 66
		ind <- which(vcf$Sample== Sample_to_plot & !is.na(vcf$Driver) & vcf$Clonality=="subclonal")
		ind <- ind[match(unique(vcf[ind,"Driver"]),vcf[ind,"Driver"])]
		if(length(ind)){
			for(i in ind){
				text(xpos,ypossub,vcf[i,"Driver"],pos=4,col=msigcol[vcf[i,"Sig.max"]])  #vcf[i,"MutSig"]
				xpos <- xpos+strwidth(vcf[i,"Driver"])*spacechar
				if(i!=max(ind)){
					text(xpos,yposclon,",",pos=4,col="black")
					xpos <- xpos+1*spacechar
					if(xpos > 90){xpos <- 66;ypossub <- ypossub-ystep}
				}
			}
		}
		yposclon <- ypossub <- min(c(yposclon,ypossub))

		# Pie charts mutational signatures
		sig <- clonsig[Sample_to_plot,]
		sig <- unlist(sig[which(sig> 0)])
		Rgraphviz::pieGlyph(sig,xpos=50,ypos=12,col=msigcol[names(sig)],border=msigcol[names(sig)],radius=4,labels=sub("Signature.","",names(sig)),cex=0.8)
		text(68,18,"Mutational signatures",pos=1)

		sig <- subsig[Sample_to_plot,]
		sig <- unlist(sig[which(sig> 0)])
		Rgraphviz::pieGlyph(sig,xpos=85,ypos=12,col=msigcol[names(sig)],border=msigcol[names(sig)],radius=4,labels=sub("Signature.","",names(sig)),cex=0.8)

		# Write CNAs
		#duplications (with timing)
		for(i in which(point.mut.time$Sample== Sample_to_plot)){
			segments(point.mut.time[i,"Moltime"]*(posClon-3)/100,-2, point.mut.time[i,"Moltime"]*(posClon-3)/100,-3,col="indianred1")#
			text(point.mut.time[i,"Moltime"]*(posClon-3)/100,-3,labels=sub("gain ","+",point.mut.time[i,"CNA"]),adj=c(1.2,1.2),srt=45,col="indianred1")#adj=c(0.5,1),point.mut.time[i,"Col"]
		}
		#deletions
		xpos <- xposstart;yposclon <- -14
		cna_data$length <- cna_data[,"POS_END"]-cna_data[,"POS_START"]
		ind <- which(cna_data$Sample== Sample_to_plot & cna_data[, nMaj.col] <= 1 & cna_data[, nMin.col] == 0 & cna_data$length >= 1e6)
		dels <- unique(as.character(as.matrix(cna_data[ind,c("startarm","endarm")])))
    	mychrs <- sub("p|q", "",dels)
	    pq <- names(which(table(mychrs) > 1))
		if(length(pq)){
			dels <- c(setdiff(dels,c(paste0(pq,"p"),paste0(pq,"q"))),pq)
			dels <- dels[order(as.numeric(gsub("p|q","",(sub("X",23,sub("Y",24, dels))))))]
		}
		if(length(dels)){
			for(d in dels){
				mut <- paste0("-",d)
				text(xpos,yposclon,mut,pos=4,col="royalblue")#colclust[match(gsub("p|q","",d),c(1:22,"X","Y"))]
				xpos <- xpos+strwidth(mut)*spacechar
				if(d!=tail(dels,1)){
					text(xpos,yposclon,",",pos=4,col="black")
					xpos <- xpos+1*spacechar
					if(xpos > 55){xpos <- xposstart;yposclon <- yposclon-ystep}
				}
			}
		}

		# Write SVs
		ypossub <- yposclon <- yposclon-ystep
		# clonal
		xpos <- xposstart
		ind <- which(sv.vcf$Sample== Sample_to_plot & !is.na(sv.vcf$Driver))# & sv.vcf$Clonality=="clonal"
		ind <- ind[match(unique(sv.vcf[ind,"Driver"]),sv.vcf[ind,"Driver"])]
		if(length(ind)){
			mut <- "SVs:"
			text(xpos,yposclon,mut,pos=4)  #vcf[i,"MutSig"]
			xpos <- xpos+(strwidth(mut)+1)*spacechar
			for(i in ind){
				text(xpos,yposclon,sv.vcf[i,"Driver"],pos=4,col=msigcol.sv[sv.vcf[i,"Sig.max"]])  #vcf[i,"MutSig"]
				xpos <- xpos+strwidth(sv.vcf[i,"Driver"])*spacechar
				if(i!=max(ind)){
					text(xpos,yposclon,",",pos=4,col="black")
					xpos <- xpos+1*spacechar
					if(xpos > 50){xpos <- xposstart;yposclon <- yposclon-ystep}
				}
			}
		}
		#subclonal
		#xpos <- 66
		#ind <- which(sv.vcf$Sample== Sample_to_plot & !is.na(sv.vcf$Driver) & sv.vcf$Clonality=="subclonal")
		#ind <- ind[match(unique(sv.vcf[ind,"Driver"]),sv.vcf[ind,"Driver"])]
		#if(length(ind)){
		#	mut <- "SVs:"
		#	text(xpos,ypossub,mut,pos=4)  #vcf[i,"MutSig"]
		#	xpos <- xpos+(strwidth(mut)+1)*spacechar
		#	for(i in ind){
		#		text(xpos,ypossub,sv.vcf[i,"Driver"],pos=4,col=msigcol.sv[sv.vcf[i,"Sig.max"]])  #vcf[i,"MutSig"]
		#		xpos <- xpos+strwidth(sv.vcf[i,"Driver"])*spacechar
		#		if(i!=max(ind)){
		#			text(xpos,ypossub,",",pos=4,col="black")
		#			xpos <- xpos+1*spacechar
		#			if(xpos > 50){xpos <- xposstart;ypossub <- ypossub-ystep}
		#		}
		#	}
		#}
	dev.off()
	}
}



#' palimpsest_clonalitySigsCompare
#'
#' Function to plot comparison between mutational signatures exposures within clonal vs. subclonal mutations in the series.
#' @param clonsig Data frame of sample x clonal mutational signature exposure format in numbers
#' @param subsig Data frame of sample x subclonal mutational signature exposure format in numbers
#' @param msigcol list of colors to for plotting signature contribution. Must match to the total number of provided signatures
#' @param resdir Result directory
#'
#' @return
#' @export
#' @examples
palimpsest_clonalitySigsCompare <- function(clonsig=clonsig,subsig=subsig,msigcol=msigcol,resdir=resdir){
  clonsigp <- clonsig/apply(clonsig,1,sum)
  subsigp <- subsig/apply(subsig,1,sum)
  labs <- c("clonal","subclonal")
  pdf(file.path(resdir,"Clonal_vs_sublclonal_signature_proportions.pdf"),width=10,height=3.5)
  bp <- plot(-10,xlim=c(0,ncol(clonsigp)*2+4),ylim=c(-0.1,1.2),axes=F,bty="n",xlab="",ylab="Proportion of mutations",main="Clonal vs.Subclonal Signature Exposures")
  for(j in 1:ncol(clonsigp)){
    sig <- colnames(clonsigp)[j]
    ind <- which(clonsigp[,sig] > 0 | subsigp[,sig] > 0)
    tmp <- data.frame(clon=clonsigp[ind,sig],sub <- subsigp[ind,sig])
    for(i in 1:nrow(tmp)){
      points(((j-1)*2.5):((j-1)*2.5+1),tmp[i,],type="b",pch=c(19),col=msigcol[sig])
      axis(1,at=((j-1)*2.5):((j-1)*2.5+1),labels=labs,pos=0,las=2,cex.axis=1.5,tick=FALSE,cex.axis=1)
    }
    text(mean(((j-1)*2.5):((j-1)*2.5+1)),1.2,sub("Signature.","",sig),col=msigcol[sig])
  }
  axis(2,at=c(0,0.5,1),las=2)
  abline(h=0)
  dev.off()
  avclon <- apply(clonsigp,2,mean);avclon
  avsub <- apply(subsigp,2,mean);avsub
  pdf(file.path(resdir,"Clonal_vs_sublclonal_signature_series.pdf"),width=12,height=5)
  par(mfrow=c(1,2),xpd=NA)
  pie(avclon,col=msigcol[names(avclon)],border = msigcol[names(avclon)],labels=sub("Signature.","",names(avclon)),main = "CLONAL")
  pie(avsub,col=msigcol[names(avsub)],border = msigcol[names(avsub)],labels=sub("Signature.","",names(avsub)),main = "SUBCLONAL")
  legend(-3.8,-1.35,legend = names(mycol),ncol = 5,fill=mycol,border = mycol,bty = "n",cex=1)
  dev.off()
}


#' RCircos.Heatmap.Plot_1
#'
#' @param heatmap.data
#' @param data.col
#' @param track.num
#' @param side
#' @param min.value
#' @param max.value
#' @param inside.pos
#' @param outside.pos
#' @param genomic.columns
#' @param is.sorted
#'
#' @return
#' @export
#' @import RCircos
#' @examples
RCircos.Heatmap.Plot_1 <- function(heatmap.data=NULL, data.col=NULL,
                                   track.num=NULL, side=c("in", "out"), min.value=NULL, max.value=NULL,
                                   inside.pos=NULL, outside.pos=NULL, genomic.columns=3, is.sorted=TRUE)
{
  tumor = "Sample";chr1 = "CHROM_1";chr2 = "CHROM_2";pos1 = "POS_1";pos2 = "POS_2";event = "Type";
  requireNamespace("RCircos", quietly = TRUE)
  if(is.null(heatmap.data))
    stop("Genomic data missing in RCircos.Heatmap.Plot().\n");
  if(is.null(genomic.columns))
    stop("Missing number of columns for genomic position.\n");
  if( is.null(data.col) || data.col <= genomic.columns)
    stop("Data column must be ", genomic.columns+1, " or bigger.\n");
  boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos,outside.pos, FALSE);
  outerPos <- boundary[1];
  innerPos <- boundary[2];
  RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
  RCircos.Pos <- RCircos.Get.Plot.Positions();
  RCircos.Par <- RCircos.Get.Plot.Parameters();
  colorMap <- RCircos.Get.Heatmap.Color.Scale(RCircos.Par$heatmap.color);
  if(is.null(min.value) || is.null(max.value))
  {
    columns <- (genomic.columns+2):ncol(heatmap.data);
    min.value <- min(as.matrix(heatmap.data[, columns]));
    max.value <- max(as.matrix(heatmap.data[, columns]));
  }
  colorLevel  <- seq(min.value, max.value, length=length(colorMap));
  heatmap.data <- RCircos.Get.Single.Point.Positions(heatmap.data,genomic.columns);
  plotLocations <- RCircos.Get.Start.End.Locations(heatmap.data,RCircos.Par$heatmap.width);
  chromosomes <- unique(as.character(RCircos.Cyto$Chromosome));
  outlineColors <- rep("white", length(chromosomes));
  RCircos.Track.Outline(outerPos, innerPos, num.layers=1,
                        chrom.list=chromosomes, track.colors=outlineColors);
  heatmapValues <- as.numeric(heatmap.data[, data.col]);
  for(aPoint in 1:length(heatmapValues))
  {
    theLevel <- which(colorLevel >= heatmapValues[aPoint]);
    cellColor <- colorMap[min(theLevel)];
    if(heatmapValues[aPoint] == 1){
      cellColor = "olivedrab3"
    }else if (heatmapValues[aPoint] == 2){
      cellColor = "red"
    }else{
      cellColor = "royalblue1"
    }
    theStart <- plotLocations[aPoint, 1];
    theEnd   <- plotLocations[aPoint, 2];
    polygonX <- c(RCircos.Pos[theStart:theEnd,1]*outerPos,
                  RCircos.Pos[theEnd:theStart,1]*innerPos);
    polygonY <- c(RCircos.Pos[theStart:theEnd,2]*outerPos,
                  RCircos.Pos[theEnd:theStart,2]*innerPos);
    polygon(polygonX, polygonY, col=cellColor, border=NA);
  }
}


#' RCircos.Link.Plot_1
#'
#' @param link.data
#' @param track.num
#' @param by.chromosome
#' @param start.pos
#' @param genomic.columns
#' @param is.sorted
#' @param lineWidth
#'
#' @return
#' @export
#' @import RCircos
#' @examples
RCircos.Link.Plot_1 <- function(link.data=NULL, track.num=NULL,
                                by.chromosome=FALSE, start.pos=NULL,
                                genomic.columns=3, is.sorted=TRUE,
                                lineWidth=rep(1, nrow(link.data)))
{
  requireNamespace("RCircos", quietly = TRUE)
  tumor = "Sample";chr1 = "CHROM_1";chr2 = "CHROM_2";pos1 = "POS_1";pos2 = "POS_2";event = "Type";
  if(is.null(link.data)) stop("Link data missing in RCircos.Link.Plot().\n");
  if(by.chromosome!=TRUE && by.chromosome!=FALSE)
    stop("Error: by.chromosome must be either TRUE or FALSE.\n");
  if(length(which(lineWidth < 1)) > 1)
    stop("Line width must be positive integer.")
  RCircos.Par <- RCircos.Get.Plot.Parameters();
  if(is.null(start.pos)) {
    locations <- RCircos.Get.Plot.Boundary(track.num, side="in",inside.pos=NULL, outside.pos=NULL, FALSE);
    line.start <- locations[1];
  } else {
    if(start.pos>=1) stop("Link line must be inside chromosome ideogram");
    line.start <- RCircos.Par$chr.ideo.pos * start.pos;
  }
  if(is.null(genomic.columns) || genomic.columns < 3)
    stop("Incorrect number of columns for genomic position.\n");
  link.data <- RCircos.Get.Paired.Points.Positions(link.data,genomic.columns, plot.type="link");
  link.colors <- RCircos.Get.Link.Colors(link.data, genomic.columns,by.chromosome);
  RCircos.Pos <- RCircos.Get.Plot.Positions();
  base.positions <- RCircos.Pos[, 1:2]*line.start;
  for(a.link in seq_len(nrow(link.data)))
  {
    point.one <- as.numeric(link.data$LinkStart[a.link]);
    point.two <- as.numeric(link.data$LinkEnd[a.link]);
    if(point.one > point.two)
    {
      point.one <- link.data$LinkEnd[a.link];
      point.two <- link.data$LinkStart[a.link];
    }
    P0 <- as.numeric(base.positions[point.one, ]);
    P2 <- as.numeric(base.positions[point.two, ]);
    links <- RCircos.Link.Line(P0, P2);
    lines(links$pos.x, links$pos.y, type="l",
          lwd=lineWidth[a.link], col="gray70" );
  }
}


#' generate_table
#'
#' @param table
#'
#' @return
#' @export
#'
#' @examples
generate_table <- function(table){
  tumor = "Sample";chr1 = "CHROM_1";chr2 = "CHROM_2";pos1 = "POS_1";pos2 = "POS_2";event = "Type";
  translocs = table[which(table$Type == "BND"),]
  translocs = translocs[,c(chr1, pos1, pos1, chr2, pos2, pos2, tumor)]
  translocs[,2] = as.numeric(translocs[,2])
  translocs[,3] = as.numeric(translocs[,3])
  translocs[,5] = as.numeric(translocs[,5])
  translocs[,6] = as.numeric(translocs[,6])
  others = table[which(table$Type != "BND"),]
  others = others[, c(chr1, pos1, pos2, event, tumor)]
  others[,pos1] = as.numeric(others[,pos1])
  others[,pos2] = as.numeric(others[,pos2])
  return(list(translocs, others))
}

#' make_plot
#'
#' @param translocs
#' @param others
#' @param Sample_to_plot
#' @param resdir
#'
#' @return
#' @export
#' @import RCircos
#' @examples
make_plot <- function(translocs, others, Sample_to_plot,resdir){
  tumor = "Sample";chr1 = "CHROM_1";chr2 = "CHROM_2";pos1 = "POS_1";pos2 = "POS_2";event = "Type";
  requireNamespace("RCircos", quietly = TRUE)
  transl = translocs[which(translocs[,tumor] == Sample_to_plot),c(1:6)]
  othe = others[which(others[,tumor] == Sample_to_plot),c(1:4)]
  othe = cbind(othe,c(rep(NA, nrow(othe))))
  othe[,5][which(othe[,4] == "DUP")] = 1
  othe[,5][which(othe[,4] == "DEL")] = 2
  othe[,5][which(othe[,4] == "INV")] = 3
  chr.exclude <- NULL
  cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
  tracks.inside <- 4
  tracks.outside <- 0
  RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)
  png(file.path(resdir,paste0("Circos_",Sample_to_plot , ".png")), width = 2880, height = 2880, res = 400)
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
  track.num <- 1
  side <- "in"
  if(nrow(othe) >0) RCircos.Heatmap.Plot_1(othe, data.col = 5,track.num, side, inside.pos=0.75, outside.pos=1)
  track.num <- 4
  if(nrow(transl) >0) RCircos.Link.Plot_1(transl, track.num, TRUE)
  dev.off()
}


