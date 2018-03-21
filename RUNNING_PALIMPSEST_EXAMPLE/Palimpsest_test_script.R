##################################################################################################
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
##  Palimpsest - Comprehensive analysis and visualization of mutational processes in cancer genomes
##  For details on the implementation visit
##  https://github.com/FunGeST/Palimpsest
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
##################################################################################################
# Load Palimpsest library
library(Palimpsest)
#-------------------------------------------------------------------------------------------------

# Define working directories:
datadir <- "Palimpsest/RUNNING_PALIMPSEST_EXAMPLE/LiC1162/" # Path to directory containing data files
resdir <- "Results";if(!file.exists(resdir))	dir.create(resdir) # Path to directory where to export results


#-------------------------------------------------------------------------------------------------
# 1] Load genomic data and reference genome
#-------------------------------------------------------------------------------------------------
load(file.path(datadir,"mut_data.RData")) # Example mutation data form Letouzé L, Shinde J et al. (44 liver cancer genomes)
load(file.path(datadir,"cna_data.RData")) # Example CNA data from Letouzé L, Shinde J et al.
load(file.path(datadir,"annot_data.RData")) # Example Sample annotation data from Letouzé L, Shinde J et al.
load(file.path(datadir,"sv_data.RData")) # Example structural variants data from Letouzé L, Shinde J et al.

library(BSgenome.Hsapiens.UCSC.hg19) # Reference genome of choice
ref_genome <- BSgenome.Hsapiens.UCSC.hg19
load("Palimpsest/data/ensgene_hg19.RData") # Ensembl genes table (ensgene)
load("Palimpsest/data/cytoband_hg19.RData") # cytoband table (cyto)

#-------------------------------------------------------------------------------------------------
# 2] De novo mutational signature analysis
#-------------------------------------------------------------------------------------------------
resdir. <- file.path(resdir,"Signatures_DeNovo");if(!file.exists(resdir.)){dir.create(resdir.)} # Defining the results directory

vcf <- preprocessInput_snv(input_data = mut_data,ensgene=ensgene,reference_genome = ref_genome)

propMutsByCat <- palimpsestInput(vcf = vcf,type = "SNV",sample.col = "Sample",mutcat.col = "mutcat3",proportion = TRUE)
denovo_signatures <- deconvolution_nmf(input_data = propMutsByCat,type = "SNV",range_of_sigs = 2:10,nrun = 20,method = "brunet",resdir = resdir.)

# Compare with existing signatures from COSMIC database:
pdf(file.path(resdir., "Cosine_Similarity.pdf"))
cosine_similarities <- deconvolution_compare(denovo_signatures,COSMIC_Signatures)# missing lsa
dev.off()

# Define color codes for signatures
library(RColorBrewer);qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
mycol <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
mycol<- mycol[sample.int(length(mycol),nrow(denovo_signatures))];names(mycol) <- rownames(denovo_signatures)

# Calculating contributions (exposures) of signatures in each sample and generate tumor-wise graphical outputs:
signatures_exp <- deconvolution_fit(vcf=vcf,type = "SNV",input_data = propMutsByCat,threshold = 5,input_signatures = denovo_signatures,sig_cols = mycol,plot = T,resdir = resdir.)

# Plotting the exposures of signatures across the series:
pdf(file.path(resdir.,"signature_content_plot.pdf"),width=10,height=7)
signature_content_plot <- deconvolution_exposure(signatures_exp$sig_nums,signatures_exp$sig_props,sig_cols = mycol)
print(signature_content_plot)
dev.off()


#-------------------------------------------------------------------------------------------------
# 3] Extract previously known signatures (COSMIC)
#-------------------------------------------------------------------------------------------------
resdir. <- file.path(resdir,"Signatures_COSMIC");if(!file.exists(resdir.)){dir.create(resdir.)}# Defining the results directory

# restricting COSMIC signatures to mutational signatures previously identified in Liver cancers:
liver_signature_names <- c("Signature.1","Signature.4","Signature.5","Signature.6","Signature.12","Signature.16","Signature.17","Signature.22","Signature.23","Signature.24")
liver_signatures <- COSMIC_Signatures[liver_signature_names,]

# define list of colors for visualizing mutational signatures. Selecting default colors
mycol <- c("darkgreen","deepskyblue4","grey","orangered1","darkred","goldenrod1","deeppink4","royalblue4","darkolivegreen3","purple4");names(mycol) <- liver_signature_names

# Calculating contributions (exposures) of signatures in each sample:

signatures_exp <- deconvolution_fit(vcf=vcf,type = "SNV",input_data = propMutsByCat,threshold = 6,input_signatures = liver_signatures,sig_cols = mycol,plot = T,resdir = resdir.)

# Plotting the exposures of signatures across the series:
pdf(file.path(resdir.,"signature_content_plot.pdf"),width=10,height=7)
signature_content_plot <- deconvolution_exposure(signatures_exp$sig_nums,signatures_exp$sig_props,sig_cols = mycol)
print(signature_content_plot)
dev.off()


#-------------------------------------------------------------------------------------------------
# 4] Assign the probability of each individual mutation being due to each process
#-------------------------------------------------------------------------------------------------
# This step is quite computation-intensive for whole genome data. For this example we restricted the analysis to coding mutations and TERT promoter mutations
vcf.cod <- vcf[(!is.na(vcf$Gene_Name) | !is.na(vcf$Driver) & vcf$Type=="SNV"),]
vcf.cod <- palimpsestOrigin(vcf=vcf.cod,type = "SNV",sample.col = "Sample",mutcat.col = "mutcat3",signature_contribution=signatures_exp$sig_nums,input_signatures=liver_signatures)

# Estimate and represent the cumulative contribution of signatures to each driver gene
drivers <- c("CTNNB1","TP53","ARID2","NFE2L2","ACVR2A","ARID1A","AXIN1","RB1","RPS6KA3","KEAP1","ALB","CDKN2A","CDKN1A","RPL22")
matprob <- matrix(nrow=length(drivers),ncol=length(liver_signature_names),dimnames=list(drivers, liver_signature_names))
sig.cols <- grep("prob",colnames(vcf.cod))
for(i in 1:nrow(matprob)){
  g <- rownames(matprob)[i]
  ind <- which(vcf.cod$Gene_Name==g)
  matprob[i,] <- apply(vcf.cod[ind,sig.cols],2,sum,na.rm=T)
}
barplot(t(matprob),col = mycol,border = mycol,las=2)
legend("top",names(mycol),fill=mycol,ncol=5,cex=0.75,bty ="n",inset = c(0,-0.3),xpd = T)

# Compare signature 16 contribution between CTNNB1 mutations and others
library(ggplot2)
vcf.cod$CTNNB1 <- (vcf.cod$Gene_Name=="CTNNB1");vcf.cod[is.na(vcf.cod$Gene_Name),"CTNNB1"] <- FALSE
ggplot(vcf.cod, aes(x = CTNNB1,y = Signature.16.prob,color = CTNNB1,fill = CTNNB1)) + geom_violin(alpha = 0.4,adjust = 1)

# Add signature probability columns to the original vcf table
vcf <- merge(vcf,vcf.cod,all=TRUE,sort=FALSE)


#-------------------------------------------------------------------------------------------------
# 5] Clonality analysis
#-------------------------------------------------------------------------------------------------
resdir. <- file.path(resdir,"Clonality");if(!file.exists(resdir.)){dir.create(resdir.)}# Defining the results directory

# Calculate the Cancer Cell Fraction (CCF) of each mutation.
# This step is a bit long... Be patient!
vcf <- cnaCCF_annot(vcf=vcf,annot_data = annot_data,cna_data = cna_data,CCF_boundary=0.95)

# Generate graphical representations of clonality analysis
cnaCCF_plots(vcf=vcf,resdir=resdir.)


#-------------------------------------------------------------------------------------------------
# 6] Compare mutational signatures between early clonal and late subclonal mutations in each tumor
#-------------------------------------------------------------------------------------------------
resdir. <- file.path(resdir,"Signatures_early_vs_late");if(!file.exists(resdir.)){dir.create(resdir.)}# Defining the results directory

# Estimate the contribution of each signature to clonal and subclonal mutations in each tumor
vcf.clonal <- vcf[which(vcf$Clonality=="clonal"),]
propMutsByCat.clonal <- palimpsestInput(vcf = vcf.clonal,type="SNV",sample.col = "Sample", mutcat.col = "mutcat3", proportion = TRUE)
signatures_exp_clonal <- deconvolution_fit(vcf = vcf.clonal,type = "SNV",input_data = propMutsByCat.clonal,threshold = 6,input_signatures = liver_signatures,sig_cols = mycol,plot = F,resdir = resdir.)

vcf.subclonal <- vcf[which(vcf$Clonality=="subclonal"),]
propMutsByCat.subclonal <- palimpsestInput(vcf = vcf.subclonal,type="SNV",sample.col = "Sample",mutcat.col = "mutcat3",proportion = TRUE)
signatures_exp_subclonal <- deconvolution_fit(vcf = vcf.subclonal,type = "SNV",input_data = propMutsByCat.subclonal,threshold = 6,input_signatures = liver_signatures,sig_cols = mycol,plot = F,resdir = resdir.)

# Generate tumor-wise comparisons of clonal and subclonal mutations
palimpsest_DissectSigs(vcf=vcf, signatures_exp_clonal = signatures_exp_clonal, signatures_exp_subclonal = signatures_exp_subclonal,sig_cols = mycol,resdir=resdir.)

# Generate across the series comparisons of signature assigned to clonal and subclonal mutations
palimpsest_clonalitySigsCompare(clonsig = signatures_exp_clonal$sig_nums, subsig = signatures_exp_subclonal$sig_nums, msigcol = mycol, resdir = resdir.)

#-------------------------------------------------------------------------------------------------
# 7] Timing Chromosomal Gains
#-------------------------------------------------------------------------------------------------
resdir. <- file.path(resdir,"ChromosomeDups_timing");if(!file.exists(resdir.)){dir.create(resdir.)}# Defining the results directory

# Annotate vcf with chromomal gain timings
chrom_dup_time <- chrTime_annot(vcf=vcf,cna_data = cna_data,cyto=cyto)
vcf <- chrom_dup_time$vcf;point.mut.time <- chrom_dup_time$point.mut.time;cna_data <- chrom_dup_time$cna_data


# Visualizing timing plots
chrTime_plot(vcf = vcf, point.mut.time = point.mut.time, resdir = resdir.)


#-------------------------------------------------------------------------------------------------
# 8] Structural variant (SV) signature analysis:
#-------------------------------------------------------------------------------------------------
resdir. <- file.path(resdir,"SV_signatures");if(!file.exists(resdir.)){dir.create(resdir.)}# Defining the results directory

library(bedr);library(RCircos) # Loading dependencies necessary for SV annotation and CIRCOS plots
# To gain the functionality of bedr package you will need to have the BEDTools program installed and in your default PATH

# Preprocessing SV inputs and annotating for further analysis:
sv.vcf <- preprocessInput_sv(input_data =  sv_data,ensgene = ensgene,resdir = resdir.)
propSVsByCat <- palimpsestInput(vcf = sv.vcf,sample.col = "Sample",type="SV",mutcat.col = "Category",proportion = FALSE)
denovo_signatures <- deconvolution_nmf(input_data = propSVsByCat,type = "SV",range_of_sigs = 2:12,nrun =20,method = "brunet",resdir = resdir.)

# define list of colors for visualizing mutational signatures. Selecting default colors
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
mycol.sv <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
mycol.sv <- mycol.sv[sample.int(length(mycol.sv),nrow(denovo_signatures))];names(mycol.sv) <- rownames(denovo_signatures)

# Calculating contributions(exposures) of signatures in each sample:
SVsignatures_exp <- deconvolution_fit(vcf=sv.vcf,type = "SV",input_data = propSVsByCat,threshold = 6,input_signatures = denovo_signatures,sig_cols = mycol.sv,plot = T,resdir = resdir.)

# Plotting the exposures of signatures across the series:
pdf(file.path(resdir.,"signature_content_plot.pdf"),width=10,height=7)
signature_content_plot <- deconvolution_exposure(SVsignatures_exp$sig_nums,SVsignatures_exp$sig_props,sig_cols = mycol.sv)
print(signature_content_plot)
dev.off()

# Estimate the probability of each event being due to each process
sv.vcf <- palimpsestOrigin(vcf=sv.vcf, type = "SV", sample.col = "Sample", mutcat.col = "Category", signature_contribution=SVsignatures_exp$sig_nums, input_signatures=denovo_signatures)


#-------------------------------------------------------------------------------------------------
# 9] Visualize the natural history of tumor samples:
#-------------------------------------------------------------------------------------------------
resdir. <- file.path(resdir,"Natural_history");if(!file.exists(resdir.)){dir.create(resdir.)}# Defining the results directory

palimpsest_plotTumorHistories(vcf = vcf,sv.vcf = NULL, cna_data, point.mut.time, clonsig=signatures_exp_clonal$sig_props, subsig=signatures_exp_subclonal$sig_props, msigcol=mycol,msigcol.sv=mycol.sv,resdir=resdir.)


##################################################################################################
#-------------------------------------------------------------------------------------------------
# In **Introduction to Palimpsest** you can find comprehensive examples and explanations for the functions.

#-------------------------------------------------------------------------------------------------
##################################################################################################
