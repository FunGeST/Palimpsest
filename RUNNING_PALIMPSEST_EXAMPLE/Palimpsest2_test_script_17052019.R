##################################################################################################
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
##  Palimpsest 2.0 - Comprehensive analysis and visualization of mutational processes in cancer genomes
##  For details on the implementation visit
##  https://github.com/FunGeST/Palimpsest
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
##################################################################################################
# Load Palimpsest library
library(Palimpsest)
library(BSgenome.Hsapiens.UCSC.hg19) # Reference genome of choice

#-------------------------------------------------------------------------------------------------

# Define working directories:
palimpdir <- "~/Documents/Code/Palimpsest2/"  ### Path Palimpsest master directory, folder structure must remain unchanged for Palimpsest to function.
resdir_parent <- paste0("~/Documents/Results/delet/Palimpsest2_test_script_",gsub("-","",Sys.Date()),"/");if(!file.exists(resdir_parent))	dir.create(resdir_parent) ## Path to parent results directory


#-------------------------------------------------------------------------------------------------
# 1] Load and annotate mutation data
#-------------------------------------------------------------------------------------------------

vcf <- load2object("/Volumes/HEPATO PARTAGE/GENOMIC_DATA/WGS/RData/Somatic_mutations/140s_113p_Oncotated_Lifted_Somatic_variants_filtered_main_columns_fully_annotated.RData") %>% 
  dplyr::select(Sample, CHROM, POS = POS_hg19, REF, ALT, Type = Variant_Type) %>% 
  filter(Sample %in% unique(Sample)[1:30]) %>% 
  mutate(Type = ifelse(Type == "SNP","SNV",Type))

vcf <- annotate_VCF(vcf = vcf, ref_genome = BSgenome.Hsapiens.UCSC.hg19,palimpdir = palimpdir,
                    ref_fasta = "~/Documents/Data/Genomes/Homo_sapiens_assembly19.fasta")


#-------------------------------------------------------------------------------------------------
# 2] Single base substitution (SBS) de novo mutational signature extraction  (N.B. de novo extraction functions work for DBS & Indel too!)
#-------------------------------------------------------------------------------------------------
resdir <- paste0(resdir_parent,"SBS_denovo/"); if(!file.exists(resdir))	dir.create(resdir) ## set results directory

# Extract de novo SBS signatures
SBS_input <- Palimpsest_input2(vcf = vcf, Type = "SBS")
SBS_denovo_sigs <- NMF_Extraction(input_matrices = SBS_input,Type = "SBS",range_of_sigs = 1:10,nrun = 10,resdir = resdir)

# Compare the de novo signatures with published COSMIC signatures
compare_tab <- compare_results(reference_sigs = SBS_cosmic,extraction_1 = SBS_denovo_sigs)
readr::write_delim(compare_tab,path = paste0(resdir,"Comparison_table.txt"))

pdf(file.path(resdir, "Cosine_Similarity_Heatmap.pdf"), width = 11, height = 10)
SBS_cosine_similarities <- deconvolution_compare2(SBS_denovo_sigs,SBS_cosmic)
dev.off()
save(SBS_cosine_similarities, file = paste0(resdir,"Cosine_Similarity_matrix.RData"))


# Define signature colours for plotting
SBS_col <- signature_colour_generator(rownames(SBS_denovo_sigs))

# Calculate and plot the exposure of the signatures across the series
SBS_signatures_exp <- deconvolution_fit2(Type = "SBS",input_matrixes = SBS_input,input_signatures = SBS_denovo_sigs,threshold = 6,resdir = resdir, sig_cols = SBS_col)

pdf(file.path(resdir,"signature_content_plot.pdf"),width=15,height=10)
deconvolution_exposure_2(sig_cols = SBS_col,signature_contribution = SBS_signatures_exp)
dev.off()


#-------------------------------------------------------------------------------------------------
# 3] Extract with published SBS COSMIC signatures
#-------------------------------------------------------------------------------------------------
resdir <- paste0(resdir_parent,"SBS_COSMIC_Extraction/");if(!file.exists(resdir))	dir.create(resdir) 

# select desired COSMIC SBS reference signatures 
SBS_liver_names <- c("SBS1","SBS4","SBS5","SBS6","SBS12","SBS16","SBS17","SBS18","SBS22","SBS23","SBS24","SBS35")
for(new_name in c(SBS_liver_names[SBS_liver_names %!in% names(sig_cols)])) sig_cols[new_name] <- signature_colour_generator(new_name)  ## generate colours for signatures currently without defined colours 

SBS_liver_sigs <- SBS_cosmic[rownames(SBS_cosmic) %in% SBS_liver_names,]

# calculate and plot the exposure of the signatures across the series
SBS_signatures_exp <- deconvolution_fit2(Type = "SBS",input_matrixes = SBS_input, input_signatures = SBS_liver_sigs, threshold = 6, sig_cols = sig_cols,resdir = resdir)

pdf(file.path(resdir,"SBS_signature_content_plot.pdf"),width=15,height=10)
deconvolution_exposure_2(signature_contribution = SBS_signatures_exp,sig_cols = sig_cols)
dev.off()


#-------------------------------------------------------------------------------------------------
# 4] Extract with published double base substitution (DBS) COSMIC signatures
#-------------------------------------------------------------------------------------------------
resdir <- paste0(resdir_parent,"DBS_COSMIC_Extraction/");if(!file.exists(resdir))	dir.create(resdir) 

DBS_input <- Palimpsest_input2(vcf = vcf, Type = "DBS")

# select desired COSMIC DBS reference signatures 
DBS_liver_names <- c("DBS2","DBS4","DBS5","DBS7","DBS11")
for(new_name in c(DBS_liver_names[DBS_liver_names %!in% names(sig_cols)])) sig_cols[new_name] <- signature_colour_generator(new_name)  ## generate colours for signatures currently without defined colours 

DBS_liver_sigs <- DBS_cosmic[rownames(DBS_cosmic) %in% DBS_liver_names,]

# calculate and plot the exposure of the signatures across the series
DBS_signatures_exp <- deconvolution_fit2(Type = "DBS",input_matrixes = DBS_input, input_signatures = DBS_liver_sigs, threshold = 6, sig_cols = sig_cols,resdir = resdir)

pdf(file.path(resdir,"DBS_signature_content_plot.pdf"),width=15,height=10)
deconvolution_exposure_2(signature_contribution = DBS_signatures_exp,sig_cols = sig_cols)
dev.off()


#-------------------------------------------------------------------------------------------------
# 5] Extract with published insertion and deletion (ID) COSMIC signatures
#-------------------------------------------------------------------------------------------------
resdir <- paste0(resdir_parent,"ID_COSMIC_Extraction/");if(!file.exists(resdir))	dir.create(resdir) 

ID_input <- Palimpsest_input2(vcf = vcf, Type = "ID")

# select desired COSMIC indel reference signatures 
ID_liver_names <- c("ID1","ID2","ID3","ID5","ID8")
for(new_name in c(ID_liver_names[ID_liver_names %!in% names(sig_cols)])) sig_cols[new_name] <- signature_colour_generator(new_name)  ## generate colours for signatures currently without defined colours 

ID_liver_sigs <- ID_cosmic[rownames(ID_cosmic) %in% ID_liver_names,]

# calculate and plot the exposure of the signatures across the series
ID_signatures_exp <- deconvolution_fit2(Type = "ID",input_matrixes = ID_input, input_signatures = ID_liver_sigs, threshold = 6, sig_cols = sig_cols,resdir = resdir)

pdf(file.path(resdir,"ID_signature_content_plot.pdf"),width=15,height=10)
deconvolution_exposure_2(signature_contribution = ID_signatures_exp,sig_cols = sig_cols)
dev.off()


#-------------------------------------------------------------------------------------------------
# 6] Assign the probability of each individual mutation being due to each process
#-------------------------------------------------------------------------------------------------
# This step is quite computationally-intensive for whole genome data. For this example we restrict the analysis to coding mutations and TERT promoter mutations
vcf.cod <- vcf[(!is.na(vcf$gene_name) & vcf$Type=="SNV"),]
vcf.cod <- signature_origins(input = vcf.cod, Type = "SBS",input_signatures = SBS_liver_sigs,signature_contribution = SBS_signatures_exp)

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
# 7] Clonality analysis
#-------------------------------------------------------------------------------------------------
resdir. <- file.path(resdir,"Clonality");if(!file.exists(resdir.)){dir.create(resdir.)}# Defining the results directory

# Calculate the Cancer Cell Fraction (CCF) of each mutation.
# This step is a bit long... Be patient!
vcf <- cnaCCF_annot(vcf=vcf,annot_data = annot_data,cna_data = cna_data,CCF_boundary=0.95)

# Generate graphical representations of clonality analysis
cnaCCF_plots(vcf=vcf,resdir=resdir.)


#-------------------------------------------------------------------------------------------------
# 8] Compare mutational signatures between early clonal and late subclonal mutations in each tumor
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
# 8] Timing Chromosomal Gains
#-------------------------------------------------------------------------------------------------
resdir. <- file.path(resdir,"ChromosomeDups_timing");if(!file.exists(resdir.)){dir.create(resdir.)}# Defining the results directory

# Annotate vcf with chromomal gain timings
chrom_dup_time <- chrTime_annot(vcf=vcf,cna_data = cna_data,cyto=cyto)
vcf <- chrom_dup_time$vcf;point.mut.time <- chrom_dup_time$point.mut.time;cna_data <- chrom_dup_time$cna_data


# Visualizing timing plots
chrTime_plot(vcf = vcf, point.mut.time = point.mut.time, resdir = resdir.)


#-------------------------------------------------------------------------------------------------
# 9] Structural variant (SV) signature analysis:
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
# 10] Visualize the natural history of tumor samples:
#-------------------------------------------------------------------------------------------------
resdir. <- file.path(resdir,"Natural_history");if(!file.exists(resdir.)){dir.create(resdir.)}# Defining the results directory

palimpsest_plotTumorHistories(vcf = vcf,sv.vcf = NULL, cna_data, point.mut.time, clonsig=signatures_exp_clonal$sig_props, subsig=signatures_exp_subclonal$sig_props, msigcol=mycol,msigcol.sv=mycol.sv,resdir=resdir.)


##################################################################################################
#-------------------------------------------------------------------------------------------------
# In **Introduction to Palimpsest** you can find comprehensive examples and explanations for the functions.

#-------------------------------------------------------------------------------------------------
##################################################################################################
