##################################################################################################
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
##  Palimpsest 2.0 - Comprehensive analysis and visualisation of mutational processes in cancer genomes
##  For more details on the implementation visit
##  https://github.com/FunGeST/Palimpsest
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
##################################################################################################
# Load Palimpsest package
library(Palimpsest)
library(BSgenome.Hsapiens.UCSC.hg19) # Reference genome of choice

#-------------------------------------------------------------------------------------------------

# define input directory containing example dataset
# D/L link: https://github.com/FunGeST/Palimpsest/tree/master/RUNNING_PALIMPSEST_EXAMPLE/LiC1162
datadir <- "/Palimpsest/RUNNING_PALIMPSEST_EXAMPLE/LiC1162/" 

# define parent results directory 
resdir_parent <- "/Results/"
if(!file.exists(resdir_parent)) dir.create(resdir_parent)

#-------------------------------------------------------------------------------------------------
# 1] Load and annotate mutation data (VCF)
#-------------------------------------------------------------------------------------------------

vcf = load2object(paste0(datadir,"vcf.RData"))

vcf <- annotate_VCF(vcf = vcf, ref_genome = BSgenome.Hsapiens.UCSC.hg19,
                    ref_fasta = "/Genomes/Homo_sapiens_assembly19.fasta")


#-------------------------------------------------------------------------------------------------
# 2] Single base substitution (SBS) de novo mutational signature extraction  
#    (N.B. these de novo extraction functions work for DBS, Indel & SV too!)
#-------------------------------------------------------------------------------------------------
resdir <- paste0(resdir_parent,"SBS_denovo/"); if(!file.exists(resdir))	dir.create(resdir) ## set results directory

# Extract de novo SBS signatures
SBS_input <- palimpsest_input(vcf = vcf, Type = "SBS")
SBS_denovo_sigs <- NMF_Extraction(input_matrices = SBS_input, range_of_sigs = 2:10, nrun = 10,
                                 resdir = resdir)

# Compare the de novo signatures with published COSMIC signatures
compare_tab <- compare_results(reference_sigs = SBS_cosmic,extraction_1 = SBS_denovo_sigs); compare_tab
readr::write_delim(compare_tab,path = paste0(resdir,"Comparison_table.txt"))

pdf(file.path(resdir, "Cosine_Similarity_Heatmap.pdf"), width = 11, height = 10)
SBS_cosine_similarities <- deconvolution_compare(SBS_denovo_sigs,SBS_cosmic)
dev.off()
save(SBS_cosine_similarities, file = paste0(resdir,"Cosine_Similarity_matrix.RData"))


# Define signature colours for plotting
SBS_col <- signature_colour_generator(rownames(SBS_denovo_sigs))

# Calculate and plot the exposure of the signatures across the series
SBS_signatures_exp <- deconvolution_fit(input_matrices = SBS_input,input_signatures = SBS_denovo_sigs,
                                        threshold = 6,resdir = resdir, signature_colours = SBS_col,
                                        input_vcf = vcf)

pdf(file.path(resdir,"signature_content_plot.pdf"),width=15,height=10)
deconvolution_exposure(signature_colours = SBS_col, signature_contribution = SBS_signatures_exp)
dev.off()


#-------------------------------------------------------------------------------------------------
# 3] Extract with published SBS COSMIC signatures
#-------------------------------------------------------------------------------------------------
resdir <- paste0(resdir_parent,"SBS_COSMIC_Extraction/");if(!file.exists(resdir))	dir.create(resdir) 

# select desired COSMIC SBS reference signatures 
SBS_liver_names <- c("SBS1","SBS4","SBS5","SBS6","SBS12","SBS16","SBS17","SBS18","SBS22","SBS23",
                     "SBS24")
for(new_name in c(SBS_liver_names[SBS_liver_names %!in% names(sig_cols)]))
  sig_cols[new_name] <- signature_colour_generator(new_name)  ## generate colours for new signatures 

SBS_liver_sigs <- SBS_cosmic[rownames(SBS_cosmic) %in% SBS_liver_names,]

# calculate and plot the exposure of the signatures across the series
SBS_signatures_exp <- deconvolution_fit(input_matrices = SBS_input, input_signatures = SBS_liver_sigs,
                                        threshold = 6, signature_colours = sig_cols,resdir = resdir, 
                                        input_vcf = vcf)

pdf(file.path(resdir,"SBS_signature_content_plot.pdf"),width=15,height=10)
deconvolution_exposure(signature_contribution = SBS_signatures_exp,signature_colours = sig_cols)
dev.off()


#-------------------------------------------------------------------------------------------------
# 4] Extract with published double base substitution (DBS) COSMIC signatures
#-------------------------------------------------------------------------------------------------
resdir <- paste0(resdir_parent,"DBS_COSMIC_Extraction/");if(!file.exists(resdir))	dir.create(resdir) 

DBS_input <- palimpsest_input(vcf = vcf, Type = "DBS")

# select desired COSMIC DBS reference signatures 
DBS_liver_names <- c("DBS2","DBS4","DBS5","DBS7","DBS11")
for(new_name in c(DBS_liver_names[DBS_liver_names %!in% names(sig_cols)])) 
  sig_cols[new_name] <- signature_colour_generator(new_name)  ## generate colours for signatures currently without defined colours 

DBS_liver_sigs <- DBS_cosmic[rownames(DBS_cosmic) %in% DBS_liver_names,]

# calculate and plot the exposure of the signatures across the series
DBS_signatures_exp <- deconvolution_fit(input_matrices = DBS_input, input_signatures = DBS_liver_sigs, 
                                        threshold = 6, signature_colours = sig_cols,resdir = resdir)

pdf(file.path(resdir,"DBS_signature_content_plot.pdf"),width=15,height=10)
deconvolution_exposure(signature_contribution = DBS_signatures_exp,signature_colours = sig_cols)
dev.off()


#-------------------------------------------------------------------------------------------------
# 5] Extract with published insertion and deletion (ID) COSMIC signatures
#-------------------------------------------------------------------------------------------------
resdir <- paste0(resdir_parent,"ID_COSMIC_Extraction/");if(!file.exists(resdir))	dir.create(resdir) 

ID_input <- palimpsest_input(vcf = vcf, Type = "ID")

# select desired COSMIC indel reference signatures 
ID_liver_names <- c("ID1","ID2","ID3","ID5","ID8")
for(new_name in c(ID_liver_names[ID_liver_names %!in% names(sig_cols)])) 
  sig_cols[new_name] <- signature_colour_generator(new_name)  ## generate colours for signatures currently without defined colours 

ID_liver_sigs <- ID_cosmic[rownames(ID_cosmic) %in% ID_liver_names,]

# calculate and plot the exposure of the signatures across the series
ID_signatures_exp <- deconvolution_fit(input_matrices = ID_input, input_signatures = ID_liver_sigs, 
                                       threshold = 6, signature_colours = sig_cols,resdir = resdir)

pdf(file.path(resdir,"ID_signature_content_plot.pdf"),width=15,height=10)
deconvolution_exposure(signature_contribution = ID_signatures_exp,signature_colours = sig_cols)
dev.off()


#-------------------------------------------------------------------------------------------------
# 6] Assign the probability of each individual mutation being due to each process 
# (Again, this function works for DBS, Indel and SV signatures too)
#-------------------------------------------------------------------------------------------------
# This step is quite computationally-intensive for whole genome data. 
# For this example we restrict the analysis to coding mutations 
vcf.cod <- vcf[(!is.na(vcf$gene_name) & vcf$Type=="SNV"),]
vcf.cod <- signature_origins(input = vcf.cod, Type = "SBS",input_signatures = SBS_liver_sigs,
                             signature_contribution = SBS_signatures_exp)

# Estimate and represent the cumulative contribution of signatures to each driver gene
drivers <- c("CTNNB1","TP53","ARID2","NFE2L2","ACVR2A","ARID1A","AXIN1","RB1","RPS6KA3","KEAP1",
             "ALB","CDKN2A","CDKN1A","RPL22")
matprob <- matrix(nrow=length(drivers),ncol=length(SBS_liver_names),dimnames=list(drivers, SBS_liver_names))
sig.cols <- grep("prob",colnames(vcf.cod))
for(i in 1:nrow(matprob)){
  g <- rownames(matprob)[i]
  ind <- which(vcf.cod$gene_name==g)
  matprob[i,] <- apply(vcf.cod[ind,sig.cols],2,sum,na.rm=T)
}
barplot(t(matprob),col = sig_cols,border = sig_cols,las=2)
legend("top",names(sig_cols)[names(sig_cols) %in% rownames(SBS_liver_sigs)],fill=sig_cols,ncol=5,
       cex=0.75,bty ="n",inset = c(0,-0.3),xpd = T)

# Compare signature 16 contribution between CTNNB1 mutations and others
library(ggplot2)
vcf.cod$CTNNB1 <- (vcf.cod$gene_name=="CTNNB1");vcf.cod[is.na(vcf.cod$gene_name),"CTNNB1"] <- FALSE
ggplot(vcf.cod, aes(x = CTNNB1,y = SBS5.prob,color = CTNNB1,fill = CTNNB1)) + geom_violin(alpha = 0.4,adjust = 1)
wilcox.test(c(vcf.cod$SBS5.prob) ~ c(vcf.cod$CTNNB1))

# Add signature probability columns to the original vcf
vcf <- merge(vcf,vcf.cod,all=TRUE,sort=FALSE)


#-------------------------------------------------------------------------------------------------
# 7] Clonality analysis
#-------------------------------------------------------------------------------------------------
resdir <- file.path(resdir_parent,"Clonality");if(!file.exists(resdir)){dir.create(resdir)}

# Load copy number analysis (CNA) and annotation (annot) data
cna_data = load2object(paste0(datadir,"cna_data.RData"))
annot = load2object(paste0(datadir,"annot_data.RData"))

# Calculate the Cancer Cell Fraction (CCF) of each mutation.
vcf_cna <- cnaCCF_annot(vcf= vcf, annot_data = annot,cna_data = cna_data, CCF_boundary = 0.95)

# Generate graphical representations of clonality analysis
cnaCCF_plots(vcf=vcf_cna,resdir=resdir)


#-------------------------------------------------------------------------------------------------
# 8] Compare mutational signatures between early clonal and late subclonal mutations in each tumour
#-------------------------------------------------------------------------------------------------
resdir <- file.path(resdir_parent,"Signatures_early_vs_late/");if(!file.exists(resdir)){dir.create(resdir)}

# Estimate the contribution of each signature to clonal and subclonal mutations in each tumour
vcf.clonal <- vcf_cna[which(vcf_cna$Clonality=="clonal"),]
SBS_input_clonal <- palimpsest_input(vcf = vcf.clonal,Type = "SBS")
signatures_exp_clonal <- deconvolution_fit(input_matrices = SBS_input_clonal,
                                           input_signatures = SBS_liver_sigs, resdir =  resdir,save_signatures_exp = F)

vcf.subclonal <- vcf_cna[which(vcf_cna$Clonality=="subclonal"),]
SBS_input_subclonal <- palimpsest_input(vcf = vcf.subclonal,Type = "SBS")
signatures_exp_subclonal <- deconvolution_fit(input_matrices = SBS_input_subclonal,
                                              input_signatures = SBS_liver_sigs, resdir =  resdir,save_signatures_exp = F)

# Generate per tumour comparisons of clonal and subclonal mutations
palimpsest_DissectSigs(vcf=vcf_cna, signatures_exp_clonal = signatures_exp_clonal,
                       signatures_exp_subclonal = signatures_exp_subclonal,sig_cols = sig_cols,resdir=resdir)

# Generate across the series comparisons of signature assigned to clonal and subclonal mutations
palimpsest_clonalitySigsCompare(clonsig = signatures_exp_clonal$sig_nums,
                                subsig = signatures_exp_subclonal$sig_nums, msigcol = sig_cols, resdir = resdir)
dev.off()


#-------------------------------------------------------------------------------------------------
# 8] Timing Chromosomal Gains
#-------------------------------------------------------------------------------------------------
resdir <- file.path(resdir_parent,"ChromosomeDups_timing/");if(!file.exists(resdir)){dir.create(resdir)}

# Annotate vcf with chromomal gain timings
chrom_dup_time <- chrTime_annot(vcf=vcf_cna,cna_data = cna_data,cyto=cytoband_hg19)
vcf_cna <- chrom_dup_time$vcf;point.mut.time <- chrom_dup_time$point.mut.time;cna_data <- chrom_dup_time$cna_data

# Visualising timing plots
chrTime_plot(vcf = vcf_cna, point.mut.time = point.mut.time, resdir = resdir,cyto = cytoband_hg19)


#-------------------------------------------------------------------------------------------------
# 9] Structural variant (SV) signature analysis:
#-------------------------------------------------------------------------------------------------
resdir <- file.path(resdir_parent,"SV_signatures/");if(!file.exists(resdir)){dir.create(resdir)}

library(bedr);library(RCircos) # Loading dependencies necessary for SV annotation and CIRCOS plots
# For the bedr package to function you will need to have the BEDTools programme installed and in your default PATH

# Load SV Data
SV_data = load2object(paste0(datadir,"sv_data.RData"))

# Preprocess SV inputs and annotate for further analysis:
SV.vcf <- preprocessInput_sv(input_data =  SV_data,resdir = resdir)
SV_input <- palimpsest_input(vcf = SV.vcf,Type = "SV")
SV_denovo_sigs <- NMF_Extraction(input_matrices =  SV_input,range_of_sigs = 1:10,nrun = 10,
                                 resdir = resdir)

# define list of colors for visualizing mutational signatures. Selecting default colors
SV_cols <- signature_colour_generator(rownames(SV_denovo_sigs))

# Calculate contribution of signatures in each sample:
SVsignatures_exp <- deconvolution_fit_SV(vcf = SV.vcf,input_data = SV_input$mut_props,
                                         input_signatures = SV_denovo_sigs,sig_cols = SV_cols,
                                         resdir = resdir)

# Plotting the exposures of signatures across the series:
pdf(file.path(resdir,"signature_content_plot.pdf"),width=10,height=7)
deconvolution_exposure(signature_contribution = SVsignatures_exp,signature_colours = SV_cols)
dev.off()

# Estimate the probability of each event being due to each process
SV.vcf <- signature_origins(input = SV.vcf,Type = "SV",signature_contribution = SVsignatures_exp,
                            input_signatures = SV_denovo_sigs)


#-------------------------------------------------------------------------------------------------
# 10] Visualise the natural history of tumour samples:
#-------------------------------------------------------------------------------------------------
resdir <- file.path(resdir_parent,"Natural_history/");if(!file.exists(resdir)){dir.create(resdir)}

palimpsest_plotTumorHistories(vcf = vcf_cna, sv.vcf = SV.vcf,cna_data =  cna_data, 
                              point.mut.time = point.mut.time, 
                              clonsig = signatures_exp_clonal$sig_props, 
                              subsig = signatures_exp_subclonal$sig_props, 
                              msigcol = sig_cols, msigcol.sv = SV_cols, resdir = resdir)


##################################################################################################
#-------------------------------------------------------------------------------------------------
# In [**Introduction to Palimpsest** you can find comprehensive examples and explanations for the functions.
# Link: https://github.com/FunGeST/Palimpsest/tree/master/Files/vignette_palimpsest.pdf
#-------------------------------------------------------------------------------------------------
##################################################################################################
