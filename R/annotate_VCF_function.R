#' annotate_VCF
#'
#' Function to add strand, gene and COSMIC mutation category annotations to a VCF. This function can take a long time depending on the number of mutations and how many of the annotation options you have selected. Please be patient!
#' @param vcf Input VCF to annotate.
#' @param add_strand_and_SBS_cats Logical indicating whether or not strand, gene and SBS category annotations are to be added (defaults to TRUE).
#' @param add_DBS_cats Logical indicating whether or not DBS category annotations are to be added (defaults to TRUE).
#' @param add_ID_cats Logical indicating whether or not Indel category annotations are to be added (defaults to TRUE). Unfortunately Indel mutation categories cannot be added to the VCF in Windows, as this R function calls a python script. Please run this step in a unix environment (Mac/Linux etc.).
#' @param genome_build Enter either "hg19" or "hg38".
#' @param ref_fasta File path to FASTA file compatable with input VCF positions and chromosomes.
#' @param ref_genome Name of reference genome object. For hg19 data we use the BSgenome.Hsapiens.UCSC.hg19 object, which is loaded into the local environment by library(BSgenome.Hsapiens.UCSC.hg19).
#'
#' @return vcf
#' @export
#' @import gtools
#' @import VariantAnnotation
#' @examples
#'vcf <- annotate_VCF(vcf = vcf, ref_genome = BSgenome.Hsapiens.UCSC.hg19, ref_fasta = "~/Documents/Data/Genomes/Homo_sapiens_assembly19.fasta")

annotate_VCF <- function(vcf = vcf, add_strand_and_SBS_cats = T, add_DBS_cats = T, add_ID_cats = T, genome_build = "hg19",
                         ref_fasta = NULL, ref_genome = BSgenome.Hsapiens.UCSC.hg19){
  
  if(length(colnames(vcf)[colnames(vcf) %in% c("Sample","CHROM","POS","ALT","REF")]) < 5) stop("VCF must contain columns named: 'Sample', 'CHROM', 'POS', 'REF', 'ALT' for Palimpsest functions to work")
  if(!all(unique(vcf$Type) %in% c("SNV","INS","DEL"))) stop("The column vcf$Type must contain single base substitutions marked 'SNV', deletions marked 'DEL' an/or insertions marked 'INS', please change accordingly")
  if(genome_build %!in% c("hg19","hg38")) stop("genome_build must be either hg19 or hg38")
  vcf <- order_vcf(vcf)
  
  chroms <- unique(vcf$CHROM)
  if (1 %in% chroms == TRUE)  vcf$CHROM <- paste0("chr", vcf$CHROM)
  
  remove_chrs <- c("chrM", "chrMT")
  vcf <- vcf[which(vcf$CHROM %!in% remove_chrs), ]
  if(add_strand_and_SBS_cats == T){
    vcf$strand.mut <- "+"
    vcf$strand.gene <- NA
    vcf$gene_name <- NA
    print("Adding gene, strand and SBS category annotations ..", quote = F)
    if(genome_build == "hg19") ensgene <- ensgene_hg19
    if(genome_build == "hg38") ensgene <- ensgene_hg38
    ensgene_split <- split(ensgene,ensgene$Chromosome.Name)
    vcf_split <- split(vcf,vcf$CHROM)
    chr_listy <- gtools::mixedsort(intersect(names(ensgene_split),names(vcf_split)))
    total <- length(chr_listy)
  
    for(chr in chr_listy){
      ind <- unlist(sapply(vcf_split[[chr]]$POS,function(pos){
        tmp <- which(ensgene_split[[chr]]$Gene.Start..bp. <= pos & ensgene_split[[chr]]$Gene.End..bp. >= pos)
        if(length(tmp)==1)	tmp else{NA}
      }))
      if (all(is.na(ind))) {
        vcf_split[[chr]]$gene_name <- NA 
      } else {
        vcf_split[[chr]]$gene_name <- ensgene_split[[chr]][ind,"Associated.Gene.Name"]
      }
    }
  
    vcf <- unsplit(vcf_split,vcf$CHROM)
    
    vcf$strand.gene <- c("-", NA, "+")[ensgene[match(vcf$gene_name, ensgene$Associated.Gene.Name), "Strand"] + 2]
    
    requireNamespace("VariantAnnotation", quietly = TRUE)
    nostrand <- which(is.na(vcf$strand.gene) | vcf$strand.gene == "*")
    vcf$strand.gene[nostrand] <- "+"
  
    vr <- VRanges(seqnames = vcf$CHROM, ranges = IRanges(start = vcf$POS, end = vcf$POS), ref = vcf$REF, 
                                     alt = vcf$ALT, sampleNames = vcf$Sample)
    vr@strand <- Rle(strand(vcf$strand.mut))
    delet_this <- nrow(filter(vcf, Type == "DEL"))
    vr3 <- palimpsest_addMutationContextToVR(vr = vr, ref = ref_genome, k = 3, num_of_DELs = delet_this)
    vr5 <- palimpsest_addMutationContextToVR(vr =vr, ref = ref_genome, k = 5, num_of_DELs = delet_this)

    vcf$strand.mut <- as.character(vr3@strand)
    vcf$strand.ts <- NA
  
    vcf$strand.ts[which(vcf$strand.gene == "+" & vcf$strand.mut == "-")] <- "ts"
    vcf$strand.ts[which(vcf$strand.gene == "+" & vcf$strand.mut ==  "+")] <- "nt"
    vcf$strand.ts[which(vcf$strand.gene == "-" & vcf$strand.mut  == "+")] <- "ts"
    vcf$strand.ts[which(vcf$strand.gene == "-" & vcf$strand.mut  == "-")] <- "nt"
    
    vcf$strand.gene[nostrand] <- NA
    vcf$strand.ts[nostrand] <- NA
    
    vcf$substype <- as.character(vr3$alteration)
    
    vcf$context3 <- as.character(vr3$context)
    vcf$context3[vcf$Type != "SNV"] <- NA
    vcf$SBS_cat3 <- paste(vcf$substype, vcf$context3, sep = "_")
    vcf$SBS_cat3[which(vcf$SBS_cat3 == "NA_NA")] <- NA
    vcf$SBS_cat3[vcf$Type != "SNV"] <- NA
    
    vcf$context5 <- as.character(vr5$context)
    vcf$context5[vcf$Type != "SNV"] <- NA
    vcf$SBS_cat5 <- paste(vcf$substype, vcf$context5, sep = "_")
    vcf$SBS_cat5[which(vcf$SBS_cat5 == "NA_NA")] <- NA
    vcf$SBS_cat5[vcf$Type != "SNV"] <- NA
    
    vcf <- dplyr::select(vcf,-c(substype,context3,context5))
  }
  if(add_DBS_cats == TRUE){
    vcf <- add_DBS_cats_ToVCF(vcf = vcf,DBS_mutations_only = F)
  }
  if(add_ID_cats == TRUE & .Platform$OS.type != "windows"){
    vcf <- add_ID_cats_ToVCF(vcf = vcf, ref_fasta = ref_fasta)
  }
  if(add_ID_cats == TRUE & .Platform$OS.type == "windows") warning("Unfortunately Indel mutation categories cannot be added to the VCF in Windows, as this R function calls a python script. Please run this step in a unix environment (Mac/Linux etc.). All other Palimpsest functions work on windows.")
  return(vcf)
}



#' palimpsest_addMutationContextToVR
#'
#' Function used within annotate_VCF to add mutation context to SBS mutations. 
#' @return vcf
#' @import gtools
#' @import Biostrings
#' @import GenomicRanges
#' @examples
#' vr3 <- palimpsest_addMutationContextToVR(vr = vr, ref_genome = ref_genome, k = 3, num_of_DELs = delet_this)




palimpsest_addMutationContextToVR <- function (vr =NULL, ref_genome = NULL, k = 3,  check.strand = FALSE, 
                                                  num_of_DELs = NULL) 
{
  requireNamespace("Biostrings", quietly = TRUE)
  requireNamespace("GenomicRanges", quietly = TRUE)
  if (any(width(vr)) != 1) 
    stop("SNVs must have width of 1.")
  if (k%%2 != 1) 
    stop("'k' must be odd.")
  mid = (k + 1)/2
  gr = granges(vr)
  strand.mut = strand(gr)
  if (any(strand.mut == "*")) 
    stop("The strand must be explicit in order for the correct strand to be read.")
  ranges = GenomicRanges::resize(gr, k, fix = "center")
  context = Biostrings::getSeq(ref_genome, ranges)
  ref_base = DNAStringSet(VariantAnnotation::ref(vr))
  alt_base = DNAStringSet(VariantAnnotation::alt(vr))
  
  ref0 = subseq(context, mid, mid)
  idx_invalid = (ref0 != ref_base)
  wronguns <- sum(idx_invalid) - num_of_DELs
  if (wronguns > 0) 
    warning(sprintf("References do not match in %d cases", 
                    wronguns))
  
  if (check.strand) {
    idx_minus = (strand.mut == "-")
    context[idx_minus] = Biostrings::reverseComplement(context[idx_minus])
    ref_base[idx_complement] = Biostrings::reverseComplement(ref_base[idx_complement])
    alt_base[idx_complement] = Biostrings::reverseComplement(alt_base[idx_complement])
    strand.mut[idx_minus] = "+"
  }
 
    idx_complement = as.character(ref_base) %in% c("A", 
                                                   "G")
    context[idx_complement] = Biostrings::reverseComplement(context[idx_complement])
    ref_base[idx_complement] = Biostrings::reverseComplement(ref_base[idx_complement])
    alt_base[idx_complement] = Biostrings::reverseComplement(alt_base[idx_complement])
    strand.mut[idx_complement] = "-"
  
  alteration = as.character(xscat(ref_base, alt_base))
  alteration[idx_invalid] <- NA
  context = as.character(context)
  context[idx_invalid] <- NA
  vr$alteration = alteration
  vr$context = context
  vr@strand = Rle(strand.mut)
  return(vr)
}




#' extract_dbs_from_vcf
#'
#' returns lines of a VCF corresponding to DBS mutations. The VCF must be ordered by sample, CHROM and position for this function to work (can be performed by the "order_vcf()" function.
#' @param vcf The input VCF from which DBS mutations are to be extracted
#' @keywords Signatures
#' @export
#' @examples
#' vcf_dbs <- extract_dbs_from_vcf(vcf = vcf)

extract_dbs_from_vcf <- function(vcf=NULL){
  vcf <- vcf[vcf$Type == "SNV",]
  res = as.data.frame(matrix(ncol = ncol(vcf)))
  colnames(res) = colnames(vcf)
  for(i in 1:nrow(vcf)){
    if(i==1){
      current_pos = vcf$POS[i]
      current_chr = vcf$CHROM[i]
      next
    }else{
      pos = vcf$POS[i]
      chr = vcf$CHROM[i]
      Type = vcf$Type[i]
      if(pos == current_pos+1 & chr ==current_chr & Type == "SNV"){
        res = rbind(res, vcf[(i-1):i,])
      }
      current_pos = pos
      current_chr = chr
    }
  }
  return(res[2:nrow(res),])
}



#' order_vcf
#'
#' Orders a VCF by Sample name, then genomic position within each sample and project (if project argument given).
#' @param vcf The VCF to be ordered.
#' @param Project_col Name of the Project column in the VCF (if any) (e.g. may contain "ICGC"). If a value is given, the samples are sorted by project first, then sample etc..
#' @keywords Signatures
#' @export
#' @examples
#' vcf <- order_vcf(vcf, Project_col = "Project")

order_vcf <- function(vcf, Project_col = NA){
  namecols <- colnames(vcf)
  vcf <- arrange(vcf,vcf$Sample)
  
  input_split <- split(vcf, vcf$Sample)
  
  nsamp <- length(input_split)
  
  for(i in 1:length(input_split)){
    tmp <- as.data.frame(input_split[i])
    colnames(tmp) <- namecols
    tmp <- arrange(tmp,tmp$POS)
    tmp <- arrange(tmp,tmp$CHROM)
    if(i == 1){
      res <- tmp
      next
    }else{
      res <- rbind(res,tmp)
    }
  }
  if(!is.na(Project_col)){
    res <- arrange(res,res[,Project_col])
  }
  res <- as.data.frame(res)
  return(res)
}






#' add_DBS_cats_ToVCF
#'
#' Adds DBS mutation categories to a VCF containing SNVs. N.B. The VCF must be ordered by sample, CHROM and position for this function to work.
#' @param input Input VCF to which DBS mutation categories are to be added.
#' @param DBS_mutations_only Logical, TRUE if all lines in input VCF correspond to double base substitutions, FALSE if it contains other SNVs and/or Indels.
#' @keywords Signatures
#' @import spgs
#' @examples
#' vcf <- add_DBS_cats_ToVCF(input = vcf, DBS_mutations_only = FALSE)

add_DBS_cats_ToVCF <- function(vcf = NULL, DBS_mutations_only = NA){
  requireNamespace("spgs", quietly = TRUE)
  DBS_contexts <- c("AC>CA", "AC>CG", "AC>CT", "AC>GA", "AC>GG", "AC>GT", "AC>TA", "AC>TG", "AC>TT", "AT>CA", "AT>CC", "AT>CG", "AT>GA", "AT>GC", "AT>TA", "CC>AA", "CC>AG",
                           "CC>AT", "CC>GA" ,"CC>GG" ,"CC>GT", "CC>TA", "CC>TG", "CC>TT", "CG>AA", "CG>AC", "CG>AT", "CG>GA" ,"CG>GC" ,"CG>TA", "CT>AA", "CT>AC", "CT>AG", "CT>GA",
                           "CT>GC" ,"CT>GG" ,"CT>TA" ,"CT>TC", "CT>TG", "GC>AA", "GC>AG", "GC>AT" ,"GC>CA" ,"GC>CG" ,"GC>TA" ,"TA>AC", "TA>AG", "TA>AT", "TA>CC", "TA>CG", "TA>GC",
                           "TC>AA", "TC>AG", "TC>AT", "TC>CA", "TC>CG", "TC>CT", "TC>GA", "TC>GG", "TC>GT", "TG>AA", "TG>AC", "TG>AT", "TG>CA", "TG>CC", "TG>CT", "TG>GA", "TG>GC",
                           "TG>GT", "TT>AA", "TT>AC", "TT>AG", "TT>CA", "TT>CC", "TT>CG", "TT>GA", "TT>GC", "TT>GG")

  vcf$unique <- c(1:nrow(vcf))
  add_cats = vcf %>% 
    mutate("Ref_a" = NA,"Alt_a" = NA)
  print ("Adding DBS categories..",quote = F)
  if(DBS_mutations_only == FALSE){
    add_cats <- extract_dbs_from_vcf(vcf=add_cats)
  }
  
  for (i in 1:(nrow(add_cats)-1)){
    if(add_cats$POS[i] == add_cats$POS[i+1] - 1 & add_cats$CHROM[i] == add_cats$CHROM[i+1] & add_cats$Type[i] == "SNV"){
      add_cats$Ref_a[i] <- paste0(add_cats$REF[i], add_cats$REF[i+1])
      add_cats$Alt_a[i] <- paste0(add_cats$ALT[i], add_cats$ALT[i+1])
      next
    } else {
      add_cats$Ref_a[i] <- NA; add_cats$Alt_a[i] <- NA
    }
  }
  
  for (i in 1:nrow(add_cats)){
    if(i %% 2 == 0){
      add_cats$Ref_a[i] <- add_cats$Ref_a[i-1]
      add_cats$Alt_a[i] <- add_cats$Alt_a[i-1]
    }
  }

  add_cats = add_cats %>% 
    mutate(Ref_b =reverseComplement(Ref_a,case = "as is"),
          Alt_b = reverseComplement(Alt_a,case = "as is"),
          Category_a = paste0(Ref_a, ">", Alt_a),
           Category_b = paste0(Ref_b, ">", Alt_b),
          DBS_cat = case_when(Category_a %in% DBS_contexts ~ Category_a,
                              Category_b %in% DBS_contexts ~ Category_b,
                              T ~ NA_character_))

  add_cats <- add_cats[,colnames(add_cats) %!in% c("Category_a","Category_b","Ref_a","Ref_b","Alt_a","Alt_b")]
 
   if(length(add_cats$DBS_cat[is.na(add_cats$DBS_cat)])>0) stop("Error in add_DBS_cats_ToVCF")
  
  if(DBS_mutations_only == FALSE){
    output <- vcf
    output$DBS_cat <- add_cats$DBS_cat[match(output$unique,add_cats$unique)]
  }
  if(DBS_mutations_only == TRUE) output <- add_cats
  output <- output[,colnames(output)!="unique"]
  if(nrow(filter(output, is.na(DBS_cat))) == nrow(output)) warning("No DBS mutations were detected.")

  return(output)
}









#' add_ID_cats_ToVCF
#'
#' Adds Indel mutation categories to a VCF.
#' @param vcf VCF to which Indel mutation categories are to be added.
#' @param tool_dir Path to folder containing PCAWG7-data-preparation-version-1.5 tool.
#' @param ref_fasta Path to fasta file for reference genome of choice (e.g. hg19 genome).
#' @keywords Signatures
#' @examples
#' vcf <- add_ID_cats_ToVCF(vcf = vcf,ref_fasta = ref_fasta)

add_ID_cats_ToVCF <- function(vcf = NULL, ref_fasta = NA){
  if((Sys.which("python")=="")==TRUE) stop("python must be installed on this device and accessible to R to allow indel categories to be added.
                                           (must be performed in a Unix environment)")
  if(nrow(filter(vcf, Type != "SNV")) == 0) warning("No rows of the VCF corresponding to insertions or deletions were detected.")
  palimpdir = NA
  for(i in length(.libPaths)){
    if("Palimpsest" %in% c(list.files(.libPaths()[i]))){
      palimpdir = paste0(.libPaths()[i],"/Palimpsest/")
    }
  }
  for(i in 1:1000000){
    if("exec" %!in% list.files(palimpdir)) palimpdir = NA
    if(is.na(palimpdir)){
      if(i ==1){
        print(" ", quote = F)
        print("ERROR: The Palimpsest package directory could not be located in the following default R library/libraries:", quote = F)
        print(.libPaths())
        print(" ", quote = F)
        print("This function needs the location of the up-to-date Palimpsest directory (containing the 'exec' folder) to launch the 
              indel category extraction in python, please find and enter the file path manually", quote = F)
      }
      if(i > 1){
        print(paste0("ERROR: The filepath entered is not the Palimpsest package directory downlaoded from github, please try again 
                     making sure that you are not using quotation marks"), quote = F)
      }
      print("example filepath: '/Users/joe_bloggs/Library/R/3.5/library/Palimpsest/'", quote = F)
      palimpdir = readline(prompt="Enter filepath (without quotation marks): ")
    } 
    if("exec" %in% list.files(palimpdir)) break	
  }
  tmpdir <- paste0(palimpdir,"Temporary/"); if(!file.exists(tmpdir))  dir.create(tmpdir) 
  heure <- Sys.time(); heure <- gsub(" ","_",heure); heure <- gsub("-","",heure); heure <- gsub(":",".",heure)
  
  nums <- c(1:nrow(vcf))
  vcf = vcf %>% 
    mutate(Unique = paste0(Sample,"_",nums))
  
  python_vcf = vcf %>% 
    filter(Type != "SNV") %>% 
    mutate(Start = POS, End = POS,col0="x", col3 = "y", Genome = "GRCh37",col11 = "1",col12 = "2", CHROM = sub("chr", "", CHROM)) %>% 
    order_vcf() %>% 
    dplyr::select(col0,Unique,col3,Genome,Type,CHROM,Start,End,REF,ALT,col11,col12)

  if("-" %!in% vcf$REF[vcf$Type == "INS"]){
    python_vcf[python_vcf$Type=="INS",] = python_vcf[python_vcf$Type=="INS",] %>% 
      mutate(REF = "-", ALT =  substr(ALT,2,nchar(ALT)))
  }
  python_vcf[python_vcf$Type=="INS",] = python_vcf[python_vcf$Type=="INS",] %>% mutate(End = End + 1)

  if("-" %!in% vcf$ALT[vcf$Type == "DEL"]){
    python_vcf[python_vcf$Type=="DEL",] = python_vcf[python_vcf$Type=="DEL",] %>% 
      mutate(REF =  substr(REF,2,nchar(REF)), ALT = "-", End = End + 1, Start = Start + 1)
  }
  
  vcf_output <- paste0(tmpdir,"python_vcf_indel.simple")
  
  write.table(format(python_vcf,scientific=FALSE),file=vcf_output,col.names = F,row.names=F,sep="\t",quote=F)
  
  
  ### RUN PCAWG7-data-preparation-version-1.5 IN PYTHON ###
  print("Runnng PCAWG7-data-preparation-version-1.5 in python to extract Indel categories..",quote = F)
  tool <- paste0(palimpdir,"exec/make_spectra_indels.py")

  cachedir <- paste0(tmpdir,"cache_",heure,"/");if(!file.exists(cachedir))  dir.create(cachedir) 
  
  argus <- paste("--cachedir",cachedir,"--genome --fasta",ref_fasta,"--output",paste0(tmpdir,"indel_output"),vcf_output)

  # system2("python",c(tool,argus))
  system(paste(tool,c(argus)))
  warning("Indel category extraction with PCAWG7-data-preparation-version-1.5 python script is finished 
  (if there are error messages above it has not been successful)")
  

  # load output and add categories to vcf
  cachename <- list.files(paste0(cachedir,"annotated/"))
  indel_cache <- read.delim(file = paste0(cachedir,"annotated/",cachename), sep = "\t", header = F) %>% 
    mutate_all(as.character())
  
  vcf$ID_cat <- as.character(indel_cache$V7[match(vcf$Unique,indel_cache$V1)])
  
  output <- indel_category_correct(input = vcf) %>% 
    dplyr::select(-c(Unique))
  
  return(output)
  unlink(tmpdir,recursive = T,force = T) ## delete the temporary folder and files
}




#' indel_category_correct
#'
#' Corrects Indel category output from PCAWG7-data-preparation-version-1.5 in the add_ID_cats_ToVCF function. Palimpsest uses the COSMIC indel categories (taken from the Sanger Institute's SigProfiler tool), whereas the PCAWG7-data-preparation-version-1.5 (from the Broad institute) uses slightly different categories.
#' @param input VCF of indel mutations
#' @keywords Signatures
#' @examples
#' vcf_ID <- indel_category_correct(input = vcf)

indel_category_correct <- function(input = NULL){
  
  true_types = c("DEL_C_1_0","DEL_C_1_1","DEL_C_1_2","DEL_C_1_3","DEL_C_1_4","DEL_C_1_5+", "DEL_T_1_0","DEL_T_1_1","DEL_T_1_2",
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
  
  for( i in 1:nrow(input)){
    if(!is.na(input$ID_cat[i])){
      if( grepl("DEL,C,1,",input$ID_cat[i]) && substr(input$ID_cat[i], 9,9) %in% (1:5) ){
        numero <- as.numeric(substr(input$ID_cat[i], 9,9))
        input$ID_cat[i] <- paste0("DEL,C,1,",(numero-1))
      }
      if( input$ID_cat[i] ==  "DEL,C,1,6+"){
        input$ID_cat[i] <- paste("DEL,C,1,5+")
      }
      
      if( grepl("DEL,T,1,",input$ID_cat[i]) && substr(input$ID_cat[i], 9,9) %in% (1:5) ){
        numero <- as.numeric(substr(input$ID_cat[i], 9,9))
        input$ID_cat[i] <- paste0("DEL,T,1,",(numero-1))
      }
      if( input$ID_cat[i] ==  "DEL,T,1,6+"){
        input$ID_cat[i] <- paste("DEL,T,1,5+")
      }
      
      if( grepl("DEL,repeats,",input$ID_cat[i]) && substr(input$ID_cat[i], 15,15) %in% (1:5) ){
        numero <- as.numeric(substr(input$ID_cat[i], 15,15))
        input$ID_cat[i] <- paste0(substr(input$ID_cat[i],1,14),(numero-1))
      }
      if( grepl("DEL,repeats,",input$ID_cat[i]) && substr(input$ID_cat[i], 15,15) == 6 ){
        input$ID_cat[i] <- paste0(substr(input$ID_cat[i],1,14),"5+")
      }
      if( grepl("DEL,repeats,",input$ID_cat[i]) && substr(input$ID_cat[i], 16,16) %in% (1:5) ){
        numero <- as.numeric(substr(input$ID_cat[i], 16,16))
        input$ID_cat[i] <- paste0(substr(input$ID_cat[i],1,15),(numero-1))
      }
      if( input$ID_cat[i] ==  "DEL,repeats,5+,6+"){
        input$ID_cat[i] <- paste("DEL,repeats,5+,5+")
      }
    }
    
  }
  if(length(setdiff(input$ID_cat[!is.na(input$ID_cat)],gsub("_",",",true_types)))>0) stop("Error in correcting indel categories from PCAWG7-data-preparation-version-1.5 output")
  return(input)
}

