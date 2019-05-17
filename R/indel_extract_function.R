#' add_ID_cats_ToVCF
#'
#' Adds Indel mutation categories to a VCF.
#' @param vcf VCF to which Indel mutation categories are to be added.
#' @param tool_dir Path to folder containing PCAWG7-data-preparation-version-1.5 tool.
#' @param tmpdir Path to folder where temporary files are stored.
#' @param ref_fasta Path to fasta file for reference genome of choice (e.g. hg19 genome).
#' @keywords Signatures
#' @examples
#' vcf <- add_ID_cats_ToVCF(vcf = vcf,palimpdir = palimpdir, ref_fasta = ref_fasta)

add_ID_cats_ToVCF <- function(vcf = NULL,  palimpdir = palimpdir, ref_fasta = NA){
  if((Sys.which("python")=="")==TRUE) stop("python must be installed on this device and accessible to R for this function to work")
  tmpdir <- paste0(palimpdir,"Temporary/"); if(!file.exists(tmpdir))	dir.create(tmpdir) 
  heure <- Sys.time(); heure <- gsub(" ","_",heure); heure <- gsub("-","",heure); heure <- gsub(":",".",heure)
  nums <- c(1:nrow(vcf))
  vcf = vcf %>% 
    mutate(Unique = paste0(Sample,"_",nums))
  
  python_vcf = vcf %>% 
    filter(Type != "SNV") %>% 
    mutate(Start = POS +1, End = POS+1,col0="x", col3 = "y", Genome = "GRCh37",col11 = "1",col12 = "2", CHROM = sub("chr", "", CHROM)) %>% 
    order_vcf() %>% 
    dplyr::select(col0,Unique,col3,Genome,Type,CHROM,Start,End,REF,ALT,col11,col12)
  
  python_vcf$Start[python_vcf$Type=="INS"] = python_vcf$Start[python_vcf$Type=="INS"] %>% -1
  
  python_vcf[python_vcf$Type=="INS",] = python_vcf[python_vcf$Type=="INS",] %>% 
    mutate(REF = "-", ALT =  substr(ALT,2,nchar(ALT)))
  
  python_vcf[python_vcf$Type=="DEL",] = python_vcf[python_vcf$Type=="DEL",] %>% 
    mutate(REF =  substr(REF,2,nchar(REF)), ALT = "-")
  
  vcf_output <- paste0(tmpdir,"python_vcf_indel.simple")
  
  write.table(python_vcf,file=vcf_output,col.names = F,row.names=F,sep="\t",quote=F)
  
  
  ### RUN PCAWG7-data-preparation-version-1.5 IN PYTHON ###
  print("Runnng PCAWG7-data-preparation-version-1.5 in python to extract Indel categories..",quote = F)
  tool_dir <- paste0(palimpdir, "Indel_category_extraction/PCAWG7-data-preparation-version-1.5_for_denovo/")
  tool <- paste0(tool_dir,"make_spectra_indels.py")
  cachedir <- paste0(tmpdir,"cache_",heure,"/");if(!file.exists(cachedir))	dir.create(cachedir) 
  
  argus <- paste("--cachedir",cachedir,"--genome --fasta",ref_fasta,"--output",paste0(tmpdir,"indel_output"),vcf_output)
  system2("python",c(tool,argus))
  print("Indel category extraction complete (if there are error messages above it has not been successful)",quote = F)
  

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




