#' annotate_VCF
#'
#' Function to add strand, gene and COSMIC mutation category annotations to a VCF. This function can take a long time depending on the number of mutations and how many of the annotation options you have selected. Please be patient!
#' @param vcf Input VCF to annotate.
#' @param add_strand_and_SBS_cats Logical indicating whether or not strand, gene and SBS category annotations are to be added.
#' @param add_DBS_cats Logical indicating whether or not DBS category annotations are to be added.
#' @param add_ID_cats Logical indicating whether or not Indel category annotations are to be added.
#' @param genome_build Enter either "hg19" or "hg38".
#' @param palimpdir File path to Palimpsest master folder.
#' @param ref_fasta File path to FASTA file compatable with input VCF positions and chromosomes.
#' @param ref_genome Name of reference genome object.
#'
#' @return vcf
#' @export
#' @import gtools
#' @import VariantAnnotation
#' @examples
#'vcf <- annotate_VCF(vcf = vcf, ref_genome = BSgenome.Hsapiens.UCSC.hg19, ref_fasta = "~/Documents/Data/Genomes/Homo_sapiens_assembly19.fasta")

annotate_VCF <- function(vcf = vcf, add_strand_and_SBS_cats = T, add_DBS_cats = T, add_ID_cats = T, genome_build = "hg19",
                         palimpdir = NULL, ref_fasta = NULL, ref_genome = BSgenome.Hsapiens.UCSC.hg19){
  
  if(length(colnames(vcf)[colnames(vcf) %in% c("Sample","CHROM","POS","ALT","REF")]) < 5) stop("VCF must contain columns named: 'Sample', 'CHROM', 'POS', 'REF', 'ALT' for Palimpsest functions to work")
  vcf <- order_vcf(vcf)
  
  chroms <- unique(vcf$CHROM)
  if (1 %in% chroms == TRUE)  vcf$CHROM <- paste("chr", vcf$CHROM, sep = "")
  
  
  remove_chrs <- c("chrM", "chrMT")
  vcf <- vcf[which(vcf$CHROM %!in% remove_chrs), ]
  if(add_strand_and_SBS_cats == T){
    vcf$strand.mut <- "+"
    vcf$strand.gene <- NA
    vcf$gene_name <- NA
    print("Adding gene, strand and SBS category annotations ..", quote = F)
    if(genome_build == "hg19") ensgene <- Palimpsest:::ensgene_hg19
    if(genome_build == "hg38") ensgene <- Palimpsest:::ensgene_hg38
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
    vr3 <- palimpsest_addMutationContextToVR_2(vr = vr, ref = ref_genome, k = 3, num_of_DELs = delet_this)
    vr5 <- palimpsest_addMutationContextToVR_2(vr =vr, ref = ref_genome, k = 5, num_of_DELs = delet_this)

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
  if(add_ID_cats == TRUE){
    vcf <- add_ID_cats_ToVCF(vcf = vcf,palimpdir = palimpdir, ref_fasta = ref_fasta)
  }
  return(vcf)
}



#' palimpsest_addMutationContextToVR_2
#'
#' Function used within annotate_VCF to add mutation context to SBS mutations. 

#' @return vcf
#' @import gtools
#' @import Biostrings
#' @import GenomicRanges
#' @examples
#' vr3 <- palimpsest_addMutationContextToVR_2(vr = vr, ref_genome = ref_genome, k = 3, unify = TRUE,num_of_DELs = delet_this)




palimpsest_addMutationContextToVR_2 <- function (vr =NULL, ref_genome = NULL, k = 3,  check.strand = FALSE, 
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

