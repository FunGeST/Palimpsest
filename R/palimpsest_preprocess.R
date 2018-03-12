#' preprocessInput_snv
#'
#' Annotating the mutation data with necessary fields for further analysis
#' @param input_data Mutation data in vcf format.
#' @param ensgene Gene table for annotations. A table of Ensembl genes is provided with the package.
#' @param reference_genome Reference genome (e.g. BSgenome.Hsapiens.UCSC.hg19)
#'
#' @return
#' @export
#' @import VariantAnnotation
#' @importFrom gtools mixedsort
#' @examples
preprocessInput_snv <- function(input_data = NULL, ensgene = ensgene, reference_genome=ref_genome)
{
  Sample.col = "Sample";CHROM.col = "CHROM";POS.col = "POS";REF.col = "REF";ALT.col = "ALT"
  chroms <- unique(input_data[, CHROM.col])
  if (1 %in% chroms == TRUE)
    input_data[, CHROM.col] <- paste("chr", input_data[, CHROM.col], sep = "")
  "%ni%" <- Negate("%in%")
  remove_chrs <- c("chrM", "chrMT")
  input_data <- input_data[which(input_data[, CHROM.col] %ni% remove_chrs), ]
  input_data$strand.mut <- "+"
  input_data$strand.gene <- NA
  print("Annotating mutation data:")
  vcf <- palimpsest_dfPosXSegm(input_data, dfPos.chrom.col = CHROM.col,
                               dfPos.pos.col = POS.col, ensgene, dfSegm.chrom.col = "Chromosome.Name",
                               dfSegm.start.col = "Gene.Start..bp.", dfSegm.end.col = "Gene.End..bp.",
                               colsToAdd = "Associated.Gene.Name", namesColsToAdd = "Associated.Gene.Name")

  vcf$strand.gene <- c("-", NA, "+")[ensgene[match(vcf$Associated.Gene.Name, ensgene$Associated.Gene.Name), "Strand"] + 2]
  print("Adding mutation categories:")
  vcf.snv <- palimpsest_addMutationContextToVcf(vcf[which(vcf$Type=="SNV"),], reference_genome,
                                            chrom.col = CHROM.col, start.col = POS.col, end.col = POS.col,
                                            ref.col = REF.col, alt.col = ALT.col, strand.mut.col = "strand.mut",
                                            strand.gene.col = "strand.gene", sample.col = Sample.col)
  vcf <- merge(vcf,vcf.snv,by=setdiff(intersect(names(vcf),names(vcf.snv)),"strand.mut"),all=TRUE,sort=FALSE)
  vcf$strand.mut <- vcf$strand.mut.x ; ind <- which(vcf$Type=="SNV") ; vcf[ind,"strand.mut"] <- vcf[ind,"strand.mut.y"] ; vcf <- vcf[,setdiff(names(vcf),c("strand.mut.x","strand.mut.y"))]
  return(vcf)
}



#' preprocessInput_sv
#'
#' Annotating the mutation data with necessary fields for further analysis
#' @param input_data Table describing somatic structural rearrangements.
#' @param ensgene Gene table for annotations. A table of Ensembl genes is provided with the package.
#' @param resdir Results directory where graphical outputs should be exported.
#'
#' @return
#' @export
#'
#' @examples
preprocessInput_sv <- function(input_data = NULL,ensgene = ensgene, resdir = resdir)
{
  Sample.col = "Sample"; CHROM_1.col = "CHROM_1"; CHROM_2.col = "CHROM_2"; POS_1.col = "POS_1"; POS_2.col = "POS_2"; type.col = "Type";
  chroms <- unique(input_data[, CHROM_1.col])
  if (1 %in% chroms == TRUE) {
    input_data[, CHROM_1.col] <- paste("chr", input_data[,CHROM_1.col], sep = "")
    input_data[, CHROM_2.col] <- paste("chr", input_data[,CHROM_2.col], sep = "")
  }
  input_data <- input_data[mixedorder(input_data[, CHROM_1.col]),
                           ]
  input_data$strand.mut <- "+"
  input_data$strand.gene <- NA
  print("Annotating mutation data:")
  vcf <- palimpsest_dfPosXSegm(input_data, dfPos.chrom.col = CHROM_1.col,
                               dfPos.pos.col = POS_1.col, ensgene, dfSegm.chrom.col = "Chromosome.Name",
                               dfSegm.start.col = "Gene.Start..bp.", dfSegm.end.col = "Gene.End..bp.",
                               colsToAdd = "Associated.Gene.Name", namesColsToAdd = "Associated.Gene.Name")
  print("Adding mutation categories:")
  vcf$strand.gene <- c("-", NA, "+")[ensgene[match(vcf$Associated.Gene.Name,
                                                   ensgene$Associated.Gene.Name), "Strand"] + 2]
  vcf <- palimpsest_addSVcategoriesToVcf(vcf, type.col = type.col,
                                             sample.col = Sample.col, CHROM_1.col = CHROM_1.col, CHROM_2.col = CHROM_2.col,
                                             POS_1.col = POS_1.col, POS_2.col = POS_2.col, resdir = resdir)
  return(vcf)
}
