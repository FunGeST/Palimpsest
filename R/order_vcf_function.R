#' order_vcf
#'
#' Orders a VCF by position, chromosome, sample and project (if project argument given).
#' @param vcf The VCF to be ordered.
#' @param Project_col Name of the Project column in the VCF (e.g. may contain "ICGC")
#' @keywords Signatures
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
    tmp <- arrange(tmp,tmp$POS)
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
