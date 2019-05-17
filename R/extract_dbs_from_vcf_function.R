#' extract_dbs_from_vcf
#'
#' Extracts lines corresponding to DBS mutations from a VCF. N.B. The VCF must be ordered by sample, CHROM and position for this function to work.
#' @param vcf The nput VCF from which DBS mutations are to be extracted
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
