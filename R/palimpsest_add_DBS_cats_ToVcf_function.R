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
  return(output)
}