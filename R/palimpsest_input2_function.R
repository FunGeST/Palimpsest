#' Palimpsest_input2
#'
#' Creates the input matrices of the numbers and proportions of each mutation category (of one mutation type, e.g. SBS mutation categories) in each sample.
#' @param vcf The input VCF file from which a matrix of mutation categories is to be made.
#' @param mutypes leave blank to use default COSMIC categories (recommeneded), or add your own. 
#' @param Type Mutation type (SBS, DBS, ID or SV).
#' @keywords Signatures
#' @export
#' @examples
#' DBS_matrix <- Palimpsest_input2(vcf=vcf, Type = "DBS", proportion = TRUE)



Palimpsest_input2 <- function(vcf = NULL, Type = NULL, mutypes = NA){
  
  if(Type == "SBS") mutcat.col <- "SBS_cat3"; if(Type == "DBS") mutcat.col <- "DBS_cat"; if(Type == "ID") mutcat.col <- "ID_cat"
  if(mutcat.col %!in% colnames(vcf)) stop(paste("vcf is missing the",Type,"mutation category column. Use the << annotate_VCF >> function to add the appropriate column."))
  
  ordering <- unique(vcf$Sample)
  if(is.na(mutypes)){
     if(Type == "SBS"){
      bases <- c("A", "C", "G", "T")
      ctxt16 <- paste(rep(bases, each = 4), rep(bases, 4), 
                      sep = ".")
      mt <- c("CA", "CG", "CT", "TA", "TC", "TG")
      Types96 <- paste(rep(mt, each = 16), rep(ctxt16, 6), 
                       sep = "_")
      mutypes <- paste0(substr(Types96,1,4), substr(Types96,1,1), substr(Types96,6,6))
      
    }
    if(Type == "DBS"){
      mutypes <- c("AC>CA", "AC>CG", "AC>CT", "AC>GA", "AC>GG", "AC>GT", "AC>TA", "AC>TG", "AC>TT", "AT>CA", "AT>CC", "AT>CG", "AT>GA", "AT>GC", "AT>TA", "CC>AA", "CC>AG",
                   "CC>AT", "CC>GA" ,"CC>GG" ,"CC>GT", "CC>TA", "CC>TG", "CC>TT", "CG>AA", "CG>AC", "CG>AT", "CG>GA" ,"CG>GC" ,"CG>TA", "CT>AA", "CT>AC", "CT>AG", "CT>GA",
                   "CT>GC" ,"CT>GG" ,"CT>TA" ,"CT>TC", "CT>TG", "GC>AA", "GC>AG", "GC>AT" ,"GC>CA" ,"GC>CG" ,"GC>TA" ,"TA>AC", "TA>AG", "TA>AT", "TA>CC", "TA>CG", "TA>GC",
                   "TC>AA", "TC>AG", "TC>AT", "TC>CA", "TC>CG", "TC>CT", "TC>GA", "TC>GG", "TC>GT", "TG>AA", "TG>AC", "TG>AT", "TG>CA", "TG>CC", "TG>CT", "TG>GA", "TG>GC",
                   "TG>GT", "TT>AA", "TT>AC", "TT>AG", "TT>CA", "TT>CC", "TT>CG", "TT>GA", "TT>GC", "TT>GG")
    }
    if(Type == "ID"){
      mutypes = c("DEL_C_1_0","DEL_C_1_1","DEL_C_1_2","DEL_C_1_3","DEL_C_1_4","DEL_C_1_5+", "DEL_T_1_0","DEL_T_1_1","DEL_T_1_2",
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
        mutypes <- gsub("_",",",mutypes)
      
    }
  }
  if(Type == "SV"){
    tt <- as.matrix(table(vcf[, mutcat.col], vcf$Sample))
    mat <- matrix(as.numeric(tt), nrow = nrow(tt), ncol = ncol(tt))
    rownames(mat) <- rownames(tt)
    colnames(mat) <- colnames(tt)
    data <- mat
    categs <- c("DEL_0-1kb_clust_1", "DEL_1-10kb_clust_1", 
                "DEL_10-100kb_clust_1", "DEL_100kb-1Mb_clust_1", 
                "DEL_1-10Mb_clust_1", "DEL_>10Mb_clust_1", "DUP_0-1kb_clust_1", 
                "DUP_1-10kb_clust_1", "DUP_10-100kb_clust_1", "DUP_100kb-1Mb_clust_1", 
                "DUP_1-10Mb_clust_1", "DUP_>10Mb_clust_1", "INV_0-1kb_clust_1", 
                "INV_1-10kb_clust_1", "INV_10-100kb_clust_1", "INV_100kb-1Mb_clust_1", 
                "INV_1-10Mb_clust_1", "INV_>10Mb_clust_1", "BND_clust_1", 
                "DEL_0-1kb_clust_0", "DEL_1-10kb_clust_0", "DEL_10-100kb_clust_0", 
                "DEL_100kb-1Mb_clust_0", "DEL_1-10Mb_clust_0", "DEL_>10Mb_clust_0", 
                "DUP_0-1kb_clust_0", "DUP_1-10kb_clust_0", "DUP_10-100kb_clust_0", 
                "DUP_100kb-1Mb_clust_0", "DUP_1-10Mb_clust_0", "DUP_>10Mb_clust_0", 
                "INV_0-1kb_clust_0", "INV_1-10kb_clust_0", "INV_10-100kb_clust_0", 
                "INV_100kb-1Mb_clust_0", "INV_1-10Mb_clust_0", "INV_>10Mb_clust_0", 
                "BND_clust_0")
    missings_cat = setdiff(categs, rownames(data))
    for (i in missings_cat) {
      data = rbind(data, assign(i, rep(0, ncol(data))))
      rownames(data)[nrow(data)] = i
    }
    if (proportion) {
      data <- prop.table(data, 2)
    }
    return(as.data.frame(data))
  }else{
    tmp <- split(vcf, vcf$Sample)
    truenames <- names(tmp)
    tmp <- lapply(tmp, function(d) {
      sapply(mutypes, function(m) {
        sum(d[, mutcat.col] == m, na.rm = T)
      })
    })
    names(tmp) <- truenames; tmp <- tmp[ordering]
    nums <- as.matrix(as.data.frame(tmp))
    
    if(Type == "DBS"){
      trip_nums <- count_triplets(vcf, ordering, mutypes)
      
      nums <- nums - trip_nums
      nums <- nums/2
      nums <- nums+trip_nums
    }
    
    nums <- nums[,(colSums(nums)>0)]
    props <- nums
    for (j in 1:ncol(props)) props[, j] <- props[, j]/sum(props[,j])
    
    res = list(mut_nums = nums, mut_props = props)
    return(res)
  }
}








#' count_triplets
#'
#' Counts the number of DBS mutations that are in fact triple base substituions. Used inside "palimpsest_input2", not designed to be used outside this function.
#' @param vcf The input vcf file.
#' @param ordering the ordering of samples from palimpsest_input2 function.
#' @param mutypes The DBS mutypes from palimpsest_input2 function.
#' @keywords Signatures
#' @export
#' @examples
#' trip_nums <- count_triplets(vcf, ordering, mutypes)



count_triplets <- function(vcf, ordering, mutypes){
  
  res <- as.data.frame(matrix(ncol = 4,nrow = nrow(vcf)))
  lizt <- c()
  for(i in 1:nrow(vcf)){
    lizt[i] <- vcf$POS[i]
    if(i>2){  
      if(lizt[i]==(lizt[i-1])+1){
        if(lizt[i]==(lizt[i-2])+2){
          res[i,] <- vcf[i,c("Sample","CHROM","POS","DBS_cat")]
          res[i-1,] <- vcf[i-1,c("Sample","CHROM","POS","DBS_cat")]
          res[i-2,] <- vcf[i-2,c("Sample","CHROM","POS","DBS_cat")]
        }
      }
    }
  }
  colnames(res) <- c("Sample","CHROM","POS","DBS_cat")
  res <- res[-which(is.na(res[,1])),]
  newres <- as.data.frame(matrix(ncol = 4,nrow = nrow(res)));   colnames(newres) <- c("Sample","CHROM","POS","DBS_cat")
  for(i in 1:(nrow(res)/3)){
    tmp <- res[((3*i)-2):(3*i),]
    if(length(unique(tmp$CHROM))==1) newres[((3*i)-2):(3*i),] <- tmp
  }
  res <- newres[-which(is.na(newres[,1])),]
  res$rm <- NA
  for(i in 2:(nrow(res)-1)){
    if(res$DBS_cat[i] == res$DBS_cat[i-1] & res$Sample[i]==res$Sample[i-1] &  res$CHROM[i]==res$CHROM[i-1]){
      if(i==2){ 
        res$rm[i] <- "rm"
        res$rm[i-1] <- "rm"
      }else{
        if(res$DBS_cat[i] != res$DBS_cat[i-2] & res$POS[i] != (res$POS[i-2])+2){ 
          res$rm[i] <- "rm"
          res$rm[i-1] <- "rm"
        }
        if(res$DBS_cat[i] == res$DBS_cat[i-1] & 
           res$POS[i] == (res$POS[i-1])+1 & 
           res$POS[i] == (res$POS[i+1])-1 & 
           res$Sample[i] == res$Sample[i+1] & 
           res$Sample[i-1] == res$Sample[i]){ 
          res$rm[i] <- "rm"
          res$rm[i-1] <- "rm"
        }
      }
    }
  }
  
  
  res <- res[-which(res$rm == "rm"),] %>% 
    dplyr::select(Sample, DBS_cat)
  
  
  tmp <- split(res, res$Sample)
  tripnames <- names(tmp)
  tmp <- lapply(tmp, function(d) {
    sapply(mutypes, function(m) {
      sum(d[, "DBS_cat"] == m, na.rm = T)
    })
  })
  names(tmp) <- tripnames
  
  trip_nums <- as.data.frame(tmp)
  for(i in 1:length(ordering)){
    if(ordering[i] %!in% colnames(trip_nums)){
      new_col_num <- ncol(trip_nums)+1
      trip_nums[,new_col_num] <- 0; colnames(trip_nums)[new_col_num] <- ordering[i]
    }
  }
  trip_nums <- trip_nums[,ordering]
  return(trip_nums)
}