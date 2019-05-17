compare_results <- function(reference_sigs = NA, extraction_1 = NA, extraction_2 = NULL, extraction_1_name ="Palimp", extraction_2_name = NA, lower_threshold = 0.6, upper_threshold = 0.9){

  #ifelse(!is.na(extraction_2),print("good to go"),stop("function currently needs 2 denovo sets to work"))
  if(!missing(extraction_2)){
    rownames(extraction_1) <- c(paste0(extraction_1_name,"_",c(1:nrow(extraction_1))))
    rownames(extraction_2) <- c(paste0(extraction_2_name,"_",c(1:nrow(extraction_2))))
    
    refs <- rownames(reference_sigs)
    denovs1 <- rownames(extraction_1)
    denovs2 <- rownames(extraction_2)
    
    longest <- max(c(length(denovs1),length(denovs2)))
    
    res <- as.data.frame(matrix(nrow = nrow(reference_sigs)+longest,ncol = 0)) %>% 
      mutate(References = NA, denovo1 = NA, ref_denovo1_cos = NA, denovo2 = NA, ref_denovo2_cos = NA, denovo1_denovo2_cos = NA, keep_version = NA)
  
    res$References[1:nrow(reference_sigs)] <- c(rownames(reference_sigs))
                                                
                                                
    # perfrom cosine similarity analysis
    mutmat1 <- t(rbind(extraction_1, reference_sigs))
    m1 <- as.matrix(palimpsest_distCosine(t(mutmat1)))
    m1 <- 1-m1
    mutmat2 <- t(rbind(extraction_2, reference_sigs))
    m2 <- as.matrix(palimpsest_distCosine(t(mutmat2)))
    m2 <- 1-m2
    mutmat3 <- t(rbind(extraction_1, extraction_2))
    m3 <- as.matrix(palimpsest_distCosine(t(mutmat3)))
    m3 <- 1-m3
    
    # make results table
    for(i in 1:nrow(reference_sigs)){
      sig <- res$References[i]
      if(max(m1[sig,colnames(m1) %in% denovs1]) > lower_threshold){
        res$ref_denovo1_cos[i] <- max(m1[sig,colnames(m1) %in% denovs1])
        res$denovo1[i] <- rownames(m1)[m1[sig,]==res$ref_denovo1_cos[i]]
      }
      if( max(m2[sig,colnames(m2) %in% denovs2]) > lower_threshold){
        res$ref_denovo2_cos[i] <- max(m2[sig,colnames(m2) %in% denovs2])
        res$denovo2[i] <- rownames(m2)[m2[sig,]==res$ref_denovo2_cos[i]]
      }
    }
    
    manque1 <- denovs1[denovs1 %!in% res$denovo1]
    manque2 <- denovs2[denovs2 %!in% res$denovo2]
    
    if(length(manque1) > 0){
      for(i in 1:length(manque1)){
          res$denovo1[length(refs)+i] <- manque1[i]
          sig <- manque1[i]
          res$denovo1_denovo2_cos[length(refs)+i] <- max(m3[sig,colnames(m3) %in% denovs2])
          res$denovo2[length(refs)+i] <- rownames(m3)[m3[sig,]==res$denovo1_denovo2_cos[length(refs)+i]]
      }
    }
    for(i in 1:(length(refs)+length(manque1))){
      if(!is.na(res$denovo1[i]) & !is.na(res$denovo2[i])) res$denovo1_denovo2_cos[i]<- m3[res$denovo1[i],res$denovo2[i]]
    }
    
    manque2x <- manque2[manque2 %!in% res$denovo2]
    if(length(manque1) > 0){
      for(i in 1:length(manque2x)){
        if(length(manque2x)==0) break
        res$denovo2[length(refs)+i+length(manque1)] <- manque2x[i]
      }
    }
    
    for(i in (length(refs)+1):nrow(res)){
      sig1 <- res$denovo1[i]; sig2 <- res$denovo2[i]
      if(!is.na(sig1)){
        if(max(m1[sig1,colnames(m1) %in% refs]) > lower_threshold){
          res$ref_denovo1_cos[i] <- max(m1[sig1,colnames(m1) %in% refs])
          res$References[i] <- colnames(m1)[m1[sig1,]==res$ref_denovo1_cos[i]]
          if(!is.na(sig2)) res$ref_denovo2_cos[i] <- m2[sig2,res$References[i]]
        }
      }
      if(!is.na(sig2)){
        if(max(m2[sig2,colnames(m2) %in% refs]) > lower_threshold){
          res$ref_denovo2_cos[i] <- max(m2[sig2,colnames(m2) %in% refs])
          res$References[i] <- colnames(m2)[m2[sig2,]==res$ref_denovo2_cos[i]]
          if(!is.na(sig1)) res$ref_denovo1_cos[i] <- m1[sig1,res$References[i]]
        }
      }
    }
    
    # add suggestion annotations
    res_f = res %>% 
      filter(rowSums(is.na(.)) < ncol(.)) %>% 
      mutate(keep = ifelse((ref_denovo1_cos > lower_threshold & ref_denovo2_cos > lower_threshold)|ref_denovo1_cos > upper_threshold | ref_denovo2_cos > upper_threshold,"Yes",NA)) #mark those to keep
    
    res_f$keep[res_f$denovo1_denovo2_cos > upper_threshold] <- "Yes"
    res_f$keep[is.na(res_f$keep)] <- "No"
    
    res_f = res_f %>%
      mutate(keep_version = ifelse(keep == "Yes" & References %in% refs & ref_denovo1_cos > upper_threshold | ref_denovo2_cos > upper_threshold , "Keep_Reference",keep_version)) #highlight for which ref is defo best
      
    res_f = res_f %>% 
      mutate(keep_version = ifelse(keep == "Yes" & is.na(keep_version),"Check_denovo_results",keep_version)) #highlight for which user should check
    res_f$keep_version[res_f$keep == "No"] <- "Discard"
    duped <- unique(res_f$References[duplicated(res_f$References) & !is.na(res_f$References)])
    res_f$keep_version[res_f$References %in% duped] <- "Check_denovo_results"
    
    colnames(res_f) <- c("Ref_Signature",paste0(extraction_1_name,"_Equivalent"),paste0("Ref_",extraction_1_name,"_cosine_score"),paste0(extraction_2_name,"_Equivalent"),
                         paste0("Ref_",extraction_2_name,"_cosine_score"), paste0(extraction_1_name,"_",extraction_2_name,"_cosine_score"),"Suggested_Action","Keep")
    
    res_f <- res_f[,c(1:7)]
    res_f
  }
  else{
    rownames(extraction_1) <- c(paste0(extraction_1_name,"_",c(1:nrow(extraction_1))))

  refs <- rownames(reference_sigs)
  denovs1 <- rownames(extraction_1)

  
  res <- as.data.frame(matrix(nrow = nrow(reference_sigs)+length(denovs1),ncol = 0)) %>% 
    mutate(References = NA, denovo1 = NA, ref_denovo1_cos = NA, keep_version = NA)
  
  res$References[1:nrow(reference_sigs)] <- c(rownames(reference_sigs))
  
  # perfrom cosine similarity analysis
  mutmat1 <- t(rbind(extraction_1, reference_sigs))
  m1 <- as.matrix(palimpsest_distCosine(t(mutmat1)))
  m1 <- 1-m1

  
  # make results table
  for(i in 1:nrow(reference_sigs)){
    sig <- res$References[i]
    if(max(m1[sig,colnames(m1) %in% denovs1]) > lower_threshold){
      res$ref_denovo1_cos[i] <- max(m1[sig,colnames(m1) %in% denovs1])
      res$denovo1[i] <- rownames(m1)[m1[sig,]==res$ref_denovo1_cos[i]]
    }
  }
  
  manque1 <- denovs1[denovs1 %!in% res$denovo1]
  
  if(length(manque1) > 0){
    for(i in 1:length(manque1)){
      res$denovo1[length(refs)+i] <- manque1[i]
    }
  }

  
  
  # add suggestion annotations
  res_f = res %>% 
    filter(rowSums(is.na(.)) < ncol(.)) %>% 
    mutate(keep = ifelse((ref_denovo1_cos > lower_threshold)|ref_denovo1_cos > upper_threshold,"Yes",NA)) #mark those to keep
  
  res_f$keep[is.na(res_f$keep)] <- "No"
  
  res_f = res_f %>%
    mutate(keep_version = ifelse(keep == "Yes" & References %in% refs & ref_denovo1_cos > upper_threshold, "Keep_Reference",keep_version)) #highlight for which ref is defo best
  
  res_f = res_f %>% 
    mutate(keep_version = ifelse(keep == "Yes" & is.na(keep_version),"Check_denovo_results",keep_version)) #highlight for which user should check
  res_f$keep_version[res_f$keep == "No"] <- "Discard"
  duped <- unique(res_f$References[duplicated(res_f$References) & !is.na(res_f$References)])
  res_f$keep_version[res_f$References %in% duped] <- "Check_denovo_results"
  
  colnames(res_f) <- c("Ref_Signature",paste0(extraction_1_name,"_Equivalent"),paste0("Ref_",extraction_1_name,"_cosine_score"),
                       "Suggested_Action","Keep")
  
  res_f <- res_f[,c(1:4)]
  res_f
    
}
}




