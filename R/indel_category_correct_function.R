#' indel_category_correct
#'
#' Corrects Indel category output from PCAWG7-data-preparation-version-1.5
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
