#' plot_signatures
#'
#' Function to plot signature profiles of all mutation types (SBS, DBS or Indel - there is a different function for plotting SV signatures). 
#' N.B. recommended output width to height ratio of 24:7 must be maintained to keep labels at the correct size.
#' @param input_data Matrix of signature(s) to plot, or an individual sample's proportion of each mutation category. If you receive an input_data related error retry with t(input_data).
#' @param Title List (must be same length as number of input signatures) or single title to use. Leave blank to use rownames from input matrix. Use " " for no title. 
#' @param label "Full" for all labels (default), "Top" for just coloured bars atop the plot; "None" for just bar graph.
#' @keywords Signatures
#' @export
#' @import scales
#' @import ggplot2
#' @examples
#' pdf(file=paste0(resdir,"SBS_liver_signature_profiles.pdf"),width=24,height=7)
#' plot_signatures(input_data = SBS_liver,Title = rownames(SBS_liver))
#' dev.off()

plot_signatures <- function (input_data = NULL, Title = NA, label = "Full") {
  requireNamespace("scales", quietly = TRUE); requireNamespace("ggplot2", quietly = TRUE)
  if(label %!in% c("Full","Top","Bottom","None")) stop("label must be 'Full', 'Top', 'Bottom' or 'None'")
  Individual <- ifelse(1 %in% dim(input_data),TRUE,FALSE)
  if (Individual == FALSE) num_of_sigs = nrow(input_data)
  if (Individual == TRUE) num_of_sigs = 1
  input_data <- as.matrix(input_data)
  
  for (i in 1:num_of_sigs) {
    
    if (Individual == FALSE){
      if(ncol(input_data) %!in% c(78,83,96)) stop("input_data format incorrect")
      if(ncol(input_data)==96) Type <- "SBS"; if(ncol(input_data)==78) Type <- "DBS"; if(ncol(input_data)==83) Type <- "ID" 
      Context <- colnames(input_data)
      plot_input <- as.data.frame(input_data[i,])
      if(length(Title)==1) plot_title <- paste(Title,i,sep=".")
      if(length(Title)>1) plot_title <- Title[i]
      if(is.na(Title[1])) plot_title <- ""
    }
    if (Individual == TRUE){
      if(nrow(input_data)==1) input_data <- t(input_data)
      if(nrow(input_data) %!in% c(38,78,83,96)) stop("input_data format incorrect")
      if(nrow(input_data)==96) Type <- "SBS"; if(nrow(input_data)==78) Type <- "DBS"; if(nrow(input_data)==83) Type <- "ID"; if(nrow(input_data)==38) Type <- "SV" 
      Context <- rownames(input_data)
      plot_input <- data.frame(input_data)
      if(!is.na(Title))plot_title <- Title
      if(is.na(Title)) plot_title <- ""
    }
    
    colnames(plot_input)[1] <- "freq"; max.y = max(plot_input$freq)
    
    if (Type == "DBS"){
      plot_input = plot_input %>% 
        mutate("Context" = Context,
               "Substype"= paste(substr(Context, 1, 2),">NN",sep=""),
               "Substype_short" = paste(substr(Context, 4, 5)),
               "Substype_blank" = NA,
               "Colours" = c(rep("lightblue",9),rep("darkblue",6),rep("darkolivegreen3",9),rep("green4",6),
                             rep("lightpink",9),rep("red",6),rep("sandybrown",6),rep("darkorange1",9),
                             rep("thistle",9),rep("darkviolet",9)),
               "Context_numerical" = c(1:nrow(plot_input)),
               "Substype_numerical" = c(as.numeric(as.factor(Substype))),
               "Start_pos" = c(rep(1,9),rep(10,6),rep(16,9),rep(25,6),rep(31,9),rep(40,6),rep(46,6),rep(52,9),rep(61,9),rep(70,9) ),
               "Substype_length" = c(rep(9,9),rep(6,6),rep(9,9),rep(6,6),rep(9,9),rep(6,6),rep(6,6),rep(9,9),rep(9,9),rep(9,9)) ) 

      
      prev <- "n'importe quoi"
      for(j in 1:nrow(plot_input)){
        if(plot_input$Substype[j] != prev){
          plot_input$Substype_blank[j] <- plot_input$Substype[j]
          prev <- plot_input$Substype[j]
          next
        }else{
          plot_input$Substype_blank[j] <- paste(" ")
        }
      }
      
      xmax_correction = 9
    }
    
    if (Type == "SBS"){
      plot_input = plot_input %>% 
        mutate("Substype"= paste0(substr(Context, 1, 1),">",substr(Context,2,2)),
               "Context" = paste0(substr(Context,4,4),substr(Context,1,1),substr(Context,6,6)),
               "Substype_blank" = NA,
               "Colours" = c(rep("skyblue3", 16), rep("black", 16), rep("red",16), rep("grey", 16), rep("green", 16), rep("pink", 16)),
               "Context_numerical" = c(1:nrow(plot_input)),
               "Start_pos" = c(rep(1,16), rep(17,16), rep(33,16), rep(49,16), rep(65,16), rep(81,16) ),
               "Substype_length" = 16,
               "Substype_numerical" = c(as.numeric(as.factor(Substype))) ) #%>% 

      
      prev <- "rien de rien"
      for(j in 1:nrow(plot_input)){
        if(plot_input$Substype[j] != prev){
          plot_input$Substype_blank[j] <- plot_input$Substype[j]
          prev <- plot_input$Substype[j]
          next
        }else{
          plot_input$Substype_blank[j] <- paste(" ")
        }
      }
      
      xmax_correction = 20
    }
    
    
    if (Type == "ID"){
      plot_input = plot_input %>% 
        mutate("Context" = as.factor(Context),
               "x_label" = c("1 ","2 ","3 ","4 ","5 ","6+","1 ","2 ","3 ","4 ","5 ","6+","0 ","1 ","2 ","3 ","4 ","5+","0 ","1 ","2 ","3 ","4 ","5+","1 ","2 ","3 ","4 ","5 ","6+","1 ","2" ,"3 ","4 ","5 ","6+","1 ","2 ","3 ","4 ","5 ","6+","1 ","2 ","3 ","4 ","5 ","6+","0 ","1 ","2 ","3 ","4 ","5+","0 ","1 ","2 ","3 ","4 ","5+","0 ","1 ","2 ","3 ","4 ","5 ","0 ","1 ","2 ","3 ","4 ","5+","1 ","1 ","2 ","1 ","2 ","3 ","1 ","2 ","3 ","4 ","5+"),
               "Substype" =  c(rep(1,each=6), rep(2,each=6), rep(3,each=6), rep(4,each=6), rep(5,each=6), 
                               rep(6,each=6), rep(7,each=6), rep(8,each=6), rep(9,each=6), rep(10,each=6), rep(11,each=6), 
                               rep(12,each=6), 13, rep(14,each=2), rep(15,each=3), rep(16,each=5)),
               "Colours" = c(rep("sandybrown",6),rep("darkorange1",6),rep("darkolivegreen3",6),rep("green4",6),
                             rep("pink",6),rep("salmon",6),rep("red2",6),rep("red4",6),rep("lightskyblue1",6),rep("lightskyblue3",6),
                             rep("steelblue3",6),rep("steelblue4",6),"mediumpurple1",
                             rep("purple2",2),rep("mediumpurple4",3),rep("purple4",5)),
               "Context_numerical" = c(1:nrow(plot_input)),
               "Start_pos" = c(rep(1,6),rep(7,6),rep(13,6),rep(19,6),rep(25,6),rep(31,6),rep(37,6),rep(43,6),rep(49,6),rep(55,6),rep(61,6),rep(67,6),73,rep(74,2),rep(76,3),rep(79,5)),#Context_numerical[1],
               "Substype_length" = c(rep(6,72),1,2,2,3,3,3,5,5,5,5,5) 
        ) 
      
      indel_lab <- as.data.frame(matrix(nrow=16,ncol = 0))
      indel_lab$Substype <- c(1:16)
      indel_lab$label_black <- c("C"," ","C"," ","2","3"," "," ","2","3"," "," "," "," "," "," ")
      indel_lab$label_white <- c(" ","T"," ","T"," "," ","4","5+"," "," ","4","5+"," "," "," ","5+")
      indel_lab$label_white2 <-c(" "," "," "," "," "," "," "," "," "," "," "," ","2"," "," "," ")
      indel_lab$label_white3 <-c(" "," "," "," "," "," "," "," "," "," "," "," "," ","3"," "," ")
      indel_lab$label_white4 <-c(" "," "," "," "," "," "," "," "," "," "," "," "," "," ","4"," ")
      plot_input <- left_join(plot_input, indel_lab, by = "Substype")
      
      xmax_correction = 5.5
    }
    
    max.y = max(plot_input$freq)
    
    ### PLot ----------------------------------------------------------------------------------------------------------------
    
    plot_output= ggplot(data=plot_input) +
      geom_bar(aes(x=as.factor(plot_input$Context_numerical),y=freq),stat="identity",fill=plot_input$Colours, width = 0.5)+
      ggtitle(plot_title) +
      ylab("") + xlab("") +
      coord_cartesian(xlim = c(1, nrow(plot_input)), ylim = c(0,max.y*1.3), clip = 'off') + #  # This keeps the labels from disappearing
      
      theme_bw() + theme(panel.border = element_blank(),
                         panel.grid.major.x = element_blank(),
                         panel.grid.major.y  = element_blank(),
                         panel.grid.minor.y = element_blank(),
                         axis.line = element_line(colour = "black"),
                         axis.text.x = element_blank(),
                         panel.background = element_blank(), plot.background = element_blank(),
                         axis.text.y = element_text(size = 15, colour = "black"),
                         axis.ticks.x = element_blank(),
                         plot.margin=unit(c(1.2,1.6,1.2,1.1),"cm"), # order = top, right, bottom, left (/never eat shredded wheat)
                         plot.title = element_text(size = 55, face = "bold", hjust = 0.02, vjust = -5)) 
    
    
    
    
    ### Add labels ----------------------------------------------------------------------------------------------------------------
    
    if (Type == "DBS" & label == "Full"){
      max.y_lab <- max.y; if((round(max.y_lab,2)*100)%%2==1) max.y_lab <- round(max.y_lab,2)+0.01
      plot_output = plot_output +
        
        xlab(" ") + ylab(" ") +
        scale_x_discrete("", labels = plot_input$Substype_short) +
        annotate("text",x = -2.7, y = max.y*0.5, label = "Proportion of Mutations", angle = 90, size = 7, fontface = "plain", colour = "black") +
        annotate("text",x = 39, y =   0-max.y*0.17, label = "Dinucleotide Base Alteration", size = 8, fontface = "plain", colour = "black") +
        theme(axis.text.x = element_text(angle = 90, vjust=0.5,  hjust = 0, size = 20, colour = "grey47", family = "mono"), 
              axis.text.y = element_text(size = 20),
              axis.ticks.y.right = element_blank()
        ) +
        
        geom_rect(mapping = aes(xmin = (plot_input$Start_pos - 0.5), ymin = max.y*1.25, xmax = (plot_input$Start_pos  + xmax_correction),
                                ymax= max.y*1.375), fill = plot_input$Colours, color = "white", size = 2.5) +
        geom_text(position = "identity", mapping = aes(x = plot_input$Start_pos + (plot_input$Substype_length/2.2), 
                                                       y = max.y*1.42, label = plot_input$Substype_blank), size = 8) +
        scale_y_continuous(expand = c(0, 0), limits = c(-1, max.y*1.42),breaks = seq(from=0,to=max.y_lab,by=max.y_lab/2),labels = scales::number_format(accuracy = 0.01),
                           sec.axis = sec_axis(~.+0,labels = NULL)) +
        geom_rect(mapping = aes(xmin = (78.5), ymin = max.y*1.25, xmax = 79,
                                ymax= max.y*1.375), fill = "white", color = "white", size = 2.5) #+
    } 
    
    if (Type == "SBS"){
      if(label == "Full"){
        max.y_lab <- max.y; if((round(max.y_lab,2)*100)%%2==1) max.y_lab <- round(max.y_lab,2)+0.01
        plot_output = plot_output +
          
          xlab(" ") + ylab(" ") +
          scale_x_discrete("", labels = plot_input$Context) +
          annotate("text",x = -3.5, y = max.y*0.5, label = "Proportion of Mutations", angle = 90, size = 7, fontface = "plain", colour = "black") +
          annotate("text",x = 48, y =   0-max.y*0.22, label = "Trinucleotide Context", size = 8, fontface = "plain", colour = "black") +
          theme(axis.text.x = element_text(angle = 90, vjust=0.5,  hjust = 0, size = 20, colour = "grey47", family = "mono"), 
                axis.text.y = element_text(size = 20),
                axis.ticks.y.right = element_blank()
          ) +
          
          geom_rect(mapping = aes(xmin = (plot_input$Start_pos - 0.5), ymin = max.y*1.25, xmax = (plot_input$Start_pos  + xmax_correction),
                                  ymax= max.y*1.375), fill = plot_input$Colours, color = "white", size = 2.5) +
          geom_text(position = "identity", mapping = aes(x = plot_input$Start_pos + (plot_input$Substype_length/2.2), 
                                                         y = max.y*1.42, label = plot_input$Substype_blank), size = 8) +
          scale_y_continuous(expand = c(0, 0), limits = c(-1, max.y*1.42),breaks = seq(from=0,to=max.y_lab,by=max.y_lab/2),labels = scales::number_format(accuracy = 0.01),
                             sec.axis = sec_axis(~.+0,labels = NULL)) +
          geom_rect(mapping = aes(xmin = (96.5), ymin = max.y*1.25, xmax = 100,
                                  ymax= max.y*1.375), fill = "white", color = "white", size = 2.5) #+
      }
    } 
    
    
    
    if(label == "Full" & Type == "ID"){
      max.y_lab <- max.y; if((round(max.y_lab,2)*100)%%2==1) max.y_lab <- round(max.y_lab,2)+0.01
      x <- 1.31
      plot_output = plot_output +
        
        scale_x_discrete("", labels = plot_input$x_label) +
        scale_y_continuous(expand = c(0, 0), limits = c(-1, max.y*1.52),breaks = seq(from=0,to=max.y_lab,by=max.y_lab/2),labels = scales::number_format(accuracy = 0.01),
                           sec.axis = sec_axis(~.+0,labels = NULL)) +  # This keeps the labels from disappearing
        theme(axis.text.x = element_text(angle = 0, vjust=0,  hjust = 0.33, size = 17.5, colour = as.character(plot_input$Colours)),
              axis.ticks.x = element_line(),
              axis.text.y = element_text(size = 20),
              axis.ticks.y.right = element_blank()
        ) +
        geom_rect(mapping = aes(xmin = (plot_input$Start_pos - 0.5), ymin = max.y*1.25, xmax = (plot_input$Start_pos  + xmax_correction),
                                ymax= max.y*1.375), fill = plot_input$Colours, color = "white", size = 2.5) +
        geom_rect(mapping = aes(xmin = (83.5), ymin = max.y*1.25, xmax = 85,
                                ymax= max.y*1.375), fill = "white", color = "white", size = 2.5) +
        
        
        geom_text(position = "identity", mapping = aes(x = plot_input$Start_pos + (plot_input$Substype_length/2.3), 
                                                       y = max.y*x, label = plot_input$label_black), size = 8) +
        geom_text(position = "identity", mapping = aes(x = plot_input$Start_pos + (plot_input$Substype_length/2.3), 
                                                       y = max.y*x, label = plot_input$label_white), size = 8, colour = "white") +
        geom_text(mapping = aes(x = plot_input$Start_pos , y = max.y*x, label = plot_input$label_white2),
                  size = 8, colour = "white", nudge_x = 0 ) +
        geom_text(mapping = aes(x = plot_input$Start_pos + (plot_input$Substype_length/4) , y = max.y*x, label = plot_input$label_white3),
                  size = 8, colour = "white", nudge_x = 0 ) +
        geom_text(mapping = aes(x = plot_input$Start_pos + (plot_input$Substype_length/3.2) , y = max.y*x, label = plot_input$label_white4),
                  size = 8, colour = "white", nudge_x = 0 ) +
        
        
        annotate("text",x = -2.7, y = max.y*0.5, label = "Proportion of Mutations", angle = 90, size = 7, fontface = "plain", colour = "black") +
        
        annotate("text",x = 6, y = max.y*1.42, label = "1bp Deletion", size = 8, fontface = "plain", colour = "black") +
        annotate("text", x = 18, y = max.y*1.42, label = "1bp Insertion", size = 8, fontface = "plain", colour = "black")  +
        annotate("text", x = 36.75, y = max.y*1.52, label = ">1bp deletions at repeats", size = 8, fontface = "plain", colour = "black") +
        annotate("text", x = 36.75, y = max.y*1.42, label = "(Deletion length)", size = 8, fontface = "plain", colour = "black") +
        annotate("text", x = 60, y = max.y*1.52, label = ">1bp insertion at repeat lengths", size = 8, fontface = "plain", colour = "black") +
        annotate("text", x = 60, y = max.y*1.42, label = "(Insertion length)", size = 8, fontface = "plain", colour = "black") +
        annotate("text", x = 78, y = max.y*1.52, label = "Deletions with microhomology", size = 8, fontface = "plain", colour = "black") +
        annotate("text", x = 78, y = max.y*1.42, label = "(Deletion length)", size = 8, fontface = "plain", colour = "black") +
        
        
        annotate("text", x = 7, y = 0-max.y*0.12, label = "Homopolymer length", size = 7) +
        annotate("text", x = 19, y = 0-max.y*0.12, label = "Homopolymer length", size = 7) +
        annotate("text", x = 37, y = 0-max.y*0.12, label = "Number of repeat units", size = 7) +
        annotate("text", x = 60, y = 0-max.y*0.12, label = "Number of repeat units", size = 7) +
        annotate("text", x = 77.5, y = 0-max.y*0.12, label = "Microhomology length", size = 7)
    }
    

    if(label=="None"){
      plot_output = plot_output +
        scale_y_continuous(expand = c(0, 0), limits = c(-1, max.y*1.42),
                           sec.axis = sec_axis(~.+0,labels = NULL)) +
        geom_hline(yintercept = max.y*1.3) +
        theme(axis.text.y = element_blank(),
              axis.ticks = element_blank())
    }
    if(label == "Top"){
      plot_output = plot_output +
        scale_y_continuous(expand = c(0, 0), limits = c(-1, max.y*1.5),
                           sec.axis = sec_axis(~.+0,labels = NULL)) +
        geom_hline(yintercept = max.y*1.3) +
        theme(axis.text.y = element_blank(),
              axis.ticks = element_blank()) +
        geom_rect(mapping = aes(xmin = (plot_input$Start_pos - 0.5), ymin = max.y*1.325, xmax = (plot_input$Start_pos  + xmax_correction),
                                ymax= max.y*1.45), fill = plot_input$Colours, color = "white", size = 2.5) +
        geom_rect(mapping = aes(xmin = nrow(plot_input)+0.5, ymin = max.y*1.325, xmax = nrow(plot_input)+5,
                                ymax= max.y*1.45), fill = "white", color = "white", size = 2.5)
    }
    if(label == "Bottom"){
      plot_output = plot_output +
        scale_y_continuous(expand = c(0, 0), limits = c(0, max.y*1.5),
                           sec.axis = sec_axis(~.+0,labels = NULL)) +
        geom_hline(yintercept = max.y*1.3) +
        theme(axis.text.y = element_blank(),
              axis.ticks = element_blank()) +
        geom_rect(mapping = aes(xmin = (plot_input$Start_pos - 0.5), ymin = max.y*-0.2, xmax = (plot_input$Start_pos  + xmax_correction),
                                ymax= 0), fill = plot_input$Colours, color = "white", size = 2.5) 
    }
    plot(plot_output)
  }
} 


#' deconvolution_exposure
#'
#' Plots the exposure of signatures across the series (bargraphs of the number and proportion of each signature in each sample).
#' @param signature_contribution List of signatures exposure numbers and proportions matrices (output from deconvolution_fit function).
#' @param rm_samples List of names of samples to remove (if any). Once the signature contribution of hypermutated samples is known, removing these from the grpahical output of this function can make the comparison of the rest of the samples easier. 
#' @param signature_colours Character vector of R-compatable colours representing each signature to be used graphical outputs. Each signature in input_signatures must have named colour in this vector for grpahical outputs to work. 
#' @keywords Signatures
#' @export
#' @import ggpubr
#' @import reshape2
#' @examples
#' pdf(file.path(resdir.,"signature_content_plot.pdf"),width=12,height=7)
#' signature_content_plot <- deconvolution_exposure(signature_contribution = SBS_signatures_exp, rm_samples = "CHC892T", signature_colours = sig_cols)
#' dev.off()

deconvolution_exposure <- function(signature_contribution = signatures_exp, rm_samples = NA,signature_colours = NA){
  
  requireNamespace("ggpubr", quietly = TRUE); requireNamespace("reshape2", quietly = TRUE)
  sig_nums <- signature_contribution$sig_nums; sig_props <- signature_contribution$sig_props
  if(!is.na(rm_samples)){
    sig_nums <- sig_nums[which(rownames(sig_nums) %!in% rm_samples),]
    sig_props <- sig_props[which(rownames(sig_props) %!in% rm_samples),]
  }
  
  nsamp <- nrow(sig_nums)
  nrows <- ceiling(nsamp/200)
  signature_content_plot <- c()
  if(nrows == 1){
    content_plot <- signature_exposure_plot(sig_nums, sig_props, signature_colours = signature_colours)
  }
  if(nrows>1){  ## for >200 samples we find it best to sort the plot into several rows to make it easier to view, the rest of the function does this
    
    nsamp_new <- nsamp
    for(j in 1:20){    ## making number of samples divisible by the number of rows 
      if( (nsamp_new/nrows)%%1 == 0){
        break
      }
      if( (nsamp_new/nrows)%%1 != 0){
        nsamp_new <- nsamp_new+1
      }
    }
    if(nsamp != nsamp_new){    ## adding "blank" samples to make number of samples up to number calculated in the previous step
      for(i in (1+nsamp):nsamp_new){
        sig_nums[i,] <- 0
        sig_props[i,] <- 0
        rownames(sig_nums)[i] <- paste(i)
        rownames(sig_props)[i] <- paste(i)
      }
    }
    
    for(i in 1:nrows){    ## sort samples into groups, one for each row
      n4me <- paste0("row_",i)
    
      if(i == 1) n4mes <- names(sort(rowSums(sig_nums)))[1:(nsamp_new/(nrows/i))]
      if(i > 1) n4mes <- names(sort(rowSums(sig_nums)))[((nsamp_new/(nrows/(i-1)))+1):(nsamp_new*i/nrows)]
      sig_nums_tmp <- sig_nums[n4mes,];sig_props_tmp <- sig_props[n4mes,]
      
      signature_content_plot[[i]] <- signature_exposure_plot(sig_nums_tmp,sig_props_tmp,signature_colours = signature_colours)

    }
    content_plot <-  ggarrange(plotlist =  signature_content_plot[nrows:1], heights = rep(7,nrows), widths = rep(10,nrows),
              ncol = 1, nrow = nrows)
  }
  print(content_plot)
}


#' signature_exposure_plot
#'
#' Function to plot the exposure of mutational signatures in samples across the entire series
#' @param mutSign_nums Matrix in sample x mutational signature exposure format in numbers
#' @param mutSign_props Matrix in sample x mutational signature exposure format in proportions
#' @param signature_colours Character vector indicating the colors representing each signature in graphical outputs. Must match to the total number of provided signatures
#'
#' @export
#' @import ggplot2
#' @importFrom ggplot2 theme_bw
#' @import gplots
#' @import reshape2
#' @importFrom reshape2 melt


signature_exposure_plot <- function (mutSign_nums, mutSign_props, signature_colours = NULL) 
{
  scale <- 1
  .theme_ss <- theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 90, 
                                                                           vjust = 0.5, size = 8 * scale, family = "mono"), axis.text.y = element_text(hjust = 0.5, 
                                                                                                                                                       size = 12 * scale, family = "mono"), axis.text = element_text(size = 12 * 
                                                                                                                                                                                                                       scale, family = "mono"))
  ordering <- order(colSums(t(mutSign_nums)), decreasing = T)
  mutSign_nums <- t(mutSign_nums)
  mutSign_props <- t(mutSign_props)
  mutSign_nums <- mutSign_nums[, ordering]
  mutSign_props <- mutSign_props[, ordering]
  sample.ordering <- colnames(mutSign_nums)
  x1 <- melt(mutSign_nums)
  x2 <- melt(mutSign_props)
  colnames(x1) <- c("Signature", "Sample", "Activity")
  colnames(x2) <- c("Signature", "Sample", "Activity")
  x1[, "class0"] <- c("Counts")
  x2[, "class0"] <- c("Proportions")
  df2 <- rbind(x1, x2)
  df2$class0 <- factor(df2$class0, c("Counts", "Proportions"))
  df2$Sample <- factor(df2$Sample, sample.ordering)
  p = ggplot(df2, aes(x = factor(Sample), y = Activity, fill = Signature))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + scale_fill_manual(values = signature_colours)
  p = p + ggtitle("Mutational Signature Exposures")
  p = p + facet_grid(class0 ~ ., scale = "free_y")
  p = p + theme(plot.title = element_text(lineheight = 1, face = "bold", 
                                          size = 15 * scale))
  p = p + xlab("Samples") + ylab("Mutational Signature Content")
  p = p + theme(axis.title.x = element_text(face = "bold", 
                                            colour = "black", size = 15 * scale))
  p = p + theme(axis.title.y = element_text(face = "bold", 
                                            colour = "black", size = 15 * scale))
  p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                           size = 12 * scale, face = "bold", colour = "black"))
  p = p + theme(axis.text.y = element_text(size = 10 * scale, 
                                           face = "bold", colour = "black"))
  p = p + theme(legend.title = element_blank())
  p = p + .theme_ss
  p = p + theme(legend.position = "top")
  return(p)
}
