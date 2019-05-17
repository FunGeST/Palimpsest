# clustering function


clustering <- function(input_matrix = NULL, distHC = "distCosine", method = "Number", graph1_input = NULL, graph1_muttype = NULL,graph1_type = NULL,
                       graph2_input = NULL, graph2_type = NULL,  graph2_muttype = NULL,graph3_input = NULL, graph3_type=NULL,graph3_muttype = NULL, sig_cols=sig_cols, Waffle = TRUE, waf_matrix= waf_matrix){
  
  require(Palimpsest)
  
  if(distHC=="distCosine")  d <- distCosine(as.matrix(input_matrix))
  if(distHC %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski"))d <- dist((input_matrix),distHC)
  methHC <- c("ward","ward.D","ward.D2","single","complete","average")[2]
  
  hc_nums <- hclust(d,method=methHC)
  
  ordering <- hc_nums$labels [hc_nums$order]; print(ordering)
  nsamp <- length(ordering)
  
  
  # make HC tree
  dhc <- as.dendrogram(hc_nums)
  maxy0 <-  round(max(get_nodes_attr(dhc, "height")))*1.1
  
  ddata <- dendro_data(dhc, type = "rectangle") ; colnames(ddata$segments)[2] <- c("Height")
  arthur_dend <- ggplot(segment(ddata)) + 
    ggtitle(label = paste("Hierachical Clustering by",distHC)) +
    geom_segment(aes(x = x, y = Height, xend = xend, yend = yend)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(-0.1, maxy0 ),breaks = seq(0,maxy0,10)) +
    # scale_y_continuous(expand = c(0, 0), limits = c(-4.99, 37)) +
    scale_x_discrete(expand = c(nsamp/97800, 0)) +
    # geom_text(mapping = aes(x = x, y = y, label = label, angle = 90, hjust = 1.1, vjust = .5), 
    #           data = label(ddata), size = 3.75) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey88", size = 0.5),
          axis.line.x = element_blank(),
          axis.title.y = element_text(),
          plot.title = element_text(size = 25, vjust = 0)
          ,
          plot.margin=unit(c(1, (nsamp/575)*2, 0, nsamp/950), "cm") # order = top, right, bottom, left
          
    ) 
  
 

  
  
  # Plot 1 
  graph1_input <- graph1_input[ordering,]
  num_per_sample <- rowSums(graph1_input)[ordering]
  
  graph1_melt <- as.data.frame(melt(t(graph1_input))) ; colnames(graph1_melt) <- c("Signature", "Sample", "Number")
  graph1_melt$Total <- rep(as.numeric(num_per_sample),each=ncol(graph1_input))
  graph1_melt$Proportion <- graph1_melt$Number / graph1_melt$Total
  if(graph1_type == "Number"){
    if(graph1_muttype =="SBS"){
      toobig <- rownames(graph1_input)[rowSums(graph1_input)> 100000]
      hereswhy <- rowSums(graph1_input)[rowSums(graph1_input)> 100000]
      graph1_melt$Total[graph1_melt$Sample%in%toobig] <- 100000
      graph1_melt$Number[graph1_melt$Sample%in%toobig] <- graph1_melt$Proportion[graph1_melt$Sample%in%toobig] *100000
      x_adj <- 875
      
    }
    maxy1 <- max(graph1_melt$Total)*1.1
    if(graph1_muttype =="DBS")  x_adj <- 537.3626
    if(graph1_muttype =="ID")  x_adj <- 537.3626
  }
  if(graph1_type == "Proportion"){
    maxy1 <- 1
    x_adj <- 537.3626
  }
  
  plot1 <- ggplot(graph1_melt, aes(x = Sample, y = graph1_melt[,graph1_type])) +
    geom_bar(stat = "identity", position = "stack", aes(fill = (Signature)),width = 1) +
    ylab(graph1_type) +
    scale_fill_manual(values = sig_cols)  +
    #facet_grid(Type_f ~ ., scale = "free_y") +
    scale_y_continuous(expand = c(0, 0), limits = c(0,maxy1)) +
    theme_bw() +
    theme(
      legend.position = "none",
      #legend.title = element_blank(),
      # plot.margin=unit(c(0, (nsamp/489)*2, 0, nsamp/489), "cm"), # order = top, right, bottom, left
      plot.margin=unit(c(0, (nsamp/489)*2, 0, nsamp/x_adj), "cm"), # order = top, right, bottom, left
      
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 15),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      strip.text.y = element_text(size = 15)
    )
  
  
  
  
  
  # plot 2
  graph2_input <- graph2_input[ordering,]
  num_per_sample2 <- rowSums(graph2_input)[ordering]
  
  graph2_melt <- as.data.frame(melt(t(graph2_input))) ; colnames(graph2_melt) <- c("Signature", "Sample", "Number")
  graph2_melt$Total <- rep(as.numeric(num_per_sample2),each=ncol(graph2_input))
  graph2_melt$Proportion <- graph2_melt$Number / graph2_melt$Total
  if(graph2_type == "Number"){
    if(graph2_muttype =="SBS"){
      toobig <- rownames(graph1_input)[rowSums(graph2_input)> 100000]
      hereswhy <- rowSums(graph2_input)[rowSums(graph2_input)> 100000]
      graph2_melt$Total[graph2_melt$Sample%in%toobig] <- 100000
      graph2_melt$Number[graph2_melt$Sample%in%toobig] <- graph2_melt$Proportion[graph2_melt$Sample%in%toobig] *100000
      x_adj <- 875
    }
    maxy2 <- max(graph2_melt$Total)*1.1
    if(graph2_muttype =="DBS")  x_adj <- 537.3626
    if(graph1_muttype =="ID")  x_adj <- 537.3626
  }
  if(graph2_type == "Proportion"){
    maxy2 <- 1
    x_adj <- 537.3626
  }
  
  plot2 <- ggplot(graph2_melt, aes(x = Sample, y = graph2_melt[,graph2_type])) +
    geom_bar(width = 1,stat = "identity", position = "stack", aes(fill = (Signature))) +
    ylab(graph2_type) +
    scale_fill_manual(values = sig_cols)  +
    #facet_grid(Type_f ~ ., scale = "free_y") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 15, hjust = 0),
          legend.spacing.x = unit(0.1, 'cm'),
          plot.margin=unit(c(0, (nsamp/489)*2, 0, nsamp/x_adj), "cm"), # order = top, right, bottom, left
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0.1, "lines"),
          strip.text.y = element_text(size = 15)
    ) +
    guides(fill = guide_legend(nrow = 1))
  
  # plot 3
  if(!missing(graph3_input)){
    print("yh")
    graph3_input <- graph3_input[ordering,]
    num_per_sample3 <- rowSums(graph3_input)[ordering]
    
    graph3_melt <- as.data.frame(melt(t(graph3_input))) ; colnames(graph3_melt) <- c("Signature", "Sample", "Number")
    graph3_melt$Total <- rep(as.numeric(num_per_sample3),each=ncol(graph3_input))
    graph3_melt$Proportion <- graph3_melt$Number / graph3_melt$Total
    if(graph3_type == "Number"){
      if(graph3_muttype =="SBS"){
        toobig <- rownames(graph1_input)[rowSums(graph3_input)> 100000]
        hereswhy <- rowSums(graph3_input)[rowSums(graph3_input)> 100000]
        graph3_melt$Total[graph3_melt$Sample%in%toobig] <- 100000
        graph3_melt$Number[graph3_melt$Sample%in%toobig] <- graph3_melt$Proportion[graph3_melt$Sample%in%toobig] *100000
        x_adj <- 875
      }
      maxy3 <- max(graph3_melt$Total)*1.1
      if(graph3_muttype =="DBS")  x_adj <- 537.3626
      if(graph3_muttype =="ID")  x_adj <- 537.3626
    }
    if(graph3_type == "Proportion"){
      maxy3 <- 1
      x_adj <- 537.3626
    }
    
    plot3 <- ggplot(graph3_melt, aes(x = Sample, y = graph3_melt[,graph3_type])) +
      geom_bar(width = 1,stat = "identity", position = "stack", aes(fill = (Signature))) +
      ylab(graph3_type) +
      scale_fill_manual(values = sig_cols)  +
      #facet_grid(Type_f ~ ., scale = "free_y") +
      scale_y_continuous(expand = c(0, 0)) +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 15, hjust = 0),
            legend.spacing.x = unit(0.1, 'cm'),
            plot.margin=unit(c(0, (nsamp/489)*2, 0, nsamp/x_adj), "cm"), # order = top, right, bottom, left
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 15),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0.1, "lines"),
            strip.text.y = element_text(size = 15)
      ) +
      guides(fill = guide_legend(nrow = 1))
  
    if(Waffle==FALSE & !missing(graph3_input)){
      ggarrange(arthur_dend, plot1, plot2,plot3, heights = c(1.7,1, 1, 1), widths = c(.9,.9,.9,.9),
                ncol = 1, nrow = 4)
    }
    }
  
  
  # combine and save plot
  if(Waffle==FALSE & missing(graph3_input)){
    ggarrange(arthur_dend, plot1, plot2, heights = c(1.7, 1, 1), widths = c(.9,.9,.9),
              ncol = 1, nrow = 3)
  }
  
  
  
  if(Waffle==TRUE){
    print(nsamp)
    nsamp <- length(ordering)
    #waf_matrix <- waf_matrix[-which(waf_matrix$Sample %!in% ordering),]
    debut <- 1
    for(i in 1:(nrow(waf_matrix)/nsamp)){
      temporary_waffle <- waf_matrix[debut:(nsamp*i),] 
      waf_matrix[debut:(nsamp*i),] <- temporary_waffle[order(match(temporary_waffle$Sample,ordering)),]
      debut <- 1+(nsamp*i)
    }
    print(nsamp); print(nrow(waf_matrix)); print(head(waf_matrix));print(tail(waf_matrix))
    waf_matrix$x <- c(rep(1:nsamp, length(unique(waf_matrix$category))))
    
    waffle <- ggplot(data = waf_matrix) + 
      geom_tile(aes(x = x, y = y), fill = waf_matrix$colour, color = "white", size = 0.1) +
      scale_x_discrete(expand = c(0, 15)) +
      scale_y_discrete(expand = c(0, 0)) +
      scale_fill_manual(values = (waf_matrix$colour)) +
      annotate("text", x = nsamp+2, y = 6, label = "Project", size = 3, hjust = 0, fontface = "plain") +
      annotate("text", x = nsamp+2, y = 5, label = "Age", size = 3, hjust = 0, fontface = "plain") +
      annotate("text", x = nsamp+2, y = 4, label = "Alcohol", size = 3, hjust = 0, fontface = "plain") +
      annotate("text", x = nsamp+2, y = 3, label = "Tobacco", size = 3, hjust = 0, fontface = "plain") +
      annotate("text", x = nsamp+2, y = 2, label = "Hepatitis.B", size = 3, hjust = 0, fontface = "plain") +
      annotate("text", x = nsamp+2, y = 1, label = "Hepatitis.C", size = 3, hjust = 0, fontface = "plain") +

      
      geom_rect(mapping = aes(xmin = 0.5, xmax = nsamp+0.5, ymin = 0.51, ymax = ncat+0.55), fill = NA, color = "black", size = 0.1) +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            plot.margin= unit(c(0.11, -0.2, 0.1, 0.1), "cm") # order = top, right, bottom, left
      )
    
    ggarrange(arthur_dend, waffle, plot1, plot2, heights = c(1.7,.5, 1, 1), widths = c(.9,.9,.9,.9),
              ncol = 1, nrow = 4)
  }
}

