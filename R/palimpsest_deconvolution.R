#' makeMutypeMatFromVcf
#'
#' Function to create a matrix in mutation type x sample format with either counts or proportions
#' @param vcf vcf data frame containing the mutations/SVs
#' @param sample.col sample.col
#' @param mutcat.col mutcat.col
#' @param mutypes mutypes
#' @param proportion proportion
#'
#' @export


makeMutypeMatFromVcf <- function (vcf,
                                  sample.col = "sample",
                                  mutcat.col = "mutcat3",
                                  mutypes = c("CA", "CG", "CT", "TA", "TC", "TG"),
                                  proportion = TRUE)
{
  tmp <- split(vcf, vcf[, sample.col])
  tmp <- lapply(tmp, function(d) {
    sapply(mutypes, function(m) {
      sum(d[, mutcat.col] == m, na.rm = T)
    })
  })
  tmp <- as.matrix(as.data.frame(tmp))
  if (proportion) {
    for (j in 1:ncol(tmp)) tmp[,j] <- tmp[, j]/sum(tmp[,j])
  }
  tmp
}

#' palimpsestInput
#'
#' Function to create a matrix in mutation type x sample format with either counts or proportions
#' @param vcf vcf data frame containing the mutations/SVs
#' @param type SNV for single nucleotide variant signatures, SV for structural variant signatures
#' @param sample.col Sample name column in vcf
#' @param mutcat.col Mutation category column name in vcf
#' @param proportion If TRUE, the output matrix will indicate mutation proportions instead of numbers
#'
#' @export
#' @import reshape2

palimpsestInput <- function(vcf,
						   type="SNV",
                           sample.col = "Sample",
                           mutcat.col = "mutcat3",
                           proportion = TRUE)
{
  if(type=="SNV")
  {
    vcf <- vcf[which(vcf$Type=="SNV"),]
  	bases <- c("A", "C", "G", "T")
  	ctxt16 <- paste(rep(bases, each = 4), rep(bases, 4), sep = ".")
  	mt <- c("CA", "CG", "CT", "TA", "TC", "TG")
  	types96 <- paste(rep(mt, each = 16), rep(ctxt16, 6), sep = "_")
  	types96 <- sapply(types96, function(z) { sub("\\.", substr(z, 1, 1), z)})
    vcf_props <- as.matrix(makeMutypeMatFromVcf(vcf = vcf,sample.col = sample.col,mutcat.col = mutcat.col,mutypes = types96,proportion = proportion))
    df = data.frame(matrix(vector(), 96, dim(vcf_props)[2],
                           dimnames=list(c(), colnames(vcf_props))),
                    stringsAsFactors=F)
    rownames(df) <- names(types96)
    data <- vcf_props[match(rownames(df),rownames(vcf_props)),]

  }
  if(type=="SV")
  {
    tt <- as.matrix(table(vcf[,mutcat.col],vcf[,sample.col]))
    mat <- matrix(as.numeric(tt),nrow=nrow(tt),ncol=ncol(tt));rownames(mat) <- rownames(tt);colnames(mat) <- colnames(tt)
    data <- mat
    categs <- c("DEL_0-1kb_clust_1", "DEL_1-10kb_clust_1", "DEL_10-100kb_clust_1",
    "DEL_100kb-1Mb_clust_1", "DEL_1-10Mb_clust_1", "DEL_>10Mb_clust_1",
    "DUP_0-1kb_clust_1", "DUP_1-10kb_clust_1", "DUP_10-100kb_clust_1",
    "DUP_100kb-1Mb_clust_1", "DUP_1-10Mb_clust_1", "DUP_>10Mb_clust_1",
    "INV_0-1kb_clust_1", "INV_1-10kb_clust_1", "INV_10-100kb_clust_1",
    "INV_100kb-1Mb_clust_1", "INV_1-10Mb_clust_1", "INV_>10Mb_clust_1",
    "BND_clust_1", "DEL_0-1kb_clust_0", "DEL_1-10kb_clust_0",
    "DEL_10-100kb_clust_0", "DEL_100kb-1Mb_clust_0", "DEL_1-10Mb_clust_0",
    "DEL_>10Mb_clust_0", "DUP_0-1kb_clust_0", "DUP_1-10kb_clust_0",
    "DUP_10-100kb_clust_0", "DUP_100kb-1Mb_clust_0", "DUP_1-10Mb_clust_0",
    "DUP_>10Mb_clust_0", "INV_0-1kb_clust_0", "INV_1-10kb_clust_0",
    "INV_10-100kb_clust_0", "INV_100kb-1Mb_clust_0", "INV_1-10Mb_clust_0",
    "INV_>10Mb_clust_0", "BND_clust_0")
    missings_cat = setdiff(categs, rownames(data))
    for(i in missings_cat){
      data = rbind(data, assign(i, rep(0, ncol(data))))
      rownames(data)[nrow(data)] = i
    }
    if(proportion){ data <- prop.table(mat,2)}
    }
    return(data)
  }

#' deconvolution_fit
#'
#' Function to identify the weighted combination of the input signatures explaining each tumor's mutational profile
#' @param vcf vcf data frame containing the mutations/SVs
#' @param type SNV for single nucleotide variant signatures, SV for structural variant signatures
#' @param input_data Matrix in mutation type x sample format in proportions
#' @param threshold Discard signatures contributing less then this percent of total mutations within each sample
#' @param input_signatures Data frame describing the mutational signatures to fit within the provided cohort of samples
#' @param sig_cols Character vector indicating the colors representing each signature in graphical outputs. Must match to the total number of provided signatures
#' @param plot Logical indicating whether graphical outputs should be generated
#' @param resdir Result directory
#'
#' @export
#' @import NMF

deconvolution_fit <- function (vcf = vcf, type = NULL, input_data = data,
                               threshold = 5, input_signatures = COSMIC_Signatures,
                               sig_cols = mycol, plot = TRUE, resdir = resdir)
{
  "%ni%" <- Negate("%in%");
  requireNamespace("NMF", quietly = TRUE)
  mutSign_props <- c()
  for (s in unique(vcf[, "Sample"])) {
    print(s)
    vcf. <- vcf[which(vcf[, "Sample"] == s), ]
    if (plot == TRUE) {
      resdir.. <- file.path(resdir, s)
      if (!file.exists(resdir..)) {
        dir.create(resdir..)
      }
      if (type == "SNV") {
      	vcf. <- vcf.[which(vcf.$Type=="SNV"),]
        plot6mutationSpectrumFromVcf(vcf., sample.col = "Sample",plot.file = file.path(resdir.., "Mean_proportion_6_substitution_types.pdf"))
        plot96mutationSpectrumFromVcf(vcf., sample.col = "Sample",plot.file = file.path(resdir.., "Mean_proportion_96_substitution_types.pdf"))
        plotStrandBias6types(vcf., plot.file = file.path(resdir..,"Strand_bias_6_substitution_types.pdf"))
        plotStrandBias96types(vcf., plot.file = file.path(resdir..,"Strand_bias_96_substitution_types.pdf"))
      }
      if (type == "SV") {
        datas_for_plot <- generate_table(vcf)
        translocs = datas_for_plot[[1]]
        others = datas_for_plot[[2]]
        df <- t(as.data.frame(input_data[, s]))
        rownames(df) <- s
        pdf(file.path(resdir.., "Mean_proportion_38_sv_types.pdf"),
            width = 24, height = 6)
        plot.SV.sigs(df)
        dev.off()
        make_plot(translocs, others, s,resdir = resdir..)
      }
    }
    res <- fcnnls(as.matrix(t(input_signatures)), input_data[,s], verbose = TRUE, pseudo = FALSE)
    sig.tot <- margin.table(res$x, 2)
    num.vec <- as.numeric(res$x)
    sig.fit <- data.frame(res$x)
    sig.fit$percent.fit <- sig.fit$res.x/sig.tot * 100
    sig.fit$sig_fit <- sig.fit$percent.fit
    sig.fit$sig_fit[sig.fit$sig_fit < threshold] <- 0
    mat.sigs <- data.frame(t(sig.fit$sig_fit))
    colnames(mat.sigs) <- rownames(res$x)
    mat.sigs$sums <- rowSums(mat.sigs)
    prop.sigs <- (mat.sigs/mat.sigs$sum) * 100
    prop.sigs <- prop.sigs[, -dim(prop.sigs)[2]]
    rownames(prop.sigs) <- s
    signature_content <- prop.sigs/100
    prop.sigs. <- data.frame(t(prop.sigs))
    prop.sigs. <- subset(prop.sigs., prop.sigs.[, 1] > 0)
    if (plot == TRUE) {
      pdf(file.path(resdir.., paste(s, "_Signature_Contribution.pdf", sep = "")), width = 12, height = 10)
      pie(t(prop.sigs.), labels = rownames(prop.sigs.),
          explode = 0.1, main = paste("Mutational Signatures Contribution in: ", colnames(prop.sigs.), sep = ""), col = sig_cols[rownames(prop.sigs.)],
          border = sig_cols[rownames(prop.sigs.)])
      dev.off()
    }
    mutSign_props <- rbind(mutSign_props, signature_content)
  }
  mutSign_nums <- mutSign_props
  tot.muts <- data.frame(table(vcf[, "Sample"]))
  mutSign_nums$total_mutations <- tot.muts[match(rownames(mutSign_nums),
                                                 tot.muts$Var1), "Freq"]
  mutSign_nums <- round(mutSign_nums * mutSign_nums$total_mutations)
  mutSign_nums <- mutSign_nums[-dim(mutSign_nums)[2]]
  outputs <- list(sig_props = mutSign_props, sig_nums = mutSign_nums)
  return(outputs)
}

#' deconvolution_nmf
#'
#' Function to perform de novo mutational signature extraction and represent the extracted signatures.
#' @param input_data Matrix in mutation type x sample format in proportions
#' @param type SNV for single nucleotide variant signatures, SV for structural variant signatures
#' @param range_of_sigs numeric Range of possible numbers of mutational signatures to extract
#' @param num_of_sigs Number of mutational signatures to extract. If left to default "auto" value, the appropriate number of signatures will be estimated from NMF metrics.
#' @param nrun Number of iterations to perform for each number of signatures
#' @param method Specification of the NMF algorithm. ‘brunet’ method corresponds to the standard NMF algorithm from Brunet,2004
#' @param plot_sigs If TRUE, will output the extracted mutational signature plots
#' @param resdir Result directory
#'
#' @export
#' @import NMF
#' @import registry
#' @import doParallel

deconvolution_nmf <- function (input_data = NULL, type = NULL, range_of_sigs = NULL, num_of_sigs ="auto",
                               nrun = nrun, method = "brunet", plot_sigs = TRUE, resdir = resdir)
{
  requireNamespace("NMF", quietly = TRUE)
  sumRows <- rowSums(input_data);sort(sumRows);zeroes <- which(sumRows==0)
  if(length(zeroes)){input_data[zeroes,1] <- 1e-10}
  
  if(num_of_sigs=="auto"){
  print("Estimating the optimal number of mutational signatures...")
  estimate <- nmfEstimateRank(x = input_data, range_of_sigs,
                              method = method, nrun = nrun, seed = 123456)
  p <- plot(estimate, y = NULL,
            what = "all",
            na.rm = FALSE, xname = "x", yname = "y",
            xlab = "Factorization rank", ylab = "",
            main = "NMF rank survey")
  pdf(file.path(resdir, "NMF_Rank_Estimates.pdf"),width = 8,height = 6)
  print(p)
  dev.off()
  z <- estimate$measures$cophenetic[which(diff(estimate$measures$cophenetic) < 0)]
  if(length(z)==1) steep_index <- which(estimate$measures$cophenetic==z)
  if(length(z)> 1) steep_index <- which(estimate$measures$cophenetic==z[ which(abs(diff(z))==max(abs(diff(z))) )]); 
  steep_index <- steep_index[length(steep_index)]
  estimated_rank <- estimate$measures$rank[steep_index]
  }
  else{
    estimated_rank <- num_of_sigs
  }
  print(paste("Calculating exposures of", estimated_rank, "signatures in the input tumors",
              sep = " "))
  res <- nmf(input_data, rank = estimated_rank, method = method,
             nrun = nrun, seed = 123456)
  sigs = t(basis(res))
  rownames(sigs) <- paste("Signature.", 1:nrow(sigs), sep = "")
  spec <- sigs/rowSums(sigs)
  if (plot_sigs == TRUE) {
    if (type == "SV") {
      pdf(file.path(resdir, "Mutational_Signatures.pdf"),width = 24, height = 5)
      plot.SV.sigs(spec)
      dev.off()
    }
    if (type == "SNV") {
      pdf(file.path(resdir, "Mutational_Signatures.pdf"),width = 24, height = 5)
      plot.SNV.sigs(spec)
      dev.off()
    }
  }
  return(spec)
}

#' deconvolution_compare
#'
#' Function to calculate cosine similarity scores between two sets of mutational signatures.
#' @param new_signatures Data frame of de-novo extracted mutational signatures
#' @param COSMIC_Signatures Data frame of already published mutational signatures (Example: Alexandrov et.al,2013)
#'
#' @export
#' @import lsa

deconvolution_compare <- function (new_signatures, COSMIC_Signatures)
{
  rownames(new_signatures) <- rep(paste("DeNovo", rownames(new_signatures),sep = "_"))
  mutmat <- t(rbind(new_signatures, COSMIC_Signatures))
  m <- as.matrix(palimpsest_distCosine(t(mutmat)))
  row_distance = as.dist(palimpsest_distCosine(t(m)))
  row_cluster = hclust(row_distance, method = "ward.D")
  col_distance = as.dist(palimpsest_distCosine(t(m)))
  col_cluster = hclust(col_distance, method = "ward.D")
  heatmap.2(1-m,
            key.title = "Cosine similarity", keysize = 0.8,
            main = " Cosine similarity Matrix",
            notecol = "black",
            density.info = "none",
            trace = "none",
            margins = c(12,9),
            col = colorRampPalette(c("white","yellow", "red"))(n = 299),
            Rowv = as.dendrogram(row_cluster),
            Colv = as.dendrogram(col_cluster))
	return(1-m)
}


#' palimpsestOrigin
#'
#' Function to estimate the probability of each individual somatic alteration being due to each mutational process
#' @param vcf data frame containing the mutations/SVs
#' @param type SNV for single nucleotide variant signatures, SV for structural variants signatures
#' @param sample.col Sample column name in vcf
#' @param mutcat.col Mutation category column name in vcf
#' @param signature_contribution Matrix in sample x mutational signature exposure format in proportions
#' @param input_signatures Data frame describing the mutational signatures to fit within the provided cohort of samples
#'

#' @export

palimpsestOrigin <- function(vcf=vcf, type = "SV",
                              sample.col="Sample", mutcat.col="Category1",
                              signature_contribution=nmf_signatures$sig_nums,
                              input_signatures=COSMIC_Signatures){
  if(type == "SNV"){
    bases <- c("A", "C", "G", "T")
    ctxt16 <- paste(rep(bases, each = 4), rep(bases, 4), sep = ".")
    mt <- c("CA", "CG", "CT", "TA", "TC", "TG")
    types <- paste(rep(mt, each = 16), rep(ctxt16, 6), sep = "_")
    types <- sapply(types, function(z) {
      sub("\\.", substr(z, 1, 1), z)})
  }else if (type == "SV"){
    types = c("BND_clust_0","BND_clust_1","DEL_0-1kb_clust_0","DEL_0-1kb_clust_1","DEL_100kb-1Mb_clust_0",
              "DEL_100kb-1Mb_clust_1","DEL_10-100kb_clust_0","DEL_10-100kb_clust_1","DEL_>10Mb_clust_0",
              "DEL_>10Mb_clust_1","DEL_1-10kb_clust_0","DEL_1-10kb_clust_1","DEL_1-10Mb_clust_0",
              "DEL_1-10Mb_clust_1","DUP_0-1kb_clust_0","DUP_0-1kb_clust_1","DUP_100kb-1Mb_clust_0",
              "DUP_100kb-1Mb_clust_1","DUP_10-100kb_clust_0","DUP_10-100kb_clust_1","DUP_>10Mb_clust_0",
              "DUP_>10Mb_clust_1","DUP_1-10kb_clust_0","DUP_1-10kb_clust_1","DUP_1-10Mb_clust_0",
              "DUP_1-10Mb_clust_1","INV_0-1kb_clust_0","INV_0-1kb_clust_1","INV_100kb-1Mb_clust_0",
              "INV_100kb-1Mb_clust_1","INV_10-100kb_clust_0","INV_10-100kb_clust_1","INV_>10Mb_clust_0",
              "INV_>10Mb_clust_1","INV_1-10kb_clust_0","INV_1-10kb_clust_1","INV_1-10Mb_clust_0","INV_1-10Mb_clust_1")
  }else{
    print("Error in palimpsestOrigin(): type of analyse must be SNV or SV")
    return(0)
  }
  contrib <- t(signature_contribution)
  sigs <- t(input_signatures)
  rownames(sigs) <- types
  vcf[,paste(rownames(contrib),"prob",sep=".")] <- NA
  print("Analyzing probabilities of the most likely mutational signatures:")
  total <- length(unique(vcf[,sample.col]))
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  counter_start <- c()
  for(samp in unique(vcf[,sample.col])){
    counter_start <- match(samp,unique(vcf[,sample.col]))
    Sys.sleep(0.1)
    setTxtProgressBar(pb, counter_start)
    for(mut_type in types){
      probs <- round(prop.table(sigs[mut_type,]*contrib[,samp]),digits=4)
      if(all(is.na(probs))){probs <- prop.table(contrib[,samp])} # special case where all active signatures in the samples have a null probability for a given mutation category
      ind <- which(vcf[,sample.col]==samp & vcf[,mutcat.col]==mut_type);length(ind)
      for(sig in names(probs)){
        vcf[ind,paste(sig,"prob",sep=".")] <- probs[sig]
      }
    }
  }
  close(pb)
  print("Assigning the origin of mutations:")
  vcf$Sig.max <- sapply(1:nrow(vcf),function(i){names(which.max(vcf[i,paste(rownames(contrib),"prob",sep=".")]))})
  vcf$Sig.max <- substr(vcf$Sig.max, 1, nchar(vcf$Sig.max)-5)
  return(vcf)
}
