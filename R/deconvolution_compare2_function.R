#' deconvolution_compare2
#'
#' Function to calculate cosine similarity scores between two sets of mutational signatures (e.g. SBS COSMIC signatures vs Palimpsest de novo signatures).
#' @param new_signatures Data frame of de novo extracted mutational signatures
#' @param COSMIC_Signatures Data frame of reference mutational signatures (The reference SBS_cosmic, DBS_cosmic & ID_cosmic matrices are taken from Alexandrov et al. (2018)).
#'
#' @export
#' @import lsa
#' @examples
#' SBS_cosine_similarities <- deconvolution_compare2(SBS_denovo_sigs,SBS_cosmic)

deconvolution_compare2 <- function(new_signatures, COSMIC_Signatures) {
  if(any(!grepl("denovo",rownames(new_signatures),ignore.case = T))) rownames(new_signatures) <- rep(paste("DeNovo", rownames(new_signatures), sep = "_"))
  mutmat <- t(rbind(new_signatures, COSMIC_Signatures))
  m <- as.matrix(palimpsest_distCosine(t(mutmat)))
  row_distance = as.dist(palimpsest_distCosine(t(m)))
  row_cluster = hclust(row_distance, method = "ward.D")
  col_distance = as.dist(palimpsest_distCosine(t(m)))
  col_cluster = hclust(col_distance, method = "ward.D")
  heatmap.2(1 - m, key.title = "Cosine similarity", keysize = 0.8, 
            main = " Cosine similarity Matrix", notecol = "black", 
            density.info = "none", trace = "none", margins = c(12, 
                                                               9), col = colorRampPalette(c("white", "yellow", "red"))(n = 299), 
            Rowv = as.dendrogram(row_cluster), Colv = as.dendrogram(col_cluster))
  return(1 - m)
}