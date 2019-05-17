#' palimpsest_distCosine
#'
#' Function to calculate cosine distance
#' @param m input
#'
#' @export
#' @importFrom lsa cosine

palimpsest_distCosine <- function(m)
{
  nsamp <- nrow(m)
  res <- matrix(NA,nrow=nsamp,ncol=nsamp)
  rownames(res) <- colnames(res) <- rownames(m)
  for(i in 1:nsamp)
  {
    for(j in 1:nsamp)
    {
      res[i,j] <- cosine(m[i,],m[j,])
    }
  }
  as.dist(1-res)
}
