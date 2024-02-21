###################################################################################################
##' @name loci_clustering
##' @author Thibaut Capblancq
##' 
##' @title Loci clustering
##' 
##' @description This function allows the user to find groups of covarying loci on the RDA space.
##' 
##' @param RDA blabla
##' @param K blabla
##' 
##' @return  
##' 
##' A \code{data.frame} containing :
##' \itemize{
##'   \item \code{p.values} : blabla
##'   \item \code{q.values} : blabla
##' }
##' 
##' 
##' @details
##' 
##' Blablabla
##' 
##' 
##' @keywords 
##' 
##' @seealso 
##' 
##' @examples
##' 
##' 
##' @importFrom stats kmeans
##' @importFrom plyr ddply
##' 
##' @export
##' 
##'
###################################################################################################

loci_clustering <- function(RDA, K)
{
  # K-means clustering with K = 2
  clusters <- kmeans(abs(RDA_outliers$CCA$v)[, 1:2], 2, iter.max = 10, nstart = 3)
  
  # Convex hulls.
  df <- data.frame(x = abs(RDA_outliers$CCA$v[, 1])
                   , y = abs(RDA_outliers$CCA$v[, 2])
                   , cluster = clusters$cluster)
  find_hull <- function(df) df[chull(df$x, df$y), ]
  hulls <- ddply(df, "cluster", find_hull)
  
  return(hulls)
}
