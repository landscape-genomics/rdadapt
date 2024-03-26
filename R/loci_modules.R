###################################################################################################
##' @name loci_modules
##' @author Thibaut Capblancq
##' 
##' @title Find groups of co-varying loci
##' 
##' @description The \code{loci_modules} function allows the user to find groups of covarying loci on the RDA space.
##' 
##' @param RDA a RDA model from which to extract loci and environmental variable scores
##' @param nb_clusters an integer specifying the number of cluster to identify
##' 
##' @return  
##' 
##' A \code{list} containing :
##' \itemize{
##'   \item \code{loci} : a \code{data.frame} with the cluster associated with each locus
##'   \item \code{polygons} : polygon coordinates to outline the different cluster on a 2D RDA space 
##' }
##' 
##' 
##' @details
##' 
##' Identify a discrete number of adaptive groups across the 2D RDA space.
##' Next version will use more than two dimensions
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

loci_modules <- function(RDA, nb_clusters)
{
  # K-means clustering 
  clusters <- kmeans(abs(RDA_outliers$CCA$v)[, 1:2], nb_clusters, iter.max = 10, nstart = 3)
  
  # Convex hulls.
  df <- data.frame(x = abs(RDA_outliers$CCA$v[, 1])
                   , y = abs(RDA_outliers$CCA$v[, 2])
                   , cluster = clusters$cluster)
  find_hull <- function(df) df[chull(df$x, df$y), ]
  hulls <- ddply(df, "cluster", find_hull)
  
  return(list(loci = df, polygons = hulls))
}
