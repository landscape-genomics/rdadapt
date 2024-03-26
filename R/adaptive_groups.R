###################################################################################################
##' @name adaptive_groups
##' @author Thibaut Capblancq
##' 
##' @title Find groups of co-varying individuals/populations
##' 
##' @description The \code{adaptive_groups} function allows the user to find groups of covarying 
##' individuals/populations based on their scores in the RDA space and a hkmean clustering procedure. 
##' \code{adaptive_groups} accommodates stack of rasters or discrete populations.
##' 
##' @param RDA a RDA model from which to extract loci and environmental variable scores
##' @param K an integer specifying the number of RDA axes to use for the clustering
##' @param env a stack of rasters or a data.frame with the environmental conditions of the sites to predict
##' @param nb_clusters an integer specifying the number of cluster to identify
##' 
##' @return  
##' 
##' A \code{list} containing :
##' \itemize{
##'   \item \code{samples} : a \code{data.frame} with the cluster associated with each individual/population
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
##' @importFrom raster rasterToPoints pointDistance
##' @importFrom factoextra hkmeans
##' @importFrom plyr ddply
##' @importFrom geometry inhulln convhulln
##' 
##' @export
##' 
##'
###################################################################################################

setGeneric("adaptive_groups", def = function(RDA, K, env, nb_clusters) { standardGeneric( "adaptive_groups") })

##'
##' @rdname adaptive_groups
##' @export
##'

setMethod('adaptive_groups', signature(RDA = "rda", env = "missing"), function(RDA, K, nb_clusters)
{
  ## CHECKS -------------------------------------------------------------------
  if (!inherits(RDA, "rda")) { stop("\n RDA must be a 'rda' object") }
  if (K > ncol(RDA$CCA$v)) {
    K = ncol(RDA$CCA$v)
    warning("\n Not enough RDA axis available, K is set to ncol(RDA$CCA$v)")
  }
  
  ## FUNCTION -----------------------------------------------------------------
  
  # K-means clustering using K number of axes and nb_clusters number of cluster
  clusters <- kmeans(RDA$CCA$u[, 1:K]
                     , nb_clusters
                     , iter.max = 10
                     , nstart = 3)
  
  # Convex hulls
  df <- data.frame(x = RDA$CCA$u[, 1]
                   , y = RDA$CCA$u[, 2]
                   , cluster = clusters$cluster)
  find_hull <- function(df) df[chull(df$x, df$y), ]
  hulls <- ddply(df, "cluster", find_hull)
  
  return(list(samples = df, polygons = hulls))
})

##'
##' @rdname adaptive_groups
##' @export
##'

setMethod('adaptive_groups', signature(RDA = "rda", env = "SpatRaster"), function(RDA, K, env, nb_clusters)
{
  ## CHECKS -------------------------------------------------------------------

  # @Maya Do we need to check again if adaptive_index below does it anyway?

  ## FUNCTION -----------------------------------------------------------------

  ## Make predictions
  AI_pres <- adaptive_index(RDA = RDA
                            , K = K
                            , env = env
                            , env_mask = NULL
                            , method = "loadings")

  # RDA1 and RDA2 scores for all the pixels of the range
  tmp_df <- data.frame(x = as.data.frame(AI_pres$RDA1, xy = TRUE)[, 1]
                    , y = as.data.frame(AI_pres$RDA1, xy = TRUE)[, 2]
                    , RDA1 = as.data.frame(AI_pres$RDA1, xy = TRUE)[, 3]
                    , RDA2 = as.data.frame(AI_pres$RDA2, xy = TRUE)[, 3])

  ## Hierarchical K-means clustering with K number of clusters
  clusters <- hkmeans(x = tmp_df
                   , k = nb_clusters
                   , hc.metric = "euclidean"
                   , hc.method = "ward.D2"
                   , iter.max = 10
                   , km.algorithm = "Hartigan-Wong")

  # Convex hulls
  tmp_df$cluster <- clusters$cluster
  find_hull <- function(tmp_df) tmp_df[chull(tmp_df$x, tmp_df$y), ]
  hulls <- ddply(tmp_df, "cluster", find_hull)

  return(list(samples = tmp_df, polygons = hulls, variance_within = clusters$withinss/clusters$totss))
})
