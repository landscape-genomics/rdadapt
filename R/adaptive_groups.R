###################################################################################################
##' @name adaptive_groups
##' @author Thibaut Capblancq
##' 
##' @title Site clustering
##' 
##' @description This function allows the user to find groups of covarying 
##' sites/individuals/populations within the RDA space.
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
##' @importFrom raster rasterToPoints pointDistance
##' @importFrom factoextra hkmeans
##' @importFrom plyr ddply
##' @importFrom geometry inhulln convhulln
##' 
##' @export
##' 
##'
###################################################################################################

setGeneric("adaptive_groups", def = function(RDA, K, env_pres, nb_clusters) { standardGeneric( "adaptive_groups") })

##'
##' @rdname adaptive_groups
##' @export
##'

setMethod('adaptive_groups', signature(RDA = "rda", env_pres = "missing"), function(RDA, K, nb_clusters)
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

# ##'
# ##' @rdname adaptive_groups
# ##' @export
# ##'
# 
# setMethod('adaptive_groups', signature(RDA = "rda", env_pres = "raster"), function(RDA, K, env_pres, nb_clusters)
# {
#   ## CHECKS -------------------------------------------------------------------
#   if (!inherits(RDA, "rda")) { stop("\n RDA must be a 'rda' object") }
#   if (K > ncol(RDA$CCA$v)) {
#     K = ncol(RDA$CCA$v)
#     warning("\n Not enough RDA axis available, K is set to ncol(RDA$CCA$v)")
#   }
#   
#   ## FUNCTION -----------------------------------------------------------------
#   
#   ## Make predictions
#   AI_pres <- adaptive_index(RDA = RDA
#                             , K = K
#                             , env = env_pres
#                             , env_mask = NULL
#                             , method = "loadings")
#   
#   # RDA1 and RDA2 scores for all the pixels of the range
#   tmp_df <- data.frame(x = rasterToPoints(res_RDA_proj_current$RDA1)[, 1]
#                     , y = rasterToPoints(res_RDA_proj_current$RDA1)[, 2]
#                     , RDA1 = rasterToPoints(res_RDA_proj_current$RDA1)[, 3]
#                     , RDA2 = rasterToPoints(res_RDA_proj_current$RDA2)[, 3])
#   
#   ## Hierarchical K-means clustering with K number of clusters
#   clust <- hkmeans(x = tmp_df
#                    , k = nb_clusters
#                    , hc.metric = "euclidean"
#                    , hc.method = "ward.D2"
#                    , iter.max = 10
#                    , km.algorithm = "Hartigan-Wong")
#   
# 
#   return(clusterNA)
# })
# 
# 
# adaptive_groups <- function(RDA, K)
# {
#   # RDA1 and RDA2 scores for all the pixels of the range
#   TAB <- data.frame(x = rasterToPoints(res_RDA_proj_current$RDA1)[, 1]
#                     , y = rasterToPoints(res_RDA_proj_current$RDA1)[, 2]
#                     , RDA1 = rasterToPoints(res_RDA_proj_current$RDA1)[, 3]
#                     , RDA2 = rasterToPoints(res_RDA_proj_current$RDA2)[, 3])
#   
#   # Hierarchical K-means clustering
#   TAB_sample <- TAB[, 3:4]
#   clust <- list()
#   for (i in 1:15) {
#     clust[[i]] <- hkmeans(x = TAB_sample, k = i, hc.metric = "euclidean"
#                           , hc.method = "ward.D2", iter.max = 10, km.algorithm = "Hartigan-Wong")
#   }
#   
#   # Select the best number of clusters
#   within_cust_var <- lapply(clust, function(x) x$withinss / x$totss)
#   plot(unlist(lapply(within_cust_var, mean)))
#   
#   # Extrapolate the clustering to all the pixels with 4 clusters
#   df <- data.frame(x = TAB_sample$RDA1, y = TAB_sample$RDA2, cluster = clust[[4]]$cluster) 
#   find_hull <- function(df) df[chull(df$x, df$y), ]
#   hulls <- ddply(df, "cluster", find_hull) # Convex hulls.
#   cluster <- rep(NA, nrow(TAB))
#   for (i in 1:4) {
#     inhull <- inhulln(convhulln(hulls[hulls$cluster == i, 1:2])
#                       , as.matrix(TAB[, 3:4])) # Intersection convex hulls and points
#     cluster[inhull] <- i
#   }
#   clusterNA <- cluster[which(is.na(cluster))]
#   for (i in 1:length(clusterNA)) {
#     clusterNA[i] <- which.min(pointDistance(TAB[which(is.na(cluster))[i], 3:4]
#                                             , clust[[4]]$centers, lonlat = FALSE)) # Fill the NA using the distance to the hull centers
#   }
#   cluster[which(is.na(cluster))] <- clusterNA
#   
#   return(clusterNA)
# }
