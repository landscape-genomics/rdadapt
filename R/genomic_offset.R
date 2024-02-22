###################################################################################################
##' @name genomic_offset
##' @author Thibaut Capblancq
##' 
##' @title Genomic offset
##' 
##' @description This function allows the user to predict genomic offset from a RDA model.
##' 
##' @param RDA blabla
##' @param K blabla
##' @param env_pres blabla
##' @param env_fut blabla
##' @param range (\emph{optional, default} \code{NULL}) \cr blabla
##' @param method  (\emph{default} \code{'loadings'}) \cr blabla
##' @param scale_env blabla
##' @param center_env blabla
##' 
##' 
##' @return  
##' 
##' A \code{list} containing blabla
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
##' @importFrom stats dist
##' @importFrom raster rasterToPoints rasterFromXYZ mask
##' @importFrom XX predict
##' 
##' @export
##' 
##'
###################################################################################################


genomic_offset <- function(rda, K, env_pres, env_fut, range = NULL, method = "loadings", scale_env = NULL, center_env = NULL)
{
  ## Checks
  # inherits(rda, "rda")
  # "CCA" %in% names(rda) ## should not be necessary ?
  # "biplot" %in% names(rda$CCA) ## should not be necessary ?
  # "eig" %in% names(rda$CCA) ## should not be necessary ?
  # K <= ncol(rda$CCA$biplot)
  # inherits(env_pres, c("SpatRaster", "RasterLayer", "RasterStack"))
  # row.names(rda_biplot) %in% names(env_pres)
  # method %in% c("loadings", "predict")
  # length(scale_env) == nrow(rda_biplot)
  # length(center_env) == nrow(rda_biplot)
  # names(scale_env) == row.names(rda_biplot)
  # names(center_env) == row.names(rda_biplot)
  # inherits(range, c("SpatRaster", "RasterLayer", "RasterStack"))
  # nlyr(range) == 1 ?
  
  ## MAKE PREDICTIONS ---------------------------------------------------------
  AI_pres <- adaptive_index(rda, K, env_pres, range, method, scale_env, center_env)
  AI_fut <- adaptive_index(rda, K, env_fut, range, method, scale_env, center_env)
  
  ## Single axis genetic offset 
  offset <- foreach(i = 1:K) %do%
    {
      ras <- abs(AI_pres[[i]] - AI_fut[[i]])
      names(ras) <- paste0("RDA", i)
      return(ras)
    }
  offset <- rast(offset)
  
  ## Weight current and future adaptive indices based on eigen values of associated axes
  AI_pres_w <- AI_pres * weights
  AI_fut_w <- AI_fut * weights
  
  ## Predict a global genetic offset
  offset_global <- offset[[1]]
  offset_global[!is.na(offset_global)] <- sapply(1:ncell(AI_pres_w), function(x) {
    tmp_mat <- rbind(AI_pres_w[x], AI_fut_w[x])
    return(dist(tmp_mat, method = "euclidean"))
  })
  names(offset_global) <- "Global_offset"
  
  
  return(list(Proj_pres = AI_pres,
              Proj_fut = AI_fut,
              Proj_offset = offset,
              Proj_offset_global = offset_global,
              weights = weights))
}
