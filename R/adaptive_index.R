###################################################################################################
##' @name adaptive_index
##' @author Thibaut Capblancq
##' 
##' @title Adaptive index
##' 
##' @description This function allows the user to project the adaptive component turnover across 
##' a landscape.
##' 
##' @param RDA blabla
##' @param K blabla
##' @param env_pres blabla
##' @param range (\emph{optional, default} \code{NULL}) \cr blabla
##' @param method  (\emph{default} \code{'loadings'}) \cr blabla
##' @param scale_env blabla
##' @param center_env blabla
##' 
##' @return  
##' 
##' A \code{matrix} containing blabla
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
##' @importFrom raster rasterToPoints rasterFromXYZ mask
##' @importFrom XX predict
##' 
##' @export
##' 
##'
###################################################################################################

adaptive_index <- function(RDA, K, env_pres, range = NULL, method = "loadings", scale_env = NULL, center_env = NULL)
{
  
  # Formatting environmental rasters for projection
  var_env_proj_pres <- as.data.frame(env_pres[[row.names(RDA$CCA$biplot)]], xy = TRUE)
  
  # Standardization of the environmental variables if necessary
  var_env_proj_RDA <- as.data.frame(var_env_proj_pres[, -c(1, 2)])
  if (!(is.null(scale_env)&is.null(center_env))) {
  var_env_proj_RDA <- as.data.frame(scale(var_env_proj_RDA, 
                                          center_env[row.names(RDA$CCA$biplot)], 
                                          scale_env[row.names(RDA$CCA$biplot)]))
  }
  
  # Predicting pixels genetic component based on RDA axes
  Proj_pres <- list()
  if (method == "loadings") {
    for (i in 1:K) {
      ras_pres <- rast(data.frame(var_env_proj_pres[, c(1, 2)], 
                                  as.vector(apply(var_env_proj_RDA[,names(RDA$CCA$biplot[, i])], 1, function(x) sum(x * RDA$CCA$biplot[, i])))), 
                       type="xyz",
                       crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Prediction with RDA model and linear combinations
  if (method == "predict") {
    pred <- predict(RDA, var_env_proj_RDA[, names(RDA$CCA$biplot[, i])],  type = "lc")
    for (i in 1:K) {
      ras_pres <- rast(data.frame(var_env_proj_pres[, c(1, 2)],
                                  as.vector(pred[, i])),
                       type="xyz",
                       crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Mask with the range if supplied
  if (!is.null(range)) {
    Proj_pres <- lapply(Proj_pres, function(x) mask(x, range))
  }
  
  # Returning projections for current climates for each RDA axis
  return(Proj_pres = Proj_pres)
}
