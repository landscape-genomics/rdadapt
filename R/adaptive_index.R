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
##' @importFrom foreach foreach %do%
##' @importFrom terra crs rast as.data.frame mask
##' @importFrom vegan predict
##' 
##' @export
##' 
##'
###################################################################################################


adaptive_index <- function(rda, K, env_pres, range = NULL, method = "loadings", scale_env = NULL, center_env = NULL)
{
  ## Checks
  # inherits(rda, "rda")
  # "CCA" %in% names(rda) ## should not be necessary ?
  # "biplot" %in% names(rda$CCA) ## should not be necessary ?
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
  
  
  ## Get RDA informations -----------------------------------------------------
  rda_biplot <- rda$CCA$biplot
  var_names <- row.names(rda_biplot)
  
  ## Transform environmental raster for prediction
  env_crs <- crs(env_pres)
  env_df <- as.data.frame(env_pres[[var_names]], xy = TRUE)
  env_xy <- env_df[, c("x", "y")]
  env_var <- env_df[, var_names]
  
  if (!(is.null(scale_env) & is.null(center_env))) { ## Standardize environmental variables
    # env_var <- as.data.frame(scale(env_var, center_env[var_names], scale_env[var_names])) ## MARCHE PAS
    env_var <- as.data.frame(scale(env_var, center_env, scale_env)) ## one value only
  }
  
  ## MAKE PREDICTIONS ---------------------------------------------------------
  if (method == "predict") {
    pred <- predict(rda, env_var, type = "lc")
  }
  Proj_pres <- foreach(i = 1:K) %do%
    {
      if (method == "loadings") { ## Predict pixels genetic component based on RDA axes
        tmp_df <- data.frame(env_xy, z = as.vector(apply(env_var,  1, function(x) sum(x * rda_biplot[, i]))))
      } else if (method == "predict") { ## Predict with RDA model and linear combinations
        tmp_df <- data.frame(env_xy, z = as.vector(pred[, i]))
      }
      ras_pres <- rast(tmp_df, type = "xyz", crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA", as.character(i))
      return(ras_pres)
    }
  Proj_pres <- rast(Proj_pres)
  
  if (!is.null(range)) { ## Mask with range
    Proj_pres <- mask(Proj_pres, range)
  }
  
  
  ## Return predictions for each RDA axis
  return(Proj_pres)
}
