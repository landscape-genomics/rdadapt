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
##' @param env blabla
##' @param env_mask (\emph{optional, default} \code{NULL}) \cr blabla
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


adaptive_index <- function(RDA, K, env, env_mask = NULL, method = "loadings", scale_env = NULL, center_env = NULL)
{
  ## CHECKS -------------------------------------------------------------------
  if (!inherits(RDA, "rda")) { stop("\n RDA must be a 'rda' object") }
  if (!("CCA" %in% names(RDA)) || !("biplot" %in% names(RDA$CCA))) {
    stop("\n RDA$CCA$biplot seems not to exist")
  }
  if (K > ncol(RDA$CCA$biplot)) {
    K = ncol(RDA$CCA$biplot)
    warning("\n Not enough RDA axis available, K is set to ncol(RDA$CCA$biplot)")
  }
  if (!inherits(env, c("SpatRaster", "RasterLayer", "RasterStack"))) {
    stop("\n env must be a 'SpatRaster', 'RasterLayer' or 'RasterStack' object")
  } else if (inherits(env, c("RasterLayer", "RasterStack"))) {
    env <- rast(env)
  }
  if (nlyr(env) != nrow(RDA$CCA$biplot) || any(! names(env) %in% rownames(RDA$CCA$biplot))) {
    stop("\n env must contain same variables as the ones used in RDA (see rownames(RDA$CCA$biplot))")
  }
  if (!is.null(env_mask)) {
    if (!inherits(env_mask, c("SpatRaster", "RasterLayer"))) {
      stop("\n env_mask must be a 'SpatRaster' or 'RasterLayer' object")
    } else if (inherits(env_mask, "RasterLayer")) {
      env_mask <- rast(env_mask)
    }
    if (nlyr(env_mask) > 1) {
      stop("\n env_mask must contain only one layer")
    }
  }
  if (!(method %in% c("loadings", "predict"))) {
    stop("\n method must be 'loadings' or 'predict'")
  }
  # length(scale_env) == nrow(RDA_biplot)
  # length(center_env) == nrow(RDA_biplot)
  # names(scale_env) == rownames(RDA_biplot)
  # names(center_env) == rownames(RDA_biplot)
  
  
  ## FUNCTION -----------------------------------------------------------------
  
  ## Get RDA informations -----------------------------------------------------
  RDA_biplot <- RDA$CCA$biplot
  var_names <- rownames(RDA_biplot)
  
  ## Transform environmental raster for prediction
  if (!is.null(env_mask)) { ## Mask with env_mask
    env <- try(mask(env, env_mask))
    if (inherits(env, "try-error")) {
      stop("\n env and env_mask must match in terms of extent, origin, crs")
    }
  }
  env_crs <- crs(env)
  env_df <- as.data.frame(env[[var_names]], xy = TRUE)
  env_xy <- env_df[, c("x", "y")]
  env_var <- env_df[, var_names]
  
  if (!(is.null(scale_env) & is.null(center_env))) { ## Standardize environmental variables
    # env_var <- as.data.frame(scale(env_var, center_env[var_names], scale_env[var_names])) ## MARCHE PAS
    env_var <- as.data.frame(scale(env_var, center_env, scale_env)) ## one value only
  }
  
  
  ## MAKE PREDICTIONS ---------------------------------------------------------
  if (method == "predict") {
    pred <- predict(RDA, env_var, type = "lc")
  }
  Proj <- foreach(i = 1:K) %do%
    {
      if (method == "loadings") { ## Predict pixels genetic component based on RDA axes
        tmp_df <- data.frame(env_xy, z = as.vector(apply(env_var, 1, function(x) sum(x * RDA_biplot[, i]))))
      } else if (method == "predict") { ## Predict with RDA model and linear combinations
        tmp_df <- data.frame(env_xy, z = as.vector(pred[, i]))
      }
      ras <- rast(tmp_df, type = "xyz", crs = crs(env))
      names(ras) <- paste0("RDA", i)
      return(ras)
    }
  Proj <- rast(Proj)
  
  
  ## Return predictions for each RDA axis
  return(Proj)
}

