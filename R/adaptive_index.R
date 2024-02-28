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

setGeneric("adaptive_index", def = function(RDA, K, env, env_mask = NULL, method = "loadings") { standardGeneric("adaptive_index") })

##'
##' @rdname adaptive_index
##' @export
##'

setMethod('adaptive_index', signature(RDA = "rda", env = "SpatRaster"), function(RDA, K, env, env_mask = NULL, method = "loadings")
{
  ## CHECKS -------------------------------------------------------------------
  if (!inherits(RDA, "rda")) { stop("\n RDA must be a 'rda' object") }
  if (!("CCA" %in% names(RDA)) || !("biplot" %in% names(RDA$CCA))) {
    stop("\n RDA$CCA$biplot seems not to exist")
  }
  if (K > ncol(RDA$CCA$biplot)) {
    K <- ncol(RDA$CCA$biplot)
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
  
  
  ## FUNCTION -----------------------------------------------------------------
  
  ## Get RDA information -----------------------------------------------------
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
  
  
  ## Make predictions ---------------------------------------------------------
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
})


##'
##' @rdname adaptive_index
##' @export
##'

setMethod('adaptive_index', signature(RDA = "rda", env = "data.frame", env_mask = "missing"), function(RDA, K, env, method = "loadings")
{
  ## CHECKS -------------------------------------------------------------------
  if (!inherits(RDA, "rda")) { stop("\n RDA must be a 'rda' object") }
  if (!("CCA" %in% names(RDA)) || !("biplot" %in% names(RDA$CCA))) {
    stop("\n RDA$CCA$biplot seems not to exist")
  }
  if (K > ncol(RDA$CCA$biplot)) {
    K <- ncol(RDA$CCA$biplot)
    warning("\n Not enough RDA axis available, K is set to ncol(RDA$CCA$biplot)")
  }
  if (ncol(env) != nrow(RDA$CCA$biplot) || any(! colnames(env) %in% rownames(RDA$CCA$biplot))) {
    stop("\n env must contain same variables as the ones used in RDA (see rownames(RDA$CCA$biplot))")
  }
  if (!(method %in% c("loadings", "predict"))) {
    stop("\n method must be 'loadings' or 'predict'")
  }
  
  
  ## FUNCTION -----------------------------------------------------------------
  
  ## Get RDA information -----------------------------------------------------
  RDA_biplot <- RDA$CCA$biplot
  var_names <- rownames(RDA_biplot)
  env_var <- env[, var_names]
  
  ## Make predictions ---------------------------------------------------------
  if (method == "predict") {
    pred <- predict(RDA, env_var, type = "lc")
  }
  AI <- foreach(i = 1:K) %do%
    {
      if (method == "loadings") { ## Predict genetic component based on RDA axes
        tmp_vec <- as.vector(apply(env_var, 1, function(x) sum(x * RDA_biplot[, i])))
      } else if (method == "predict") { ## Predict with RDA model and linear combinations
        tmp_vec <- as.vector(pred[, i])
      }
      return(tmp_vec)
    }
  names(AI) <- paste0("RDA", 1:K)
  AI <- as.data.frame(do.call(cbind, AI))
  
  ## Return predictions for each RDA axis
  return(AI)
})