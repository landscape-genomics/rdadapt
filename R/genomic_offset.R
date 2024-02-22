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
  
  ## Get RDA informations -----------------------------------------------------
  rda_biplot <- rda$CCA$biplot
  var_names <- row.names(rda_biplot)
  weights <- rda$CCA$eig / sum(rda$CCA$eig) ## Weights based on axis eigen values
  
  ## Transform environmental raster for prediction
  if (!is.null(range)) { ## Mask with range
    env_pres <- mask(env_pres, range)
    env_fut <- mask(env_fut, range)
  }
  env_crs <- crs(env_pres)
  env_df_pres <- as.data.frame(env_pres[[var_names]], xy = TRUE)
  env_df_fut <- as.data.frame(env_fut[[var_names]], xy = TRUE)
  env_xy <- env_df_pres[, c("x", "y")]
  env_var_pres <- env_df_pres[, var_names]
  env_var_fut <- env_df_fut[, var_names]
  
  if (!(is.null(scale_env) & is.null(center_env))) { ## Standardize environmental variables
    # env_var <- as.data.frame(scale(env_var, center_env[var_names], scale_env[var_names])) ## MARCHE PAS
    env_var_pres <- as.data.frame(scale(env_var_pres, center_env, scale_env)) ## one value only
    env_var_fut <- as.data.frame(scale(env_var_fut, center_env, scale_env)) ## one value only
  }
  
  
  ## MAKE PREDICTIONS ---------------------------------------------------------
  if (method == "predict") {
    pred_pres <- predict(rda, env_var_pres, type = "lc")
    pred_fut <- predict(rda, env_var_fut, type = "lc")
  }
  RES <- foreach(i = 1:K) %do%
    {
      ## Current and future predictions
      Proj <- foreach(sce = c("pres", "fut"), env = list(env_var_pres, env_var_fut), pred = list(pred_pres, pred_fut)) %do%
        {
          if (method == "loadings") {
            ras <- env_pres[[1]]
            ras[!is.na(ras)] <- as.vector(apply(env, 1, function(x) sum(x * rda_biplot[, i])))
          } else {
            tmp_df <- data.frame(env_xy, Z = as.vector(pred[, i]))
            ras <- rast(tmp_df, type = "xyz", crs = env_crs)
          }
          names(ras) <- paste0("RDA", i, "_", sce)
          return(ras)
        }
      Proj <- rast(Proj)
      
      ## Single axis genetic offset 
      Proj_offset <- abs(Proj[[1]] - Proj[[2]])
      names(Proj_offset) <- paste0("RDA", i)
      
      return(list(proj = Proj, offset = Proj_offset))
    }
  names(RES) <- paste0("RDA", 1:K)
  
  
  ## Weight current and future adaptive indices based on eigen values of associated axes
  proj_df <- do.call(cbind, lapply(1:K, function(x) as.data.frame(RES[[x]]$proj, xy = FALSE)))
  for (i in 1:K) {
    ind_i <- grep(names(weights)[i], colnames(proj_df))
    for (ii in ind_i) {
      proj_df[, ii] <- proj_df[, ii] * weights[i]
    }
  }
  
  ## Predict a global genetic offset
  offset_global <- RES[[1]]$offset
  offset_global[!is.na(offset_global)] <- sapply(1:nrow(proj_df), function(x) {
    tmp_mat <- rbind(unlist(proj_df[x, grep("pres", colnames(proj_df))])
                     , unlist(proj_df[x, grep("fut", colnames(proj_df))]))
    return(dist(tmp_mat, method = "euclidean"))
  })
  names(offset_global) <- "Global_offset"
  
  # Return projections for current and future climates for each RDA axis, prediction of genetic offset for each RDA axis and a global genetic offset 
  # return(list(Proj_pres = Proj_pres, 
  #             Proj_fut = Proj_fut, 
  #             Proj_offset = Proj_offset, 
  #             Proj_offset_global = Proj_offset_global, weights = weights[1:K]))
  return(offset_global)
}
