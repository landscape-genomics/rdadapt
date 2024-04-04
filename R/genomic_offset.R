###################################################################################################
##' @name genomic_offset
##' @author Thibaut Capblancq
##' 
##' @title Predict temporal or geographic genomic offset
##' 
##' @description The \code{genomic_offset} function allows the user to predict genomic offset 
##' from a RDA model. The \code{genomic_offset} function can estimate both temporal or spatial 
##' genomic offset and accommodates raster data or discrete populations.
##' 
##' @param RDA a \code{RDA} model from which to extract loci and environmental variable scores
##' @param K an \code{integer} specifying the number of RDA axes to use for the projection
##' @param env_pres a \code{RasterStack} object or a \code{data.frame} with the environmental 
##' conditions in the present
##' @param env_fut a \code{RasterStack} object or a \code{data.frame} with the environmental 
##' conditions in the future
##' @param env_mask (\emph{optional, default} \code{NULL}) \cr a \code{Raster} object to limit 
##' the projection to a specific area
##' @param method  (\emph{default} \code{'loadings'}) \cr a \code{character} defining whether 
##' the function is to use weighted averages (scaling type 1, \code{loadings}) or linear 
##' combinations (scaling type 2, \code{predict}) of the projected environmental variables 
##' to predict site scores (i.e., adaptive index)
##' 
##' 
##' @return  
##' 
##' A \code{list} containing :
##' \itemize{
##'   \item \code{genomic_offset} : a \code{RasterStack} or a \code{data.frame} containing 
##'   the genomic offset predictions for the \code{K} first RDA axes, as well as the overall 
##'   genomic offset prediction
##'   \item \code{weights} : the weights associated with each RDA axis used for the predictions
##' }
##' 
##' 
##' @details
##' 
##' This RDA-based method to predict \emph{genomic offset} is relatively simple. 
##' RDA is first used to predict the \emph{optimal adaptive genetic composition} for each 
##' environmental pixel under consideration (see \code{\link{adaptive_index}} function), 
##' using both current and future environmental conditions. 
##' The \emph{euclidean distance} between these two predictions in the RDA space provides an 
##' estimate of the change in genetic composition that would be required to track climate change.
##' 
##' 
##' @keywords RDA, genomic offset, adaptive index, projection
## 
##' @seealso \code{\link{adaptive_index}}
## 
## @examples
##' 
##' 
##' @importFrom stats dist
##' @importFrom raster rasterToPoints rasterFromXYZ mask
##' @include rda_class.R
##'  
##' @export
##' 
##'
###################################################################################################

setGeneric("genomic_offset", def = function(RDA, K, env_pres, env_fut, env_mask = NULL, env_gar, method = "loadings") { standardGeneric("genomic_offset") })

##'
##' @rdname genomic_offset
##' @export
##'

setMethod('genomic_offset', signature(RDA = "rda", env_pres = "SpatRaster", env_gar = "missing"), function(RDA, K, env_pres, env_fut, env_mask = NULL, method = "loadings")
{
  ## CHECKS -------------------------------------------------------------------
  if (!inherits(RDA, "rda")) { stop("\n RDA must be a 'rda' object") }
  if (!("CCA" %in% names(RDA)) || !("eig" %in% names(RDA$CCA))) {
    stop("\n RDA$CCA$eig seems not to exist")
  }
  if (K > length(RDA$CCA$eig)) {
    K <- length(RDA$CCA$eig)
    warning("\n Not enough RDA axis available, K is set to length(RDA$CCA$eig)")
  }
  
  ## FUNCTION -----------------------------------------------------------------
  
  ## Make predictions
  AI_pres <- adaptive_index(RDA = RDA
                            , K = K
                            , env = env_pres
                            , env_mask = env_mask
                            , method = method)
  AI_fut <- adaptive_index(RDA = RDA
                           , K = K
                           , env = env_fut
                           , env_mask = env_mask
                           , method = method)
  
  ## Single axis genetic offset -----------------------------------------------
  offset <- foreach(i = 1:K) %do%
    {
      ras <- abs(AI_pres[[i]] - AI_fut[[i]])
      names(ras) <- paste0("RDA", i)
      return(ras)
    }
  offset <- rast(offset)
  
  ## Weight current and future adaptive indices based on eigen values of associated axes
  weights <- RDA$CCA$eig / sum(RDA$CCA$eig)
  weights <- weights[1:K]
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
})


##'
##' @rdname genomic_offset
##' @export
##'

setMethod('genomic_offset', signature(RDA = "rda", env_pres = "data.frame", env_mask = "missing", env_gar = "missing"), function(RDA, K, env_pres, env_fut, method = "loadings")
{
  ## CHECKS -------------------------------------------------------------------
  if (!inherits(RDA, "rda")) { stop("\n RDA must be a 'rda' object") }
  if (!("CCA" %in% names(RDA)) || !("eig" %in% names(RDA$CCA))) {
    stop("\n RDA$CCA$eig seems not to exist")
  }
  if (K > length(RDA$CCA$eig)) {
    K <- length(RDA$CCA$eig)
    warning("\n Not enough RDA axis available, K is set to length(RDA$CCA$eig)")
  }
  
  ## FUNCTION -----------------------------------------------------------------
  
  ## Make predictions
  AI_pres <- adaptive_index(RDA = RDA
                            , K = K
                            , env = env_pres
                            , method = method)
  AI_fut <- adaptive_index(RDA = RDA
                           , K = K
                           , env = env_fut
                           , method = method)
  
  ## Single axis genetic offset -----------------------------------------------
  offset <- foreach(i = 1:K) %do%
    {
      vec <- abs(AI_pres[,i] - AI_fut[,i])
      return(vec)
    }
  offset <- as.data.frame(do.call(cbind, offset))
  
  ## Weight offset based on eigen values of associated axes
  weights <- RDA$CCA$eig / sum(RDA$CCA$eig)
  weights <- weights[1:K]
  offset <- data.frame(do.call(cbind, lapply(1:ncol(offset), function(x) offset[,x]*weights[x])))
  colnames(offset) <- paste0("RDA", 1:K)
  
  ## Predict a global genetic offset
  offset$Global <- apply(offset, 1, function(x) sqrt(sum(x*x)))
  
  return(list(Proj_pres = AI_pres,
              Proj_fut = AI_fut,
              genomic_offset = offset,
              weights = weights))
})


##'
##' @rdname genomic_offset
##' @export
##'

setMethod('genomic_offset', signature(RDA = "rda", env_pres = "data.frame", env_mask = "missing", env_fut = "missing"), function(RDA, K, env_pres, env_gar, method = "loadings")
{
  ## CHECKS -------------------------------------------------------------------
  if (!inherits(RDA, "rda")) { stop("\n RDA must be a 'rda' object") }
  if (!("CCA" %in% names(RDA)) || !("eig" %in% names(RDA$CCA))) {
    stop("\n RDA$CCA$eig seems not to exist")
  }
  if (K > length(RDA$CCA$eig)) {
    K <- length(RDA$CCA$eig)
    warning("\n Not enough RDA axis available, K is set to length(RDA$CCA$eig)")
  }
  
  ## FUNCTION -----------------------------------------------------------------
  
  ## Make predictions
  AI_pres <- adaptive_index(RDA = RDA
                            , K = K
                            , env = env_pres
                            , method = method)
  
  AI_gar <- adaptive_index(RDA = RDA
                           , K = K
                           , env = env_gar
                           , method = method)
  
  ## Single axis genetic offset -----------------------------------------------
  offset <- foreach(i = 1:K) %do%
    {
      vec <- abs(AI_pres[,i] - AI_gar[,i])
      return(vec)
    }
  offset <- as.data.frame(do.call(cbind, offset))
  
  ## Weight offset based on eigen values of associated axes
  weights <- RDA$CCA$eig / sum(RDA$CCA$eig)
  weights <- weights[1:K]
  offset <- data.frame(do.call(cbind, lapply(1:ncol(offset), function(x) offset[,x]*weights[x])))
  colnames(offset) <- paste0("RDA", 1:K)
  
  ## Predict a global genetic offset
  offset$Global <- apply(offset, 1, function(x) sqrt(sum(x*x)))
  
  return(list(genomic_offset = offset,
              weights = weights))
})