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
##' @param env_mask (\emph{optional, default} \code{NULL}) \cr blabla
##' @param method  (\emph{default} \code{'loadings'}) \cr blabla
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