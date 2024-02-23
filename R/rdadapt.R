###################################################################################################
##' @name rdadapt
##' @author Thibaut Capblancq
##' 
##' @title RDA based genome scan
##' 
##' @description This function allows the user to ...
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
##' @importFrom stats median qchisq pchisq
##' @importFrom robust covRob
##' @importFrom qvalue qvalue
##' @importClassesFrom vegan rda
##' 
##' @export
##' 
##'
###################################################################################################


setGeneric("rdadapt", def = function(RDA, K, zscores) { standardGeneric( "rdadapt") })

##'
##' @rdname rdadapt
##' @export
##'

setMethod('rdadapt', signature(RDA = "rda"), function(RDA, K)
{
  ## CHECKS -------------------------------------------------------------------
  if (!inherits(RDA, "rda")) { stop("\n RDA must be a 'rda' object") }
  if (!("CCA" %in% names(RDA)) || !("v" %in% names(RDA$CCA))) {
    stop("\n RDA$CCA$v seems not to exist")
  }
  if (K > ncol(RDA$CCA$v)) {
    K = ncol(RDA$CCA$v)
    warning("\n Not enough RDA axis available, K is set to ncol(RDA$CCA$v)")
  }
  
  ## FUNCTION -----------------------------------------------------------------
  zscores <- RDA$CCA$v[, 1:as.numeric(K)]
  return(rdadapt(zscores = zscores))
})

##'
##' @rdname rdadapt
##' @export
##'

setMethod('rdadapt', signature(zscores = "data.frame", RDA = "missing", K = "missing"), function(zscores)
{
  zscores <- as.matrix(zscores)
  return(rdadapt(zscores = zscores))
})

##'
##' @rdname rdadapt
##' @export
##'

setMethod('rdadapt', signature(zscores = "matrix", RDA = "missing", K = "missing"), function(zscores)
{
  K <- ncol(zscores)
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action = na.omit, estim = "pairwiseGK")$dist
  lambda <- median(resmaha) / qchisq(0.5, df = K)
  reschi2test <- pchisq(resmaha / lambda, K, lower.tail = FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt <- qval$qvalues
  return(data.frame(p.values = reschi2test, q.values = q.values_rdadapt))
})

