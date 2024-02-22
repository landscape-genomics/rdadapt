###################################################################################################
##' @name rdadapt
##' @author Thibaut Capblancq
##' 
##' @title RDA based genome scan
##' 
##' @description This function allows the user to ...
##' 
##' @param rda blabla
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
##' 
##' @export
##' 
##'
###################################################################################################

rdadapt <- function(rda, K)
{
  ## Checks
  # inherits(rda, "rda")
  # "CCA" %in% names(rda) ## should not be necessary ?
  # "v" %in% names(rda$CCA) ## should not be necessary ?
  # K <= ncol(rda$CCA$v)
  
  zscores <- rda$CCA$v[, 1:as.numeric(K)]
  return(rdadapt.z(zscores))
}

rdadapt.z <- function(zscores)
{
  K <- ncol(zscores)
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action = na.omit, estim = "pairwiseGK")$dist
  lambda <- median(resmaha) / qchisq(0.5, df = K)
  reschi2test <- pchisq(resmaha / lambda, K, lower.tail = FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt <- qval$qvalues
  return(data.frame(p.values = reschi2test, q.values = q.values_rdadapt))
}
