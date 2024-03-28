###################################################################################################
##' @name rdadapt
##' @author Thibaut Capblancq
##' 
##' @title Detect outlier loci based on their contribution to a RDA model
##' 
##' @description The \code{rdadapt} function performs redundancy analysis and computes p-values to
##' test for outliers based on loci extremeness along a distribution of Mahalanobis distances 
##' estimated between each locus and the center of the RDA space using a certain number of axes (K). 
##' \code{rdadapt} accommodates individual genotypes or allele frequencies.
##' 
##' @param RDA a RDA model from which to extract loci loadings
##' @param K an integer specifying the number of RDA axes to use for the detection
##' 
##' @return  
##' 
##' A \code{data.frame} containing :
##' \itemize{
##'   \item \code{p.values} : the p-value associated with each locus
##'   \item \code{q.values} : the q-value associated with each locus, estimated using the \code{q.values} package and allowing to use a FDR (False Discovery Rate) approach instead of a p-value threshold to identify outliers.
##' }
##' 
##' 
##' @details
##' 
##' First, Mahalanobis distances are computed for all genetic marker using a robust estimate of both mean and covariance matrix between the \code{K} RDA vectors of loadings.
##' Then, to compute p-values, Mahalanobis distances are divided by a genomic inflation factor (\code{gif}), giving a scaled statistic that should follow a chi-squared distribution with \code{K} degrees of freedom. 
##' 
## @keywords 
## 
## @seealso 
## 
## @examples
##' 
##' 
##' @importFrom stats median qchisq pchisq
##' @importFrom robust covRob
##' @importFrom qvalue qvalue
##' @include rda_class.R
##' 
##' @export
##' 
##'
###################################################################################################


setGeneric("rdadapt", def = function(RDA, K, scores) { standardGeneric( "rdadapt") })

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
  scores <- RDA$CCA$v[, 1:as.numeric(K)]
  return(rdadapt(scores = scores))
})

##'
##' @rdname rdadapt
##' @export
##'

setMethod('rdadapt', signature(scores = "data.frame", RDA = "missing", K = "missing"), function(scores)
{
  scores <- as.matrix(scores)
  return(rdadapt(scores = scores))
})

##'
##' @rdname rdadapt
##' @export
##'

setMethod('rdadapt', signature(scores = "matrix", RDA = "missing", K = "missing"), function(scores)
{
  K <- ncol(scores)
  resscale <- apply(scores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action = na.omit, estim = "pairwiseGK")$dist
  lambda <- median(resmaha) / qchisq(0.5, df = K)
  reschi2test <- pchisq(resmaha / lambda, K, lower.tail = FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt <- qval$qvalues
  return(data.frame(p.values = reschi2test, q.values = q.values_rdadapt))
})

