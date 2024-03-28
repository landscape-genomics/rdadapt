##' 
##' @export
##' 

## Class not exported by vegan, so initialized here
setClass("rda", slots = c("colsum", "tot.chi", "Ybar"
                          , "method", "call"
                          , "pCCA", "CCA", "CA", "inertia"
                          , "regularization", "terms", "terminfo"))
