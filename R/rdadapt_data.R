#' Maladaptation tutorial
#'
#' A \code{list} containing names, coordinates and genomic informations for 
#' \pkg{rdadapt} tutorial.
#' 
#'
#' @format A \code{list} object with 3 elements:
#' \describe{
#'   \item{pop}{a \code{vector} identifying population for 2036 individuals}
#'   \item{coord}{a \code{data.frame} containing coordinates for 2036 individuals}
#'   \item{Y}{a \code{data.frame} containing genomic information for 2036 
#'   individuals about 2333 genes}
#' }

"maladapt"

# ## Get population labels
# pop <- read.csv2("data/pop.csv", sep = ",", head = TRUE)[,-1]
# 
# ## Get individual coordinates
# coord <- read.csv2("data/position.csv", sep = ",", head = TRUE)[,-1]
# coord <- data.frame(apply(coord, 2, as.numeric))
# 
# ## Get genetic dataset
# Y <- read.csv2("data/genome.csv", sep = ",", head = TRUE, row.names = 1)
# 
# ##Gather verything in a list
# maladapt <- list(pop = pop, coord = coord, Y = Y)

# usethis::use_data(maladapt, overwrite = TRUE)
