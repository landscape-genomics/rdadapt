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
# ## Get the current and future climatic variables
# xx.current <- read.table("data/var_current.csv", sep = ",", head = TRUE, row.names = 1)
# xx.future <- read.table("data/var_futur.csv", sep = ",", head = TRUE, row.names = 1)
# 
# ## Add extra variables (correlated)
# VAR3.current <- 2 * xx.current[, 1] - xx.current[, 2] + rnorm(nrow(xx.current), sd = 0.07)
# VAR4.current <- xx.current[, 1] + xx.current[, 2] + rnorm(nrow(xx.current), sd = 0.08)
# 
# VAR3.future <- 2 * xx.future[, 1] - xx.future[, 2] + rnorm(nrow(xx.future), sd = 0.07)
# VAR4.future <- xx.future[, 1] + xx.future[, 2] + rnorm(nrow(xx.future), sd = 0.08)
# 
# xx.current <- data.frame(xx.current, VAR3 = VAR3.current, VAR4 = VAR4.current)
# xx.future <- data.frame(xx.future, VAR3 = VAR3.future, VAR4 = VAR4.future)
# 
# ## Get adaptive loci for env variable 1 and env variable 2
# causal1 <- as.numeric(read.csv2("data/mutationm2.csv", sep = ",", header = TRUE)[,-1])
# causal2 <- as.numeric(read.csv2("data/mutationm3.csv", sep = ",", header = TRUE)[,-1])
# 
# ## Verifying that prediction are matching with a decline in fitness on the landscape
# fitness <- read.table("data/fitness_futur.csv", sep = ",", header = TRUE)
# 
# 
# ## Gather verything in a list
# maladapt <- list(pop = pop, coord = coord, Y = Y
#                  , xx.current = xx.current, xx.future = xx.future
#                  , causal1 = causal1, causal2 = causal2
#                  , fitness = fitness)

# usethis::use_data(maladapt, overwrite = TRUE)
