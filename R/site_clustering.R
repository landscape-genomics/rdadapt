#### Function to find groups of covarying sites/indnividuals/populations on the RDA space


# Required libraries
library(factoextra)
library(plyr)
library(geometry)

# Function
loci_clustering <- function(RDA, K){

# RDA1 and RDA2 scores for all the pixels of the range
TAB <- data.frame(x = rasterToPoints(res_RDA_proj_current$RDA1)[,1], y = rasterToPoints(res_RDA_proj_current$RDA1)[,2], RDA1 = rasterToPoints(res_RDA_proj_current$RDA1)[,3], RDA2 = rasterToPoints(res_RDA_proj_current$RDA2)[,3])

# Hierarchical K-means clustering for 10,000 random pixels
TAB_sample <- TAB[sample(1:nrow(TAB), 10000),3:4]
clust <- list()
for(i in 1:15){
  clust[[i]] <- hkmeans(x = TAB_sample, k = i, hc.metric = "euclidean", hc.method = "ward.D2", iter.max = 10, km.algorithm = "Hartigan-Wong")
}

# Select the best number of clusters
within_cust_var <- lapply(clust, function(x) x$withinss/x$totss)
plot(unlist(lapply(within_cust_var, mean)))

# Extrapolate the clustering to all the pixels with 4 clusters
df <- data.frame(x = TAB_sample$RDA1, y = TAB_sample$RDA2, cluster = clust[[4]]$cluster) 
find_hull <- function(df) df[chull(df$x, df$y),]
hulls <- ddply(df, "cluster", find_hull) # Convex hulls.
cluster <- rep(NA, nrow(TAB))
for(i in 1:4){
  inhull <- inhulln(convhulln(hulls[hulls$cluster==i,1:2]), as.matrix(TAB[,3:4])) # Intersection convex hulls and points
  cluster[inhull] <- i
}
clusterNA <- cluster[which(is.na(cluster))]
for(i in 1:length(clusterNA)){
  clusterNA[i] <- which.min(pointDistance(TAB[which(is.na(cluster))[i], 3:4], clust[[4]]$centers, lonlat = F)) # Fill the NA using the distance to the hull centers
}
cluster[which(is.na(cluster))] <- clusterNA

}