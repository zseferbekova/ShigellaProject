## Load dependencies
# install.packages("ggrasp")
library(ggrasp)

## Import a nwk tree
tree <- ggrasp.load(" Results/Tree_midpoint.nwk")

## Cluster the tree into 7 clusters
tree_clustered <- ggrasp.cluster(tree, num.clusters = 7)

## See the results 
plot(tree_clustered)

## Write the results to an iTOL annotation file
> ggrasp.write(tree_clustered, type="itol", "Clusters.txt")