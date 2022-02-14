

# In[ ]:


library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
allowWGCNAThreads()
library(matrixStats)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(cluster)
library(factoextra)
library(devtools)
library(ggpubr)
library(rstatix)
library(ggprism)
library(patchwork)
library(magrittr)

tdf<-df
df<-as.data.frame(t(tdf))
df<-column_to_rownames(df, var = "Group.1")


powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(df, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# In[ ]:


softPower = 4;
adjacency = adjacency(df, power = softPower);


# In[ ]:



TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


# In[ ]:



geneTree = hclust(as.dist(dissTOM), method = "average");

sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
    labels = FALSE, hang = 0.04);


# In[ ]:


minModuleSize = 30;

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
              deepSplit = 2, pamRespectsDendro = FALSE,
              minClusterSize = minModuleSize);
table(dynamicMods)


# In[ ]:


dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                  dendroLabels = FALSE, hang = 0.03,
                  addGuide = TRUE, guideHang = 0.05,
                  main = "Gene dendrogram and module colors")


# In[ ]:



MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

MEDiss = 1-cor(MEs);

METree = hclust(as.dist(MEDiss), method = "average");

sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")
MEDissThres = 0.25

abline(h=MEDissThres, col = "red")

merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

mergedColors = merge$colors;

mergedMEs = merge$newMEs;


# In[ ]:


sizeGrWindow(12, 9)

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                  c("Dynamic Tree Cut", "Merged dynamic"),
                  dendroLabels = FALSE, hang = 0.03,
                  addGuide = TRUE, guideHang = 0.05)


# In[ ]:


moduleColors = mergedColors

colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;





tdf$row_var = rowVars(as.matrix(tdf[,]))


tdf <-tdf[order(tdf$row_var),]  #arranging in ascending order

tdf<-tdf[-c(1:20520),] #deleting other rows

# In[ ]:
tdf<-tdf[,-c(46)]

# In[ ]:
tdf<-as.data.frame(t(tdf))
k = 2
kmeans_cluster = kmeans(tdf, centers = k, nstart = 50)
fviz_cluster(kmeans, data = tdf)


# In[ ]:
kmeans_cluster_number<-as.data.frame(kmeans_cluster$cluster)

colnames(kmeans_cluster_number)<-c('number')
data<-cbind(kmeans_cluster_number,tdf)

kmeans_cluster1<-as.data.frame(data)
kmeans_cluster2<-as.data.frame(data)

# In[ ]:
kmeans_cluster1<-data[-c(1,5,6,9,10,11,12,13,14,15,18,19,21,24,26,27,28,29,32,33,34,35,37,38,39,40,44),]
kmeans_cluster1<-kmeans_cluster1[,-c(1)]

#clubbed the  cluster  one data points in one separate file

kmeans_cluster2<-data[c(1,4,6,9,10,11,12,13,14,15,18,19,21,24,26,27,28,29,32,33,34,35,37,38,39,40,44),]
kmeans_cluster2<-kmeans_cluster2[,-c(1)]

#clubbed the  cluster two data points in another separate file

MEs_cluster<-cbind(kmeans_cluster_number,MEs)
MEs_cluster.m <- melt(MEs_cluster, id.var = "number")
ggplot(data = MEs_cluster.m, aes(x=variable, y=value, color=cluster)) + geom_boxplot(aes(fill=cluster))


boxplot(MEs_cluster[,-1], boxfill = NA, border = NA) #invisible boxes - only axes and plot area
boxplot(MEs_cluster[MEs_cluster$number=="1", -1], xaxt = "n", add = TRUE, boxfill="red",
        boxwex=0.25, at = 1:ncol(MEs_cluster[,-1]) - 0.15) #shift these left by -0.15
boxplot(MEs_cluster[MEs_cluster$number=="2", -1], xaxt = "n", add = TRUE, boxfill="blue",
        boxwex=0.25, at = 1:ncol(MEs_cluster[,-1]) + 0.15)
