#!/usr/bin/env python
# coding: utf-8

# In[ ]:


library(tidyverse)
library(ggplot2)
library(dplyr)
library(cluster)
library(factoextra)
library(devtools)


# In[ ]:


df<-column_to_rownames(df, var = "Group.1")


# In[ ]:


pca = prcomp(df, center = FALSE, scale = FALSE)
summary(pca)


# In[ ]:


transform = as.data.frame(-pca$x[,1:2])


# In[ ]:


fviz_nbclust(transform, kmeans, method = 'wss')


# In[ ]:


fviz_nbclust(transform, kmeans, method = 'silhouette')
fviz_nbclust(transform, kmeans, method = 'gap_stat')


# In[ ]:


k = 2
kmeans = kmeans(transform, centers = k, nstart = 50)
fviz_cluster(kmeans, data = transform)


# In[ ]:
kmeans_tdata_clus<-as.data.frame(kmeans.cluster)

tdata_cluster1<-as.data.frame(kmeans_tdata_cluster)
tdata_cluster2<-as.data.frame(kmeans_tdata_cluster)


# In[ ]:


tdata_cluster1<-tdata_cluster1[-c(1,4,5,6,9,10,11,12,13,14,15,18,19,21,24,26,27,28,29,32,33,34,35,37,38,39,40,44),]
tdata_cluster1<-tdata_cluster1[,-c(3001)]


# In[ ]:


tdata_cluster2<-tdata_cluster2[c(1,4,5,6,9,10,11,12,13,14,15,18,19,21,24,26,27,28,29,32,33,34,35,37,38,39,40,44),]
tdata_cluster2<-tdata_cluster2[,-c(3001)]



MEs1<-cbind(kmeans_tdata_clus,MEs)
MEs1.m <- melt(MEs1, id.var = "cluster")
ggplot(data = MEs1.m, aes(x=variable, y=value, color=cluster)) + geom_boxplot(aes(fill=cluster))


boxplot(MEs1[,-1], boxfill = NA, border = NA) #invisible boxes - only axes and plot area
boxplot(MEs1[MEs1$cluster=="1", -1], xaxt = "n", add = TRUE, boxfill="red",
        boxwex=0.25, at = 1:ncol(MEs1[,-1]) - 0.15) #shift these left by -0.15
boxplot(MEs1[MEs1$cluster=="2", -1], xaxt = "n", add = TRUE, boxfill="blue",
        boxwex=0.25, at = 1:ncol(MEs1[,-1]) + 0.15)
