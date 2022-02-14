#!/usr/bin/env python
# coding: utf-8

# In[ ]:


library(GEOquery)


# In[ ]:


my_id <- "GSE4843"
gse <- getGEO(my_id)


# In[ ]:


gse <-gse[[1]]


# In[ ]:


df_old<-as.data.frame(gse)


# In[ ]:


matrix_old=t(df_old)


# In[ ]:


matrix_new<-head(matrix_old,-26)


# In[ ]:


matrix_new_t=t(matrix_new)


# In[ ]:


df<-as.data.frame(matrix_new_t)


# In[ ]:


num<-as.data.frame(as.numeric(unlist(df)))


# In[ ]:


expression = exprs(gse)


# In[ ]:


df2 <- data.matrix(df)


# In[ ]:


"download GPL file and load it"


# In[ ]:


GPL<- read.delim(file.choose(), stringsAsFactor = FALSE)


# In[ ]:


df_favoured<-as.data.frame(matrix_new)


# In[ ]:


for ( row in 1:nrow(df_favoured)){
    row.names(df_favoured)[row] <-  sub("X", "", row.names(df_favoured)[row])
}


# In[ ]:


df_favoured<-lapply(df_favoured,as.numeric)


# In[ ]:


x3<-cbind(df_favoured,GPL)


# In[ ]:


x3<-as.data.frame(x3)


# In[ ]:


agg=aggregate(x3[,1:45],by=list(x3$Gene.Symbol),FUN=mean,na.rm=TRUE)


# In[ ]:


data<-as.data.frame(agg[-c(1),])


# In[ ]:


data<-data.matrix(data)

