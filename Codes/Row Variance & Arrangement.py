#!/usr/bin/env python
# coding: utf-8

# In[ ]:


library(matrixStats)


# In[ ]:


df$row_var = rowVars(as.matrix(df[,]))


# In[ ]:


df <-df[order(df$row_var),]

