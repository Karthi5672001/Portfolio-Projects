#!/usr/bin/env python
# coding: utf-8

# In[41]:


import pandas as pd
import numpy as np
import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib
plt.style.use('ggplot')
from matplotlib.pyplot import figure

get_ipython().run_line_magic('matplotlib', 'inline')
matplotlib.rcParams['figure.figsize'] = (12,8)

pd.options.mode.chained_assignment = None



# We read in the data

df = pd.read_csv(r'C:\Users\karth\Desktop\Data Science\Corana portfolio\movies.csv')

# We look at the data

df


# In[4]:


# We check if any missing data

for col in df.columns:
    pct_missing = np.mean(df[col].isnull())
    print('{} - {}%'.format(col, round(pct_missing*100)))
    


# In[6]:


# Display data types for our columns

print(df.dtypes)


# In[7]:


#We examine a boxplot to see the distribution and check for outliers
df.boxplot(column=['gross'])


# In[8]:


#Cleaning up data by getting rid of duplicates
df.drop_duplicates()


# In[43]:


#Delete released column,as it is unnecessary for analysis
del df['released']


# In[44]:


df


# In[45]:


# Order our Data by Score

df.sort_values(by=['score'], inplace=False, ascending=False)


# In[46]:


#Seeing correlation between score and budget - We find a weak correlation if any.
sns.regplot(x="score", y="budget", data=df)


# In[57]:


#Seeing correlation between budget and gross - We find some correlation.
sns.regplot(x="budget", y="gross", data=df)


# In[58]:


#Examining correlation
df.corr()


# In[49]:


# Correlation Matrix for pearson

df.corr(method ='pearson')


# In[50]:


# Correlation Matrix for spearman

df.corr(method ='spearman')


# In[51]:


# Correlation Matrix for kendall

df.corr(method ='kendall')


# In[67]:


#Pearson correlation. Dark means low correlation, bright is high correlation.
 
correlation_matrix = df.corr()

sns.heatmap(correlation_matrix, annot = True)

plt.title("Correlation matrix for Numeric Values")

plt.xlabel("Movie features")

plt.ylabel("Movie features")

plt.show()


# In[66]:


#Spearman correlation. Dark means low correlation, bright is high correlation.
correlation_matrix = df.corr(method ='spearman')

sns.heatmap(correlation_matrix, annot = True)

plt.title("Correlation matrix for Numeric Values")

plt.xlabel("Movie features")

plt.ylabel("Movie features")

plt.show()


# In[65]:


#Kendall correlation. Dark means low correlation, bright is high correlation.
correlation_matrix = df.corr(method ='kendall')

sns.heatmap(correlation_matrix, annot = True)

plt.title("Correlation matrix for Numeric Values")

plt.xlabel("Movie features")

plt.ylabel("Movie features")

plt.show()


# In[62]:


# Using factorize - this assigns a random numeric value for each unique categorical value in pearson

df.apply(lambda x: x.factorize()[0]).corr(method='pearson')


# In[73]:


# Using factorize - this assigns a random numeric value for each unique categorical value in spearman

df.apply(lambda x: x.factorize()[0]).corr(method='spearman')


# In[74]:


# Using factorize - this assigns a random numeric value for each unique categorical value in kendall

df.apply(lambda x: x.factorize()[0]).corr(method='kendall')


# In[70]:


#Lambda Correlation matrix for Movies(Pearson)
correlation_matrix = df.apply(lambda x: x.factorize()[0]).corr(method='pearson')

sns.heatmap(correlation_matrix, annot = True)

plt.title("Correlation matrix for Movies")

plt.xlabel("Movie features")

plt.ylabel("Movie features")

plt.show()


# In[71]:


#Lambda Correlation matrix for Movies(Spearman)
correlation_matrix = df.apply(lambda x: x.factorize()[0]).corr(method='spearman')

sns.heatmap(correlation_matrix, annot = True)

plt.title("Correlation matrix for Movies")

plt.xlabel("Movie features")

plt.ylabel("Movie features")

plt.show()


# In[72]:


#Lambda Correlation matrix for Movies(Kendall)
correlation_matrix = df.apply(lambda x: x.factorize()[0]).corr(method='kendall')

sns.heatmap(correlation_matrix, annot = True)

plt.title("Correlation matrix for Movies")

plt.xlabel("Movie features")

plt.ylabel("Movie features")

plt.show()


# In[75]:


correlation_mat = df.apply(lambda x: x.factorize()[0]).corr()
corr_pairs = correlation_mat.unstack()
print(corr_pairs)


# In[76]:


sorted_pairs = corr_pairs.sort_values(kind="quicksort")
print(sorted_pairs)


# In[77]:


# We can now relationships that have a high correlation (> 0.5)

strong_pairs = sorted_pairs[abs(sorted_pairs) > 0.5]

print(strong_pairs)


# In[79]:


# Looking at the top 20 compaies by gross revenue

CompanyGrossSum = df.groupby('company')[["gross"]].sum()

CompanyGrossSumSorted = CompanyGrossSum.sort_values('gross', ascending = False)[:20]

CompanyGrossSumSorted = CompanyGrossSumSorted['gross'].astype('int64') 

CompanyGrossSumSorted


# In[ ]:




