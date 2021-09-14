# -*- coding: utf-8 -*-
"""
@author: Hjohan
"""

#import relevant packages
from sklearn.datasets import make_regression
from sklearn.linear_model import LinearRegression
import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
#from sklearn.metrics import r2_score

#generate the synthetic linear model
n_samples, n_features = 1000, 17
rng = np.random.RandomState(0)
x, y = make_regression(n_samples, n_features, noise=70, random_state=rng)
#plt.scatter(x,y)

#fit the data and get R2
lr = LinearRegression()
lr.fit(x, y)
R2=lr.score(x, y)
print(R2)

#get the model input and output data
dfx = pd.DataFrame(x)
dfy = pd.DataFrame(y)
dfx.to_excel('File path\File name.xlsx', index = False)
dfy.to_excel('File path\File name2.xlsx', index = False)
