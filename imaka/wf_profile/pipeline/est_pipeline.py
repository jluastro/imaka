##############
## Eden McEwen
## August 6th 2020
# rebuilding the estimator file

import numpy as np
import pandas as pd

import math
import time
import os
import sys

# Using Clustering
import hdbscan
import seaborn as sns

# Self-witten code
from code.Estimator import *

#### Change this
#### Output dir
out_dir = "/home/emcewen/out_est/"


### ITERATORS TO WORK THROUGH A DF ###

def df_iterator(df, xcor=False, sub_len=5, sig=3, thresh=3):
    """ 
    Takes in a df and iterates Estimator on each row
        input: data frame, xcor option, post subtraction length, sigma, threshold
        output: dataframe with estimates added
    """
    frames = []
    for index, row in df.iterrows():
        row_est = Estimator(row)
        est = row_est.estimate(xcor, sub_len, sig, thresh)
        frames.append(est)
    return pd.concat(frames, ignore_index=True)

def df_iterator_mult(df, sl_list, sig_list, thresh_list, c_list):
    """ 
    Takes in a df and iterates Estimator on each row, with variable options
        input: data frame, xcor options, post subtraction length options, sigma options, threshold options
        output: dataframe with estimates added
    """
    frames = [] 
    for index, row in df.iterrows():
        row_est = Estimator(row)        
        #TODO: figure out how to test multiple params
        for i in range(len(sl_list)):
            xc = False
            sl = sl_list[i]
            row_est.gen_wfs(xcor=xc, sub_len=sl)
            for j in range(len(sig_list)):
                for k in range(len(thresh_list)):
                    sg = sig_list[j]
                    thr = thresh_list[k]
                    row_est.gen_peaks(sig=sg, thresh=thr)
                    for m in range(len(c_list)):
                        cm = c_list[m]
                        row_est.gen_cluster(c_size=cm)
                        frames.append(row_est.get_df())
    return pd.concat(frames, ignore_index=True)
