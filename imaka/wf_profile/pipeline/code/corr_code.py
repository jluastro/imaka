# Eden McEwen
# June 2 2020
# This file contains code for WFS cross correlation
# Helper functions are developed for the imaka data file format

import math
import imageio
import pandas as pd
import numpy as np
from scipy.io import readsav
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import matplotlib.animation as animation
from astropy.table import Table
from astropy.stats import sigma_clip
from astropy.io import fits

mask_8_8_center = [[0,0,1,1,1,1,0,0],
           [0,1,1,1,1,1,1,0],
           [1,1,1,1,1,1,1,1],
           [1,1,1,0,0,1,1,1],
           [1,1,1,0,0,1,1,1],
           [1,1,1,1,1,1,1,1],
           [0,1,1,1,1,1,1,0],
           [0,0,1,1,1,1,0,0]]

# Multiplies on data matrix
# TODO: Make more robust to the wrong size matrix
# TODO: Make more options for inputs?
def mask_data(mat):
    return mat*mask_8_8_center

# Takes out the time averaged slopes from data
def sub_data_avg(mat):
    # incoming data expected as: (t, x, y)
    return mat - np.average(mat, 0)

# Takes out the time averaged slopes from data
def sub_data_t_avg(mat):
    # incoming data expected as: (t, x, y)
    t_list = np.average(mat, axis=(1,2))
    t_list = np.reshape(t_list, (len(t_list), 1, 1))
    return mat - t_list


# Two Matrix cross Correlation
def mat_c_corr(mat_1, mat_2, mask):
    #assume mask is the same shape as mat1 and mat2
    t1, x1, y1 = mat_1.shape
    t2, x2, y2 = mat_2.shape
    # shift based on mat_1 only....
    di_range = range(-x1+1, x1)
    dj_range = range(-y1+1, y1)
    A = np.zeros((len(di_range), len(dj_range)))

    for di in di_range:
        for dj in dj_range:
            # tmp added to each cell, will take time avg
            tmp = np.zeros(t1)
            # overlap count
            o_di_dj = 0
            for i in range(x1):
                for j in range(x2):
                    if i+di >= 0 and i+di < x2 and j+dj >= 0 and j+dj < y2:
                        if mask[i][j] and mask[i+di][j+dj]:
                            mult_out = np.multiply(mat_1[:, i, j], mat_2[:, i+di, j+dj])
                            tmp += mult_out
                            o_di_dj += 1
            A[di, dj] = 0 if o_di_dj==0 else np.average(tmp)/o_di_dj
    #we might need to roll A bc of how I've indexed
    A = np.roll(A, y1-1, axis=0)
    return np.roll(A, x1 - 1, axis=1)

#matrix auto correlation
def mat_a_corr(mat, mask):
    return mat_c_corr(mat, mat, mask)

# TEMPORAL cross CORRELATION
#assumes same size and time dimensions
def td_cross_corr(data_mat1, data_mat2, tmax, mask):
    t, i, j =  data_mat1.shape
    #checks to make sure tmax is in range
    if tmax >= t:
        tmax = t - 1
    C = []
    for ti in range(tmax):
        if ti==0:
            C = [mat_c_corr(data_mat1,data_mat2, mask)]
        s1 = data_mat1[:-ti-1,:,:]
        #print("no shift size: ", s1.shape)
        s2 = data_mat2[ti:-1,:,:]
        #print("shift size: ", s2.shape)
        C.append(mat_c_corr(s1, s2, mask))
    return C

# TEMPORAL auto CORRELATION
def td_auto_corr(data_mat1, tmax, mask):
    return td_cross_corr(data_mat1,data_mat1,tmax, mask)

# creating a running average
def running_avg(mat, t_avg):
    # mat should be (t, x, y)
    # we're goign to return a matrix of the same size
    dt = t_avg//2
    (t, x, y) = mat.shape
    mat_avgs = np.zeros_like(mat)
    for ti in range(t):
        t_m = ti-dt
        t_p = ti+dt+1
        t_min = t_m if t_m >= 0 else 0
        t_max = t_p if t_p <= t else t
        mat_avgs[ti] = np.average(mat[t_min:t_max, :,:], axis=0)
    return mat_avgs

def running_med(mat, t_avg):
    # mat should be (t, x, y)
    # we're goign to return a matrix of the same size
    dt = t_avg//2
    (t, x, y) = mat.shape
    mat_meds = np.zeros_like(mat)
    for ti in range(t):
        t_m = ti-dt
        t_p = ti+dt+1
        t_min = t_m if t_m >= 0 else 0
        t_max = t_p if t_p <= t else t
        mat_meds[ti] = np.median(mat[t_min:t_max, :,:], axis=0)
    return mat_meds
    
