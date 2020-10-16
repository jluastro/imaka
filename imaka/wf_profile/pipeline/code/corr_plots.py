### corr_plots.py
## Eden McEwen
## September 25, 2020
# File contains all plots used in Correlator function

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
    #from array2gif import write_gif
from celluloid import Camera
from IPython.display import Image 

## Helper Function

def unique_file(prefix, suffix):
    i = 0
    while os.path.exists(prefix + "_%s.%s" % (i, suffix)):
        i += 1
    return prefix + "_%s.%s" % (i, suffix)

# Questions: difference between gif_3_mat and gif_3_mat_vmin

##################### Graphing Functions ######################

# Plots three rows of matricies
def graph_3_rows_t_mat(mat_cube_1, mat_cube_2, mat_cube_3, t_list, title = None, label_1=None, label_2=None, label_3=None):
    # each matrix are expected to be in (t,x,y) format
    mat_list = [mat_cube_1, mat_cube_2, mat_cube_3]
    label_list = [label_1, label_2, label_3]
    min_val = np.min(mat_list)
    max_val = np.max(mat_list)
    n_len = 3
    n_width = len(t_list)
    
    fig, axes = plt.subplots(nrows=n_len, ncols=n_width, figsize=(10,6), gridspec_kw={"width_ratios":[1 for n in range(n_width)]})
    fig.subplots_adjust(wspace=0.1)
    fig.suptitle(str(title))
    
    for i in range(n_len):
        for j in range(n_width):
            ax = axes[i][j]
            t = t_list[j]
            try:
                im = ax.matshow(mat_list[i][t])
            except:
                print("error")
            ax.set_xlabel(label_list[i] + " t = " + str(t))  
            ax.set_xticklabels([])
            ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            if i !=0 or j != 0:
                ax.set_yticklabels([])
                ax.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)
    return plt


# Plotting correlations in a triangle
def graph_5_5_ccor(mat_arr, t=0, title = None):
    # each matrix are expected to be in (t,x,y) format
    min_val = np.min(np.array(mat_arr[:,:,t,:,:]))
    max_val = np.max(np.array(mat_arr[:,:,t,:,:]))
    
    n_len, n_width = mat_arr.shape[0], mat_arr.shape[1] 
    
    fig, axes = plt.subplots(nrows=n_len, ncols=n_width, figsize=(2*n_len,2*n_width))
    fig.subplots_adjust(wspace=0.3)
    fig.suptitle(str(title) + ", t=" +str(t))
    
    for i in range(n_len):
        for j in range(n_width):
            if i>=j:
                ax = axes[i][j]
                im = ax.matshow(mat_arr[i][j][t])
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                ax.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)
                if j ==0:
                    ax.set_ylabel("WFS" + str(i))  
                if i == n_len-1:
                    ax.set_xlabel("WFS"+ str(j))  
            else:
                axes[i, j].remove() 
    fig.tight_layout()
    return plt

##################### Animating Functions ######################

# animating three matricies side by side
def gif_3_mat(mat_1, mat_2, mat_3, tmax, title = None, out_file='out.gif', scale_v = False, label_1=None, label_2=None, label_3=None):
    min_val = np.min([mat_1, mat_2, mat_3])
    max_val = np.max([mat_1, mat_2, mat_3])
    fig, (ax,ax2,ax3) = plt.subplots(ncols=3,figsize=(9,3), 
                  gridspec_kw={"width_ratios":[1,1,1]})
    fig.subplots_adjust(wspace=0.3)
    fig.suptitle(str(title))
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.set_xticklabels([])
    ax3.set_yticklabels([])
    ax3.set_xticklabels([])
    
    camera = Camera(fig)

    for t in range(tmax):
        if scale_v:
            im  = ax.matshow(mat_1[t])
            im2 = ax2.matshow(mat_2[t])
            im3 = ax3.matshow(mat_3[t])
        else:
            im  = ax.matshow(mat_1[t], vmin=min_val, vmax=max_val)
            im2 = ax2.matshow(mat_2[t], vmin=min_val, vmax=max_val)
            im3 = ax3.matshow(mat_3[t], vmin=min_val, vmax=max_val)
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        ax.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        ax2.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        ax3.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        ax3.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        ax.text(0.05,1.05," t = " + str(t), bbox=dict(facecolor='white', alpha=0.5), transform=ax.transAxes)
        camera.snap()
    animation = camera.animate()  
    animation.save(out_file, writer = 'imagemagick') 
    
# animating three matricies side by side
def gif_3_mat_vmin(mat_1, mat_2, mat_3, tmax, title = None, out_file='out.gif', label_1=None, label_2=None, label_3=None):
    min_val = np.amin([mat_1, mat_2, mat_3], axis =(0, 2, 3))
    max_val = np.amax([mat_1, mat_2, mat_3], axis =(0, 2, 3))
    fig, (ax,ax2,ax3,cax) = plt.subplots(ncols=4,figsize=(9,3), 
                  gridspec_kw={"width_ratios":[1,1,1,0.05]})
    fig.subplots_adjust(wspace=0.3)
    fig.suptitle(str(title))
    ax.set_xlabel(str(label_1))
    ax2.set_xlabel(str(label_2))
    ax3.set_xlabel(str(label_3))
    #ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.set_xticklabels([])
    ax3.set_yticklabels([])
    ax3.set_xticklabels([])
    
    camera = Camera(fig)

    for t in range(tmax):
        im  = ax.matshow(mat_1[t], vmin=min_val[t], vmax=max_val[t])
        im2 = ax2.matshow(mat_2[t], vmin=min_val[t], vmax=max_val[t])
        im3 = ax3.matshow(mat_3[t], vmin=min_val[t], vmax=max_val[t])
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        ax2.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        ax3.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        ax3.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        ax3.text(0.5,-0.5," t = " + str(t), transform=ax.transAxes)
        cax.cla()
        fig.colorbar(im, cax=cax)
        camera.snap()   
    animation = camera.animate()  
    animation.save(out_file, writer = 'imagemagick')


#animating a 3 by n subarray of matricies
def gif_3_rows_mat(mat_list_1, mat_list_2, mat_list_3, tmax, title = None, out_file='out.gif', label_1=None, label_2=None, label_3=None):
    mat_list = [mat_list_1, mat_list_2, mat_list_3]
    label_list = [label_1, label_2, label_3]
    min_val = np.amin(mat_list, axis=(0, 1, 3, 4))
    max_val = np.amax(mat_list, axis=(0, 1, 3, 4))
    n_len = 3
    n_width = len(mat_list_1)
    fig, axes = plt.subplots(nrows=n_len, ncols=n_width, figsize=(10,6), gridspec_kw={"width_ratios":[1 for n in range(n_width)]})
    fig.subplots_adjust(wspace=0.3)
    fig.suptitle(str(title))
    # setup
    for i in range(n_len):
            for j in range(n_width):
                ax = axes[i][j]
                ax.set_xlabel(label_list[i] + " WFS = " + str(j))
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)
    #animate
    camera = Camera(fig)
    for t in range(tmax):
        for i in range(n_len):
            for j in range(n_width):
                ax = axes[i][j]
                im = ax.matshow(mat_list[i][j][t], vmin = min_val[t], vmax=max_val[t])
                ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        ax.text(0.5,-0.5," t = " + str(t), transform=ax.transAxes)
        camera.snap()
    animation = camera.animate()  
    animation.save(out_file, writer = 'imagemagick') 
    
# Triangular matricies
def gif_ccor_mat(mat_arr, tmax, title=None, out_file='out.gif'):
    min_val = np.min(mat_arr, axis=(0, 1, 3, 4))
    max_val = np.max(mat_arr, axis=(0, 1, 3, 4))
    n_len, n_width = mat_arr.shape[0], mat_arr.shape[1] 
    fig, axes = plt.subplots(nrows=n_len, ncols=n_width, figsize=(2*n_len,2*n_width))
    fig.subplots_adjust(wspace=0.3)
    fig.suptitle(str(title))
    # setup
    for i in range(n_len):
            for j in range(n_width):
                if i>=j:
                    ax = axes[i][j]
                    ax.set_xlabel("WFS" + str(i) + " WFS"+ str(j))  
                    ax.set_xticklabels([])
                    ax.set_yticklabels([])
                    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                    ax.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)
                else:
                    axes[i, j].remove() 
    #animate
    camera = Camera(fig)
    for t in range(tmax):
        for i in range(n_len):
            for j in range(n_width):
                if i>=j:
                    ax = axes[i][j]
                    #im = ax.matshow(mat_arr[i][j][t], vmin = min_val[t], vmax=max_val[t])
                    im = ax.matshow(mat_arr[i][j][t])
                    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                    ax.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)
                if i==0 and j==0:
                    ax.text(0, 1.1," t = " + str(t), transform=ax.transAxes)
        #fig.tight_layout()
        camera.snap()
    animation = camera.animate()  
    animation.save(out_file, writer = 'imagemagick')