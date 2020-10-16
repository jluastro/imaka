## Estimator
# Eden McEwen
# September 25, 2020
# Function takes in a correlation file and predicts speed and direction from data.

import math
import time
import os
import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import statistics as stat
from scipy import signal
from astropy.io import fits
import pandas as pd
from astropy.stats import sigma_clipped_stats
from celluloid import Camera

# Using Clustering
import hdbscan
import seaborn as sns

# Self-witten code
from tpoint import TPoint

#generating a mask for corr data
mask_data = dp.mask_8_8_center
mask_cor = signal.correlate2d(mask_data, mask_data)
mask = ~np.array(mask_cor, dtype=bool)

#### Change this
#### Output dir
out_dir = "/home/emcewen/out_est/"

############################ 
######## BASE CLASS ########
############################ 


class Estimator(object):
    columns = ['dataname','DATETIME', "cfht_wspd","cfht_wdir", "250_wspd", "250_wdir", "p_spd", 'p_spd_sdev', "p_dir", 'p_dir_sdev']
    
    def __init__(self, row):
        self.df_in = row
        self.out_fits = self.pull_row('outfits')
        self.data = Correlator("", "", "", f_file=self.out_fits)
        self.name = self.data.name
        self.date = self.data.date
        self.srates = self.pull_srates()
        self.srate = np.average([el for el in self.srates if el >20])
        self.obs_list = self.pull_obs()
        self.min_run_length = 5
        #params, changed each itter
        self.sub_len = 0 
        self.sig = 0
        self.thresh = 0 
        self.c_min = 0
        # Outputs, can be updated
        self.wfs=None
        self.peaks=None
        self.est=None
        self.estimator=None
    
    ### ITER CALL ###
    def gen_wfs(self, xcor=True, sub_len=200):
        self.xcor = xcor
        self.sub_len = sub_len
        self.wfs = self.pull_data()
        return True
        
    def gen_peaks(self, sig=2, thresh=3):
        self.sig = sig
        self.thresh = thresh
        self.peaks = self.peak_find(self.wfs)
        return True
        
    def gen_cluster(self, c_size=4):
        self.c_min = c_size
        estimates, clusterer = self.peak_cluster(self.peaks, self.wfs)
        self.est = estimates
        return True
    
    def get_df(self):
        if self.est == None:
            return False
        return self.df_create(self.est)
    
    ### SINGLE CALL ###
    def estimate(self, xcor, sub_len, sig=2, thresh=3, c_size=4):
        #should contain specifics for estimation
        self.xcor = xcor
        self.sub_len = sub_len 
        self.sig = sig
        self.thresh = thresh
        self.c_min = c_size
        wfs = self.pull_data()
        peaks = self.peak_find(wfs)
        estimates, clusterer = self.peak_cluster(peaks, wfs)
        df = self.df_create(estimates)
        return df
    
    ###### MAIN FUNCTIONS ######
    def pull_data(self):
        #TODO: mask out dead WFS
        xcor = self.xcor
        sub_len = self.sub_len
        s_rates = np.array(self.srates)
        wfs_use = [el>=20 for el in s_rates]
        if xcor:
            print("XCORR")
            x_ccor, y_ccor = self.data.data_get_cc_all(avg_sub=True, avg_len=sub_len)
            avg_ccor = (x_ccor + y_ccor)/2
            ccors_only = np.zeros_like(avg_ccor[0][0])
            count = 0
            for i, row in enumerate(avg_ccor):
                for j in range(len(row)):
                    if i < j and wfs_use[i] and wfs_use[j]:
                        ccors_only = ccors_only + row[j]
                        count = count + 1
            avg_wfs = np.divide(ccors_only, count)
            return avg_wfs
        else:
            print("AUTOCORR")
            x_acor, y_acor = self.data.data_get_ac(avg_sub=True, avg_len=sub_len)
            avg_wfs = np.average((x_acor + y_acor)/2, axis=0)
            return avg_wfs
    
    def peak_find(self, cor_mat):
        sig = self.sig
        thresh = self.thresh
        vals, indxs = max_value_sdev(cor_mat, sigma=sig, thresh=thresh)
        t_points = [[TPoint(t, indxs[t][i][1], indxs[t][i][0], vals[t][i], srate=self.srate) for i in range(len(vals[t]))] for t in range(vals.shape[0])]
        runs = gen_run_strings(t_points, 50, min_len=self.min_run_length)
        runs = clear_noise(runs)
        runs_f = flatten_tpoints(runs)
        return runs_f
    
    def peak_cluster(self, points, wfs):
        # sets up data
        c_size = len(points)//self.c_min +1
        if len(points) == 0:
            return [], None
        dirs_n = [pt.dir_d() for pt in points]
        v = [pt.vel_ms() for pt in points]
        data = np.column_stack((dirs_n, v))
        # clusters
        clusterer = hdbscan.HDBSCAN(min_cluster_size=c_size, min_samples=1).fit(data)
        #PLOT HERE
        threshold = pd.Series(clusterer.outlier_scores_).quantile(0.9)
        non_outliers = np.where(clusterer.outlier_scores_ <= threshold)[0]
        labels = clusterer.labels_[non_outliers]
        points = np.array(points)
        points_used = points[non_outliers]
        #calculating from clusters
        n_clusters = clusterer.labels_.max() + 1
        layer_stats = []
        for n in range(n_clusters):
            m = ma.masked_where(labels != n, points_used)
            masked_tp = np.ma.compressed(m)
            v_n = [pt.vel_ms() for pt in masked_tp]
            dir_n = [pt.dir_d() for pt in masked_tp]
            mean_x = np.mean(dir_n)
            mean_y = np.mean(v_n)
            sd_x =  stat.stdev(dir_n)
            sd_y =  stat.stdev(v_n)
            layer_stats.append([(mean_y, sd_y), (mean_x, sd_x)])
        prms=[self.xcor, self.sub_len, self.sig, self.thresh, self.c_min]
        # spd dir plot
        dirs = [pt.dir_d() for pt in points]
        out_d = out_dir + str(self.date) + "/"
        plot_speed_dir_est(out_d, self.name, dirs, v, clusterer, layer_stats, params=prms)
        # max peak animation
        plot_max_pts(out_d, self.name, [pt.x for pt in points_used], [pt.y for pt in points_used], [pt.vel_ms() for pt in points_used], prms)
        #gif_cor_peaks(wfs, points_used, 50, out_d=out_d, name = self.name, params=prms)
        return layer_stats, clusterer
    
    ### Helper Functions ###
    
    def df_create(self, estimates):
        n_est = len(estimates)
        df_est = pd.DataFrame(columns = self.columns)        
        for n in range(n_est):
            est_list = [estimates[n][0][0], estimates[n][0][1], estimates[n][1][0], estimates[n][1][1]]
            row_list = self.obs_list.copy()
            row_list.extend(est_list)
            new_row = pd.Series(row_list, index = self.columns)
            df_est = df_est.append(new_row, ignore_index=True)
        return df_est
    
    #pull functions
    def pull_row(self, string):
        if string in self.df_in:
            return self.df_in[string]
        else:
            return False
        
    def pull_obs(self):
        name = self.pull_row('dataname')       
        datetime = self.pull_row('DATETIME')
        cfht_spd = self.pull_row("cft_wspd") * 0.5144447
        cfht_dir = self.pull_row("cft_wdir")
        wspd_250 = self.pull_row("250_wspd") * 0.5144447
        wdir_250 = self.pull_row("250_wdir")
        return [name, datetime, cfht_spd, cfht_dir, wspd_250, wdir_250]
        
    def pull_srates(self):
        #want to return a list of all sample rates
        fts = self.out_fits
        hdul = fits.open(fts)
        n_wfs = self.pull_row('nwfs')
        # for each, return 
        srates = [hdul[n].header.get('FSAMPLE') for n in range(1, n_wfs+1)]
        hdul.close()
        return srates
            
############################        
### External functions #####
############################ 

def max_value_sdev(mat_3d, sigma, thresh):
    # n = number of max values
    max_vs = []
    max_idx = []
    for mat in mat_3d:
        mean, median, std = sigma_clipped_stats(mat, sigma=sigma, mask=mask)
        n = np.count_nonzero(np.where(mat >= thresh*std + mean, mat, 0))
        #flatten matrix
        shape = mat.shape
        mat_flat = mat.flatten()
        temp = np.argpartition(-mat_flat, n)
        result_indx = temp[:n]  # indices of highest vals
        temp = np.partition(-mat_flat, n)
        result = -temp[:n] # highest vals
        reslt_idx_real = [np.unravel_index(i,shape) for i in result_indx]
        #update the stored values
        max_vs.append(result)
        max_idx.append(reslt_idx_real)
    return np.array(max_vs), np.array(max_idx)

## Builds strings by location
def gen_run_strings(t_points, len_count, min_len=7):
    #seed the runs:
    if len_count> len(t_points): len_count = len(t_points)
    runs = [[pnt] for pnt in t_points[0]]
    for t in np.arange(1, len_count):
        for point in t_points[t]:
            #check if velocity is zero, or is >100
            if point.vel_ms() > 100 or point.vel_ms() == 0:
                continue  
            #check distances, pick pixel with min distance within 2
            comps = [point.run_comp(run) for run in runs]
            true_idxs = [i for i, x in enumerate(comps) if x]
            if sum(comps) == 0:
                runs.append([point])
            elif sum(comps) == 1:
                idx = true_idxs[0]
                runs[idx].append(point)
            else:
                #if tied, pick the stronger power of the multiple
                closest_id = point.run_closest(runs, true_idxs)
                runs[closest_id].append(point)
    return [run for run in runs if len(run) >= min_len]

def flatten_tpoints(runs):
    if runs == []:
        return runs
    funs_f = runs[0]
    for run in runs[1:]:
        funs_f.extend(run)
    return funs_f

def clear_noise(runs):
    our_runs = []
    for run in runs:
        points = [[run[0], 0]]
        for point in run:
            seen = False
            for i in range(len(points)):
                if point.dist(points[i][0]) == 0 and not seen:
                    points[i][1] = 1 + points[i][1]
                    seen = True
            if not seen:
                points.append([point, 1])
        counts = np.array([p[1] for p in points])
        check1 = np.sum(counts > len(run)*.8)
        check2 = np.sum(counts > 10)
        if not check1 and not check2:
            our_runs.append(run)
    return our_runs