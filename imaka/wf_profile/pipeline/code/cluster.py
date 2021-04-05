## Eden McEwen
## March 19, 2021
## Clustering routines to be used with the estimator_r code

import math
import numpy as np

# Using Clustering
import hdbscan
import seaborn as sns
from numpy import unique
from numpy import where
from sklearn.datasets import make_classification
from sklearn.cluster import MeanShift
from sklearn.cluster import DBSCAN
from sklearn.cluster import AffinityPropagation


############################
### Clustering routines ####
############################       
        
def detect_cluster_meanshift(data, vspd, vedgex, vedgey):
    # catch empty detections
    if vspd.shape[0] == 0:
        return [], []
    X = np.array([[vspd[i],vedgex[i],vedgey[i]] for i in range(vspd.shape[0])])
    # define the model
    model = MeanShift()
    # fit model and predict clusters
    yhat = model.fit_predict(X)
    # retrieve unique clusters
    clusters = unique(yhat)
    return clusters, yhat
    
    
def detect_cluster_dbscan(data, vspd, vedgex, vedgey):
    # catch empty detections
    if vspd.shape[0] == 0:
        return [], []
    
    X = np.array([[vspd[i],vedgex[i],vedgey[i]] for i in range(vspd.shape[0])])

    # define the model
    model = DBSCAN(eps=0.30, min_samples=5)
    # fit model and predict clusters
    yhat = model.fit_predict(X)
    # retrieve unique clusters
    clusters = unique(yhat)
    
    return clusters, yhat
    
    
def detect_cluster_affprop(data, vspd, vedgex, vedgey):
    # catch empty detections
    if vspd.shape[0] == 0:
        return [], []
    
    X = np.array([[vspd[i],vedgex[i],vedgey[i]] for i in range(vspd.shape[0])])

    # define the model
    model = AffinityPropagation(damping=0.9)
    # fit the model
    model.fit(X)
    # assign a cluster to each example
    yhat = model.predict(X)
    # retrieve unique clusters
    clusters = unique(yhat)
    
    return clusters, yhat



############### Helpers to see inside a cluster routine

def plot_clusters(data, clusters, dirs, spds, model = ""):
    # printing the 
    for cluster in clusters:
        # get row indexes for samples with this cluster
        row_ix = where(yhat == cluster)
        # scatter plot
        plt.scatter(dirs[row_ix], spds[row_ix])
    # show the plot
    plt.title(model +" \n" + data.name)
    pyplot.show()
    
def print_clusters(clusters, dirs, spds, model = ""):
    print("============ "+model+" ==========")
    # create scatter plot for samples from each cluster
    for cluster in clusters:
        # get row indexes for samples with this cluster
        row_ix = where(yhat == cluster)
        # s
        print("Dir:    ", np.median(dirs[row_ix]))
        ## try with mean as well
        print("  sdv:  ", np.std(dirs[row_ix]))
        print("Spd:    ", np.median(spds[row_ix]))
        print("  sdv:  ", np.std(spds[row_ix]))
        
        # TODO: print on one line, shorten
        print("======================")