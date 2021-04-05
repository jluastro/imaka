## estr_pipeline.py
## Eden McEwen
## Created: April 5th 2020
# a pipeline structure for the radial estimator code. 

import os
import sys
import fnmatch

# Personal code
import pipeline.code.Estimator as est


# Filepath defaults
# updated with fn set_global_paths
data_path = "../"
out_path = "../"
target_path = "../"
pipe_f = "pipe_iteration"
parallel = False


# old imports
import numpy as np
import pandas as pd

import math
import time
import os
import sys

#### Change this
#### Output dir
out_dir = "/home/emcewen/out_est/"




################
# preserving old file names
def df_iterator(df, xcor=False, sub_len=5, sig=3, thresh=3, c_size=4):
    return est.df_iterator(df, xcor=False, sub_len=5, sig=3, thresh=3, c_size=4)

def df_iterator_mult(df, sl_list, sig_list, thresh_list, c_list):
    return est.df_iterator_mult(df, sl_list, sig_list, thresh_list, c_list)
    

######################

def set_global_paths(conf_dict):
    #set data_path if available
    if conf_dict["data_path"]:
        global data_path
        data_path = conf_dict["data_path"]
    #set out_path if available
    if conf_dict["out_path"]:
        global out_path
        out_path = conf_dict["out_path"]
    #set target_path if available
    if conf_dict["target_path"]:
        global target_path
        target_path = conf_dict["target_path"]

def start_run(conf_f, run_f):
    """Sets up the pipeline run from a configuration file and a list of date names
    --------
    conf_f : string, file.txt
        Contains file paths 
    run_f : string, file.txt
        Contains dates selected to be run in pipeline
    Outputs
    ----------
    darkname : string
        The name of the output dark FITS file.
    """
    
    ##### Pipeline inputs
    names=[]   # names per each file
    d_files=[] # aocb files paths
    o_dirs=[]  # out dirs per file
    t_files=[] # target dirs per file
    
    ##### CONF FILE: sets  #####
    # Read in configuration file 
    conf_d = read_file_dic(conf_f)
    print(conf_d)
    
    # Set the paths
    set_global_paths(conf_d)
    
    ##### RUN FILE: made into list entries #####
    # look at infile, error if infile 
    entries = read_file(run_f)
    
    #search for all needed files:
    names, d_files, o_dirs, t_files = read_d(entries, data_path, out_path, target_path)
        
    ##### START PIPE
    num_files = len(d_files)
    logging.info("Number of files: %s", str(num_files))
    logging.info('Name of files: %s', names)
    
    #### Applying Pipeline code ###
    if parallel:
        # parallel pipeline init
        pipe_run_parallel(pipe_fn, names, d_files, o_dirs, t_files)
    else:
        # reg pipeline init
        pipe_run(pipe_fn, names, d_files, o_dirs, t_files)


###################################################################
########################### MAIN FUNCTION #########################
###################################################################


if __name__ == '__main__':
    """
    arg1: Conf file (txt file with paths, function, parallel option)
    arg2: text file (with it's extension)
    """
    conf_file = sys.argv[1]
    run_file = sys.argv[2]
    
    ## Starts logger
    #dir_in = run_file.replace(".txt", "_log/")
    #if not os.path.exists(dir_in):
    #    os.makedirs(dir_in)
    #config_root_logger(dir_in)
    
    ## Runs Pipeline
    start_run(conf_file, run_file)