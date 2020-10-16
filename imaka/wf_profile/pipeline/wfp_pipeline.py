## wfp_pipeline.py
## Eden McEwen
## Created: June 15th 2020
## Edited: June 25th 2020
## UH REU
# This file takes in a text file of dates, finds them on ehu, and then returns correlation outputs

import os
import sys
import fnmatch

import logging
import threading
import time
import concurrent.futures

sys.path.append('/home/emcewen/code_dev')
from pipeline.cor_pipeline import *
from pipeline.code.log_ex import *
from pipeline.code.file_reader import *

data_path = "../"
out_path = "../"
target_path = "../"
pipe_f = "pipe_iteration"
parallel = False
        
########################### PIPELINE INITS #########################
    
def pipe_run(pipe_funct, names, d_files, o_dirs, t_files):
    """ 
    Applying Pipeline code
        input: lists of names, data paths, output dirs, target paths
        generates: whatever cor_pipeline function called
        output: none
    """
    logging.info("============ STARTING RUN ============")
    start = time.time()
    
    for i, f in enumerate(d_files):
        logging.info("==|==|==|==|==|== Starting file %s of %s ==|==|==|==|==|=="%(i + 1, len(d_files)))
        name = names[i]
        data_f = f
        out_d = o_dirs[i]
        target_f = t_files[i]
        
        #this is what change based on the type of pipeline you're trying to run
        pipe_funct(i, name, data_f, out_d, target_f)
        logging.info("==|==|==|==|==|==> TIME SO FAR: %s"% (time.time() - start))
        
    logging.info("============ RUN FINISHED IN: %s ============"%(time.time() - start))

    
def pipe_run_parallel(pipe_funct, names, d_files, o_dirs, t_files):
    """ 
    Applying Pipeline code in parallel threads
    Uses executor to map 5 workers to each call of pipe iteration
        input: lists of names, data paths, output dirs, target paths
        generates: whatever cor_pipeline function called
        output: none
    """
    num_files = len(d_files)
    logging.info("============ STARTING THREAD RUN ============")
    logging.info("BEGIN THREADS")
    threaded_start = time.time()
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        executor.map(pipe_iteration, range(num_files), names, d_files, o_dirs, t_files)
    
    logging.info("END THREADS")
    time_diff = time.time() - threaded_start
    logging.info("THREAD Time: %s", str(time_diff))
    

########################### PIPE ITERS #########################


def pipe_iteration(i, name, data_f, out_d, target_f):
    """ 
    Applies the file to data_proc_all
    Generates corr fits and corr graphs
    """
    thread_log_handler = start_thread_logging()
    try:
        logging.info("File %s: starting %s"% (i, name))
        logging.info("   %s name: %s"% (i, name))
        logging.info("   %s d_file: %s"% (i, data_f))
        logging.info("   %s o_dirs: %s"% (i, out_d))
        logging.info("   %s t_file: %s"% (i, target_f))
        data_proc_all(name, data_f, out_d)
        data_proc_all(name, data_f, out_d, sub_len=5)
        data_proc_all(name, data_f, out_d, s_sub=True)
        data_proc_all(name, data_f, out_d, s_sub=True, sub_len=5)
        data_proc_all(name, data_f, out_d, s_sub=True, tt_sub=True)
        data_proc_all(name, data_f, out_d, s_sub=True, tt_sub=True, sub_len=5)
    except Exception as e:
        logging.error("Iteration %s error: %s"%(i, e))
    logging.info("File %s: ending", i)
    stop_thread_logging(thread_log_handler)
    
def pipe_data_iteration(i, name, data_f, out_d, target_f):
    """ 
    Applies the file to data_gen_all
    Generates corr fits
    """
    thread_log_handler = start_thread_logging()
    try:
        logging.info("File %s: starting %s"% (i, name))
        logging.info("   %s name: %s"% (i, name))
        logging.info("   %s d_file: %s"% (i, data_f))
        logging.info("   %s o_dirs: %s"% (i, out_d))
        logging.info("   %s t_file: %s"% (i, target_f))
        data_gen_all(name, data_f, out_d, target_f)
        data_gen_all(name, data_f, out_d, target_f, s_sub=True)
        data_gen_all(name, data_f, out_d, target_f, s_sub=True, tt_sub=True)
    except Exception as e:
        logging.error("Iteration %s error: %s"%(i, e))
    logging.info("File %s: ending", i)
    stop_thread_logging(thread_log_handler)
    
def pipe_plot_iteration(i, name, data_f, out_d, target_f):
    """ 
    Applies the file to plot_gen_all
    Generates corr graphs
    """
    thread_log_handler = start_thread_logging()
    try:
        logging.info("File %s: starting %s"% (i, name))
        logging.info("   %s name: %s"% (i, name))
        logging.info("   %s d_file: %s"% (i, data_f))
        logging.info("   %s o_dirs: %s"% (i, out_d))
        logging.info("   %s t_file: %s"% (i, target_f))
        plot_gen_all(name, data_f, out_d)
        plot_gen_all(name, data_f, out_d, sub_len=5)
        plot_gen_all(name, data_f, out_d, s_sub=True)
        plot_gen_all(name, data_f, out_d, s_sub=True, sub_len=5)
        plot_gen_all(name, data_f, out_d, s_sub=True, tt_sub=True)
        plot_gen_all(name, data_f, out_d, s_sub=True, tt_sub=True, sub_len=5)
    except Exception as e:
        logging.error("Iteration %s error: %s"%(i, e))
    logging.info("File %s: ending", i)
    stop_thread_logging(thread_log_handler)
    
########################### HELPER FUNCTS #########################

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
    
    ##### Check input files:
    #TODO
    
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
    
    # Chose pipe_fn
    dispatcher = {'all':pipe_iteration, 'data':pipe_data_iteration, 'plots':pipe_plot_iteration}
    pipe_fn = pipe_iteration
    try:
        pipe_fn = dispatcher[conf_d["pipe_fn"]]
    except KeyError:
        logging.error("No valid function input: using all")
        pipe_fn = pipe_iteration
    
    #TODO: set the file type
    #have both d or e options
    
    # Set if parallel
    if conf_d["parallel"] == "True":
        global parallel
        parallel = True
    
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
    dir_in = run_file.replace(".txt", "_log/")
    if not os.path.exists(dir_in):
        os.makedirs(dir_in)
    config_root_logger(dir_in)
    
    ## Runs Pipeline
    start_run(conf_file, run_file)
    

    
    

    
    