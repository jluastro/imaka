### cor_pipeline.py
# Eden McEwen
# June 18th 2020
# A class for holding all the functions for processing data in the wind profiling pipeline
# Major edit July 10th
# - removed target support
#    deleted set_target
# - set_ssub and set_ttsub changed to set_pre_subs
# - added optional loud arg if want some set up functions ot be quiet
# Major edit September 25th
#  chaged to cor_pipeline
#  seperated out Correlator.py

import time
import os
import sys

import logging
import logging.config

import pandas as pd
import numpy as np

# Self-witten code
from pipeline.code.Correlator import *


def data_gen_all(name, data_f, out_d, t_file=None, tmax=200, s_sub=False, tt_sub = False):
    """ 
    Generate correlation fits files, no graphs generated 
        input: datafile name, aocb file path, out directory, target file, correlation parameters
        generates: one correlation fits
        output: none
    """
    logging.info("======== %s ========" % name)
    # create an object
    curr_data = Correlator(name, data_f, out_d, tmax=tmax, s_sub=s_sub, tt_sub= tt_sub)
    if t_file and not curr_data.target_file :
        curr_data.set_target(t_file)
    try:
        # generate acor
        logging.info("====> Starting acor")
        t0 = time.time()
        curr_data.acor_gen()
        t1 = time.time()
        logging.info("== finished in %s s"% str(t1-t0))
        # generate ccor
        logging.info("====> Starting ccor")
        t2 = time.time()
        curr_data.ccor_gen()
        t3 = time.time()
        logging.info("== finished in %s s"% str(t3-t2))
        # create fits
        logging.info("====> writing fits")
        out = curr_data.fits_write()
        logging.info("== %s"% out)
        logging.info("===> complete")
    except Exception as e:
        logging.warning("====> Data Gen error: %s", e)
        
        
def plot_gen_all(name, data_f, out_d, t_file=None, tmax=200, sub_len=200, avg_sub=True,  s_sub=False, tt_sub = False):
    # create an object
    """ 
    Generates all the useful plots, no fits file saved
        input: datafile name, aocb file path, out directory, target file, correlation parameters, graphing parameters
        generates: acor_avg_gif, acor_png, ccor_all_gif, ccor_all_gif
        output: none
    """
    logging.info("======== %s ========" % name)
    curr_data = Correlator(name, data_f, out_d, tmax=tmax, s_sub=s_sub, tt_sub= tt_sub)
    #if t_file and not curr_data.target_file :
    #    curr_data.set_target(t_file)
    try:
        logging.info("====> Graphing")
        # graph acor
        g_out = curr_data.acor_graph(t_list=[0,5,10,20,30], avg_sub=avg_sub, avg_len=sub_len)
        #time.sleep(2)
        logging.info(g_out)
        g_out = curr_data.acor_animate_avg(dt_max=40, avg_sub=avg_sub, avg_len=sub_len)
        #time.sleep(2)
        logging.info(g_out)
        # graph ccor
        g_out = curr_data.ccor_graph_all(avg_sub=avg_sub, avg_len=sub_len)
        #time.sleep(2)
        logging.info(g_out)
        g_out = curr_data.cor_animate_all(dt_max=40, avg_sub=avg_sub, avg_len=sub_len) 
        #time.sleep(2)
        logging.info(g_out)
        # create fits
        logging.info("====> writing fits")
        out = curr_data.fits_write()
        logging.info("== %s"% out)
        logging.info("===> complete")
    except Exception as e:
        logging.warning("====> Data Proc error: %s", e)


def data_proc_all(name, data_f, out_d, t_file=None, tmax=200, sub_len=0, s_sub=False, tt_sub = False):
    """ 
    Generates correlation fits files and graphs
        input: datafile name, aocb file path, out directory, target file, correlation parameters, graphing parameters
        generates: corr fits file, acor_avg_gif, acor_png, ccor_all_gif, ccor_all_gif
        output: none
    """
    # create an object
    logging.info("======== %s ========" % name)
    logging.info("== sub_len = %s " % sub_len)
    logging.info("== s_sub = %s " % s_sub)
    logging.info("== tt_sub = %s " % tt_sub)
    if not sub_len: sub_len=tmax
    curr_data = Correlator(name, data_f, out_d, tmax=tmax, s_sub=s_sub, tt_sub= tt_sub)
    if t_file and not curr_data.target_file :
        curr_data.set_target(t_file)
    try:
        # generate acor
        logging.info("====> Starting acor")
        t0 = time.time()
        curr_data.acor_gen()
        t1 = time.time()
        logging.info("== finished in %s s"% str(t1-t0))
        # generate ccor
        logging.info("====> Starting ccor")
        t2 = time.time()
        curr_data.ccor_gen()
        t3 = time.time()
        logging.info("== finished in %s s"% str(t3-t2))
        # create fits
        logging.info("====> writing fits")
        out = curr_data.fits_write()
        logging.info("== %s"% out)
        logging.info("====> Graphing")
        # graph acor
        g_out = curr_data.acor_graph(t_list=[0,5,10,20,30], avg_sub=True, avg_len=sub_len)
        #time.sleep(2)
        logging.info(g_out)
        g_out = curr_data.acor_animate_avg(dt_max=40, avg_sub=True, avg_len=sub_len)
        #time.sleep(2)
        logging.info(g_out)
        # graph ccor
        g_out = curr_data.ccor_graph_all(avg_sub=True, avg_len=sub_len)
        #time.sleep(2)
        logging.info(g_out)
        g_out = curr_data.cor_animate_all(40, avg_sub=True, avg_len=sub_len) 
        #time.sleep(2)
        logging.info(g_out)
        logging.info("===> complete")
    except Exception as e:
        logging.warning("====> Data Proc error: %s", e)

        
def data_proc_acor(name, data_f, out_d, t_file=None, tmax=200, s_sub=False, tt_sub = False):
    """ 
    Generates correlation fits files and graphs for ONLY for auto correlations
        input: datafile name, aocb file path, out directory, target file, correlation parameters, graphing parameters
        generates: acorr fits file, acor_avg_gif, acor_png
        output: none
    """
    logging.info("======== %s ========" % name)
    # create an object
    curr_data = Correlator(name, data_f, out_d, tmax=tmax)
    curr_data.set_pre_subs(s_sub=s_sub, tt_sub= tt_sub)
    if t_file and not curr_data.target_file:
        curr_data.set_target(t_file)
    try:
        # generate acor
        logging.info("====> Starting acor")
        t0 = time.time()
        curr_data.acor_gen()
        t1 = time.time()
        logging.info("== finished in %s s"% str(t1-t0))
        # create fits
        logging.info("====> writing fits")
        out = curr_data.fits_write()
        logging.info("== %s"% out)
        logging.info("====> Graphing")
        # graph acor
        curr_data.acor_graph(t_list=[0,5,10,20,30], avg_sub=True, avg_len=tmax)
        curr_data.acor_animate_avg(dt_max=40, avg_sub=True, avg_len=tmax)
        logging.info("===> complete")
    except Exception as e:
        logging.error("====> Data Proc error: %s", e)

        
def data_proc_ccor(name, data_f, out_d, t_file=None, tmax=200, s_sub=False, tt_sub = False):
    """ 
    Generates correlation fits files and graphs for ONLY for cross correlations
        input: datafile name, aocb file path, out directory, target file, correlation parameters, graphing parameters
        generates: acorr fits file, ccor_all_gif, ccor_all_gif
        output: none
    """
    # create an object
    logging.info("======== %s ========" % name)
    curr_data = Correlator(name, data_f, out_d, tmax=tmax)
    curr_data.set_pre_subs(s_sub=s_sub, tt_sub= tt_sub)
    if t_file and not curr_data.target_file :
        curr_data.set_target(t_file)
    try:
        # generate ccor
        logging.info("====> Starting ccor")
        t2 = time.time()
        curr_data.ccor_gen()
        t3 = time.time()
        logging.info("== finished in %s s"% str(t3-t2))
        # create fits
        logging.info("====> writing fits")
        out = curr_data.fits_write()
        logging.info("== %s"% out)
        logging.info("====> Graphing")
        # graph ccor
        curr_data.ccor_graph_all(avg_sub=True, avg_len=tmax)
        curr_data.cor_animate_all(40, avg_sub=True, avg_len=tmax)
        logging.info("===> complete")
    except Exception as e:
        logging.error("====> Data Proc error: %s", e)
        
