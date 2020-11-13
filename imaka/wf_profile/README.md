# Wavefront Profiling Demo
Created by Eden McEwen <emcewen@berkeley.edu> <br>
Summer 2020 REU, [IfA at the University of Hawaii](https://student.ifa.hawaii.edu/reu/) <br>
Mentor: Mark Chun

## Overview

The Imaka Wavefront Profiler takes in imaka Adaptive Optics circular buffers and produces spacial-temporal correlations fits from Wavefront Sensor slopes, as well as graphical representations of those correlations. 

- **README:** This README will give a very basic pointer for Utilization and the main sections of this code. 

- **Interactive HOWTO:** [HOWTO.ipynb](https://github.com/jluastro/imaka/blob/wfp_demo/imaka/wf_profile/HOWTO.ipynb) this jupyter notebook walks through the main pipeline inputs (Note that this is viewabe and useful through the above link even without running the notebook.)

- **Detailed Descriptions** of each programming file, please see the README included in `wf_profile/pipeline`the code directory.


## Basic Utilization:
For running a large batch of files, use the profiler pipeline: `wfp_pipeline.py`
This runs the correlator on the input RUN files

**Example for running this function:**

 `python3 pipeline/wfp_pipeline.py inputs/CONF_ex.txt inputs/RUN_ex.txt`


## Package Structure:

**wf_profile/**
- HOWTO.ipynb
- **inputs/**
  - RUN_ex.py
  - CONF_ex.py
  - **conf/**
  - **runs/**
  - **targets/**
- **pipeline/**
  - cor_pipeline.py
  - est_pipeline.py
  - wfp_pipeline.py
  - **code/**
    - **Correlator.py**
    - corr_plots.py
    - corr_code.py
    - **Estimator.py**
    - est_plots.py
    - data_table.py
    - graph_code.py
    - log_ex.py
    - tpoint.py

### Directory explaination:

#### inputs/

#### inputs/conf

#### inputs/runs

#### inputs/targets


### Main Files explaination:




### [wfp_pipeline.py](https://github.com/jluastro/imaka/blob/wfp_demo/imaka/wf_profile/wp_pipeline.py)
This has the main peice of xecutable code for reading out files and creating fits and plot products from it.
Each file will be searched, verified as existing, and then cross correlations and graphs will be generated from it. <br>
**An example call:** <br>
`$ python3 wfp_pipeline.py -d example_files/ex_RUN.txt`
- `-d` indicates that the file only contains dates in the form YYYYMMDD
- `ex_RUN.txt` contains files seperated by new lines, which are fed into the pipline

### [data_pipeline.py](https://github.com/jluastro/imaka/blob/wfp_demo/imaka/wf_profile/data_pipeline.py)
This program creates a DataPipe object (pedning better name), that takes in a name, datafile, outdirectory, and tmax
You can create an empty object without error, and then read in information from a corresponding fits file. 
Once this object is created, functions in the object can be used to set data or parameters, generate correlations, generate graphs/animations or generate a fits file.
This can, and often is, called without using `wfp_pipeline.py`

### [corr_code.py](https://github.com/jluastro/imaka/blob/wfp_demo/imaka/wf_profile/corr_code.py)
This file holds core cross correlation functions.
It is mainly a supporting code file.

### [graph_code.py](https://github.com/jluastro/imaka/blob/wfp_demo/imaka/wf_profile/graph_code.py)
This file contains supporting graph code used in data_pipeline.py. 
Though it can be called on its own, the most useful implementations of it are already written into data_pipeline.py

## Example input files
If you would like to generate new files, example_files gives a breakdown of what to add. 
