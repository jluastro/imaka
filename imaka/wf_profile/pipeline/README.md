
# Pipeline Code

Detail in a broad strokes how to use each major Pipeline:
- cor_pipeline.py
- est_pipeline.py
- Wfp_pipeline.py

## cor_pipeline.py
(Correlation Code) 

**Description:** A set of functions that allows for abstraction when interacting with batch requests from Correlator.py. Simplifies interactions when requested by wfp_pipeline.py. All functions write to a python logger in the INFO category. 
For simplicity, the default correlation parameters used are: 
	tmax=200, s_sub=False, tt_sub = False 
The default graphing parameters are:
	sub_len=200, avg_sub=True

**Associated py files:**
- Correlator.py

**Functions:**
- data_gen_all
  - Generate correlation fits files 
  - input: datafile name, aocb file path, out directory, target file, correlation parameters
  - generates: one correlation fits
  - output: none
- plot_gen_all
  - Generates all the useful plots, no fits file saved
  - input: datafile name, aocb file path, out directory, target file, correlation parameters
  - generates: acor_avg_gif, acor_png, ccor_all_gif, ccor_all_gif
  - output: none
- data_proc_all
  - Generates correlation fits files and graphs
  - input: datafile name, aocb file path, out directory, target file, correlation parameters, , graphing parameters
  - generates: corr fits file, acor_avg_gif, acor_png, ccor_all_gif, ccor_all_gif
  - output: none
- data_proc_acor
  - Generates correlation fits files and graphs for ONLY for auto correlations
  -  input: datafile name, aocb file path, out directory, target file, correlation parameters, , graphing parameters
  - generates: acorr fits file, acor_avg_gif, acor_png
  - output: none
- Data_proc_ccor
  - Generates correlation fits files and graphs for ONLY for cross correlations
  - input: datafile name, aocb file path, out directory, target file, correlation parameters, graphing parameters
  - generates: acorr fits file, ccor_all_gif, ccor_all_gif
  -  output: none


## est_pipeline.py
(Estimation Code)

**Description:**  Takes in a dataframe with a list of aocb and adds wfs estimates to it with the estimator file.   

**Associated py files:**
- Estimator.py

**Functions:**
- **Df_iterator**
  - Takes in a df and iterates Estimator on each row
  - input: data frame, xcor option, post subtraction length, sigma, threshold
  - output: dataframe with estimates added
- **Df_iterator_mult**
  - Takes in a df and iterates Estimator on each row, with variable options
  - **input:** data frame, xcor options, post subtraction length options, sigma options, threshold options
  - **output:** dataframe with estimates added


## wfp_pipeline.py
(WaveFront Profiler Pipeline Code)

**Description:** The main function for creating correlation files. This interfaces with user input to call on other pipelines. 

**Associated py files:**
- Cor_pipeline
  - All calls to Correlator.py are made through this function
- code/log_ex.py
  - Sets up the logger, or the output of 
- code/file_reader.py

**Functions:**
- pipe_run
  - Applying Pipeline code
- pipe_run_parallel
  - Applying Pipeline code in parallel threads
  - Uses executor to map 5 workers to each call of pipe iteration 
- pipe_iteration
  - tbd
- pipe_data_iteration
  - Makes a call to data_gen_all
  - Will not generate plots
  - Calls raw, ttsub, stt_sub, with 200 and 5 length subs
- pipe_plot_iteration
  - Makes a call to plot_gen_all
  - Will not generate a fits file
  - Calls raw, ttsub, stt_sub, with 200 and 5 length subs
- Set_global_paths
  - Sets the path bases for data, out dir, and targets
  - Uses conf file
- start_run
  - MAIN FUNCTION
- __main__
  - Used when file is called explicitly
  - Takes in a conf file and a run file
  - Ports inputs to start_run
