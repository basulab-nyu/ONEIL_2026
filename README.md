This repository contains the code that has been used for ex-vivo data analysis and processing for the following publication: Temporal Coding rather than Circuit Wiring allows Hippocampal CA3 Neurons to Dynamically Distinguish Different Cortical Inputs

All code used to process and analyze the in-vivo data for this publication can be found in the following repository: https://github.com/lukearend/ONeil_2026

System Requirements: 
  ipf files - IgorPro (Wavemetrics) 8.04
    Necessary Files:
      TaroTools - https://sites.google.com/site/tarotoolsregister/
  
  m files - MATLAB (MATHWORKS) 2024a
    Necessary Files:
      abfload - https://www.mathworks.com/matlabcentral/fileexchange/6190-abfload
  
  mat files - MATLAB (MATHWORKS) 2024a

Description of Code Useage: 
AUC_cursors.ipf
  Getting AUC information from top graph, cursor placement required
Compile_Kinetics.m
  Compiling individual kinetics data within a specific folder, must run after indv_kinetics_2025.m
Compile_Kinetics_Binned.m
  Binning compiled kinetics data for dual opsin comparison
Display_Waves_KO.ipf
  Displaying, renaming, and moving traces to specific folders within Igor
Evoked_Events_Multi_Analysis_KO_230923.ipf
  Obtaining kinetics data/values from top graph (can have multiple traces)
OpenABF_wfold.m
  Compiling stimulation information (multi-channel) for an individual folder/cell, creates data structure for each individual cell
  must have abfload downloaded before using
STP_20p_New.m
  Make data structure for compiled STP data
VRKO_IN_Merge.m
  Merging two datasets for interneuron silencing data
Xt_topgraph.ipf
  Obtaining amplitude and timing values for the timing/integration dataset
firing_topgraph.ipf
  Obtaining firing information (rheobase, frequency, etc) from a voltage step protocol
get_waves_from_pxp_KOdataconverted.ipf
  Obtain all folders with abf waves in a specific folder and put into one workspace
get_waves_in_folders_KO.ipf
  Obtain all waves in folders in workspace containing specified phrase/characters
in_timecourse_man.m
  Process interneuron timecourse data
indv_kinetics_2025.m
 Create a data structure for individual kinetics for a single cell, save to specific folder
leak_access.ipf
  Measure the leak and access resistance from a current step protocol 
led_power.mat
  all blue and red LED power values used for this publication, binned and unbinned 
resortStruct.m
  Resort kinetics data for plotting and statistics
save_tables_KO_old.ipf
  Save and compile similar waves in folder as a csv file
spike_prob_CC.m
  Find number of cells that spike and compile data structure from only those cells, run after Compile_Kinetics_Binned.m
timing2025_compile.m
  Compile individual timing data structures in a specific folder, run after timing_2025.m
timing_2025.m
  Create a data structure for individual timing data, save to specific folder
TTX_io.m
  Make data structure for compiled ttx data, save to specific folder
ttx_bin.m
  Bin compiled ttx data for dual opsin comparison, run after TTX_io.m 
