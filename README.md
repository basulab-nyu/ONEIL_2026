This repository contains the code that has been used for ex-vivo data analysis and processing for the following publication: **Temporal Coding rather than Circuit Wiring allows Hippocampal CA3 Neurons to Dynamically Distinguish Different Cortical Inputs** <br>

All code used to process and analyze the in-vivo data for this publication can be found in the following repository: https://github.com/lukearend/ONeil_2026 <br>

**System Requirements:** <br>
&emsp; ipf files - IgorPro (Wavemetrics) 8.04 <br>
&emsp;&emsp; Necessary Files: <br>
&emsp; &emsp; &emsp; TaroTools - https://sites.google.com/site/tarotoolsregister/ <br>
  
&emsp; m files - MATLAB (MATHWORKS) 2024a <br>
&emsp;&emsp; Necessary Files: <br>
&emsp; &emsp; &emsp; abfload - https://www.mathworks.com/matlabcentral/fileexchange/6190-abfload <br>
  
&emsp; mat files - MATLAB (MATHWORKS) 2024a

**Description of Code Useage:** <br>

AUC_cursors.ipf <br>
  &emsp; Getting AUC information from top graph, cursor placement required <br>
	
Compile_Kinetics.m <br>
  &emsp; Compiling individual kinetics data within a specific folder, must run after indv_kinetics_2025.m <br>
	
Compile_Kinetics_Binned.m <br>
  &emsp; Binning compiled kinetics data for dual opsin comparison <br>
	
Display_Waves_KO.ipf <br>
  &emsp; Displaying, renaming, and moving traces to specific folders within Igor <br>
	
Evoked_Events_Multi_Analysis_KO_230923.ipf <br>
  &emsp;Obtaining kinetics data/values from top graph (can have multiple traces) <br>
	
OpenABF_wfold.m <br>
  &emsp; Compiling stimulation information (multi-channel) for an individual folder/cell, creates data structure for each individual cell <br>
  &emsp; must have abfload downloaded before using <br>
	
STP_20p_New.m <br>
  &emsp; Make data structure for compiled STP data <br>
	
VRKO_IN_Merge.m <br>
  &emsp; Merging two datasets for interneuron silencing data <br>
	
Xt_topgraph.ipf <br>
  &emsp; Obtaining amplitude and timing values for the timing/integration dataset <br>
	
firing_topgraph.ipf <br>
  &emsp; Obtaining firing information (rheobase, frequency, etc) from a voltage step protocol <br>
	
get_waves_from_pxp_KOdataconverted.ipf <br>
  &emsp; Obtain all folders with abf waves in a specific folder and put into one workspace <br>
	
get_waves_in_folders_KO.ipf <br>
  &emsp; Obtain all waves in folders in workspace containing specified phrase/characters <br>
	
in_timecourse_man.m <br>
  &emsp; Process interneuron timecourse data <br>
	
indv_kinetics_2025.m <br>
	&emsp; Create a data structure for individual kinetics for a single cell, save to specific folder <br>
	
leak_access.ipf <br>
	&emsp; Measure the leak and access resistance from a current step protocol <br>
	
led_power.mat <br>
&emsp; all blue and red LED power values used for this publication, binned and unbinned <br>

resortStruct.m <br>
&emsp; Resort kinetics data for plotting and statistics <br>

save_tables_KO_old.ipf <br>
&emsp; Save and compile similar waves in folder as a csv file <br>

spike_prob_CC.m <br>
&emsp; Find number of cells that spike and compile data structure from only those cells, run after Compile_Kinetics_Binned.m <br>

timing2025_compile.m <br>
&emsp; Compile individual timing data structures in a specific folder, run after timing_2025.m <br>

timing_2025.m <br>
&emsp; Create a data structure for individual timing data, save to specific folder <br>

TTX_io.m <br>
&emsp; Make data structure for compiled ttx data, save to specific folder <br>

ttx_bin.m <br>
&emsp; Bin compiled ttx data for dual opsin comparison, run after TTX_io.m  <br>
