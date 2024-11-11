# Temporal_Collateral_Sensitivity
Code associated with the manuscript "Dynamic collateral sensitivity profiles highlight opportunities and challenges for optimizing antibiotic sequences"


"AnalyzeIC50s.m" and its associated functions "shadedErrorBar.m" and "densityplot.m" analyze and create figures related to Figures 2C, S1, S2, S3 using "TempData.mat"

"QuickLook" and "MakeHeterogeneityFigurePanels.m" files use the time-dependent data in the *_day*.mat files to produce figure 3.

"plot_mdp_data.py" will recreate the data seen in figure 5 using the t-MDP.  In addition .pickle files are provided with data to recreate said data. "opt_policy*.pickle" files are much to large for github. Code to create those data sets can be found in "main_td.rs", "lib.rs" and the corresponding cargo file.

"enumerate_plot_all.py" will recreate the data and plots from figure S6ABCD. In addition the .mat file "HeatmapWorkspace.mat" is uploaded here and required for the python code to run.
