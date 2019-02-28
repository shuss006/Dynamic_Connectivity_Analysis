# Sana Hussain

Hello,

My name is Sana Hussain and I am a 3rd year PhD student at UC Riverside.  I currently work under the advisement of Dr. Xiaoping Hu and Dr. Megan Peters.  As part of my graduate studies I have written the above scripts in MATLAB which are coupled as described below:

The first set of scripts (dynamic_connectivity_analysis.m and rest_filter.m) analyze the dynamic functional connectivity between networks' BOLD signals.  Specifically, this code utilizes inputs of fMRI data and ROI masks created in FSL in the form of nifti files to acquire the static and dynamic connectivity measurements between the ROIs in a network.  Included in the script is also code for a series of plots to visualize the results.

The second set of scripts (Bayesian_analysis.m and Bayesian_analysis_function.m) employ Bayes' rule to determine whether a subject squeezed a squeeze-ball given his/her network connectivity patterns.

Please feel free to contact me with any questions. Thank you!

Update (Feb. 28th, 2019): I have made the code in the dynamic_connectivity_analysis.m more efficient by splitting up the analyses to be performed into different functions.  The analyses and reasons driving analyses are the same, it is now just a bit easier now to call the functions.  The m files used in this set of code are: main_script.m, extract_bold_signal.m, static_correlation.m, dynamic_connectivity.m, and LC_analysis.m






