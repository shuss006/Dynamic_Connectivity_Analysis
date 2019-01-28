%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
%                                                                         %
%                            Bayesian Analysis                            %
%                           Author: Sana Hussain                          %
%                    University of California, Riverside                  %
%                                                                         %
%-------------------------------------------------------------------------%

%This code utilizes inputs of fMRI data, specifically the dynamic
%connectivity matrices created from the dynamic_connectivity_analysis.m
%script. These connectivity matrices are inputted into 
%Bayesian_analysis_function.m where probability density functions of
%the number of times a connectivity pattern occurred is created. Then
%Bayesian analysis is used to find the probability that the subject
%squeezed a squeeze-ball (the details of the experiment are described below)
%given the before and after connectivity patterns of four networks. The 
%outputs from that analysis are used to generate the x and y values of
%an ROC curve, and to compute the area under the curve (AUC) of that curve.

%I investigated the connectivity pattern of a total of four networks (DMN, 
%FPCN, DAN, and SN) and examined them in isolation, in pairs, in triplets 
%and all together.  Upon investigating four, three, and two networks at a 
%time, multivariable kernel density estimation was used to create a four-, 
%three-, and two- dimensional probability density function, respectively, 
%using an arbitrary bandwidth value. Details of this process are described
%in Bayesian_analysis_function.m.

%The fMRI data acquired is a pseudo-resting state design were subjects
%underwent a squeezing task followed by a subsequent "resting state" block.
%The paradigm is as follows: 5 minutes of resting state, 30 seconds of
%squeezing a squeezeball, and 5 minutes of resting state. With a TR of 2s,
%each resting state block is 150 time points, and each squeezing period is
%15 time points. To investigate each resting state epoch separately, I
%divded the time series into each of these four blocks and analyzed each
%one separately.

%Commonly used phrases:
%"before"           = before the squeezing session (first resting state epoch)   
%"after"            = after the squeezing session (second resting state epoch)
%"overlapping"      = dynamic connectivity analysis using a sliding time window
%                     approach and overlapping sliding windows
%"nonoverlapping"   = dynamic connectivity analysis using a sliding time
%                     window approach and nonoverlapping sliding windows

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all

%Load the within-network dynamic connectivity values for all four networks
%for before and after the squeezing session.
%Details into how the connectivity values were analytically obtained can be
%found in the dynamic_connectivity_analysis.m script.

load avg_ROI_dynamic_DMN_before
load avg_ROI_dynamic_FPCN_before
load avg_ROI_dynamic_DAN_before
load avg_ROI_dynamic_SN_before

load avg_ROI_dynamic_DMN_after
load avg_ROI_dynamic_FPCN_after
load avg_ROI_dynamic_DAN_after
load avg_ROI_dynamic_SN_after

load avg_ROI_dynamic_DMN_before_nonoverlap
load avg_ROI_dynamic_FPCN_before_nonoverlap
load avg_ROI_dynamic_DAN_before_nonoverlap
load avg_ROI_dynamic_SN_before_nonoverlap

load avg_ROI_dynamic_DMN_after_nonoverlap
load avg_ROI_dynamic_FPCN_after_nonoverlap
load avg_ROI_dynamic_DAN_after_nonoverlap
load avg_ROI_dynamic_SN_after_nonoverlap

load raw_LC_right_ts
load raw_LC_left_ts

%%

%Now I input the before and after connectivity patterns in to the
%corresponding function for all four networks, three networks, two
%networks, and each individual network to utilize Bayes' rule. This gives
%the probability of whether the subject squeezed a squeeze-ball given the
%connectivity patterns of each network. This analysis would also provide 
%insight into whether a single network, or a combination of networks, is 
%dominating or dragging down the ability of our classifier to remain accurate.

%% ALL NETWORKS

all_before = [avg_ROI_dynamic_DMN_before avg_ROI_dynamic_FPCN_before avg_ROI_dynamic_DAN_before avg_ROI_dynamic_SN_before];
all_after = [avg_ROI_dynamic_DMN_after avg_ROI_dynamic_FPCN_after avg_ROI_dynamic_DAN_after avg_ROI_dynamic_SN_after];
all_bw = [4.0579 10.7345 4.4185 11.5466];

[output_all] = Bayesian_analysis_function(all_before,all_after,all_bw); 

all_Bayes_Rule = output_all.Bayes_Rule;
perfcurve_all_x = output_all.perfcurve_x;
perfcurve_all_y = output_all.perfcurve_y;
perfcurve_all_AUC = output_all.perfcurve_AUC;

save('all_Bayes_Rule.mat','all_Bayes_Rule');
save('perfcurve_all_x.mat','perfcurve_all_x');
save('perfcurve_all_y.mat','perfcurve_all_y');
save('perfcurve_all_AUC.mat','perfcurve_all_AUC');


%% DMN and FPCN and DAN

DMN_FPCN_DAN_before = [avg_ROI_dynamic_DMN_before avg_ROI_dynamic_FPCN_before avg_ROI_dynamic_DAN_before];
DMN_FPCN_DAN_after = [avg_ROI_dynamic_DMN_after avg_ROI_dynamic_FPCN_after avg_ROI_dynamic_DAN_after];
three_bw = [10.7345 4.4185 11.5466];

[output_DMN_FPCN_DAN] = Bayesian_analysis_function(DMN_FPCN_DAN_before,DMN_FPCN_DAN_after, three_bw); 

DMN_FPCN_DAN_Bayes_Rule = output_DMN_FPCN_DAN.Bayes_Rule;
perfcurve_DMN_FPCN_DAN_x = output_DMN_FPCN_DAN.perfcurve_x;
perfcurve_DMN_FPCN_DAN_y = output_DMN_FPCN_DAN.perfcurve_y;
perfcurve_DMN_FPCN_DAN_AUC = output_DMN_FPCN_DAN.perfcurve_AUC;

save('DMN_FPCN_DAN_Bayes_Rule.mat','DMN_FPCN_DAN_Bayes_Rule');
save('perfcurve_DMN_FPCN_DAN_x.mat','perfcurve_DMN_FPCN_DAN_x');
save('perfcurve_DMN_FPCN_DAN_y.mat','perfcurve_DMN_FPCN_DAN_y');
save('perfcurve_DMN_FPCN_DAN_AUC.mat','perfcurve_DMN_FPCN_DAN_AUC');


%% DMN and FPCN and SN

DMN_FPCN_SN_before = [avg_ROI_dynamic_DMN_before avg_ROI_dynamic_FPCN_before avg_ROI_dynamic_SN_before];
DMN_FPCN_SN_after = [avg_ROI_dynamic_DMN_after avg_ROI_dynamic_FPCN_after avg_ROI_dynamic_SN_after];
three_bw = [10.7345 4.4185 11.5466];

[output_DMN_FPCN_SN] = Bayesian_analysis_function(DMN_FPCN_SN_before,DMN_FPCN_SN_after, three_bw); 

DMN_FPCN_SN_Bayes_Rule = output_DMN_FPCN_SN.Bayes_Rule;
perfcurve_DMN_FPCN_SN_x = output_DMN_FPCN_SN.perfcurve_x;
perfcurve_DMN_FPCN_SN_y = output_DMN_FPCN_SN.perfcurve_y;
perfcurve_DMN_FPCN_SN_AUC = output_DMN_FPCN_SN.perfcurve_AUC;

save('DMN_FPCN_SN_Bayes_Rule.mat','DMN_FPCN_SN_Bayes_Rule');
save('perfcurve_DMN_FPCN_SN_x.mat','perfcurve_DMN_FPCN_SN_x');
save('perfcurve_DMN_FPCN_SN_y.mat','perfcurve_DMN_FPCN_SN_y');
save('perfcurve_DMN_FPCN_SN_AUC.mat','perfcurve_DMN_FPCN_SN_AUC');

%% DMN and DAN and SN

DMN_DAN_SN_before = [avg_ROI_dynamic_DMN_before avg_ROI_dynamic_DAN_before avg_ROI_dynamic_SN_before];
DMN_DAN_SN_after = [avg_ROI_dynamic_DMN_after avg_ROI_dynamic_DAN_after avg_ROI_dynamic_SN_after];
three_bw = [10.7345 4.4185 11.5466];

[output_DMN_DAN_SN] = Bayesian_analysis_function(DMN_DAN_SN_before,DMN_DAN_SN_after, three_bw); 

DMN_DAN_SN_Bayes_Rule = output_DMN_DAN_SN.Bayes_Rule;
perfcurve_DMN_DAN_SN_x = output_DMN_DAN_SN.perfcurve_x;
perfcurve_DMN_DAN_SN_y = output_DMN_DAN_SN.perfcurve_y;
perfcurve_DMN_DAN_SN_AUC = output_DMN_DAN_SN.perfcurve_AUC;

save('DMN_DAN_SN_Bayes_Rule.mat','DMN_DAN_SN_Bayes_Rule');
save('perfcurve_DMN_DAN_SN_x.mat','perfcurve_DMN_DAN_SN_x');
save('perfcurve_DMN_DAN_SN_y.mat','perfcurve_DMN_DAN_SN_y');
save('perfcurve_DMN_DAN_SN_AUC.mat','perfcurve_DMN_DAN_SN_AUC');

%% FPCN and DAN and SN

FPCN_DAN_SN_before = [avg_ROI_dynamic_FPCN_before avg_ROI_dynamic_DAN_before avg_ROI_dynamic_SN_before];
FPCN_DAN_SN_after = [avg_ROI_dynamic_FPCN_after avg_ROI_dynamic_DAN_after avg_ROI_dynamic_SN_after];
three_bw = [10.7345 4.4185 11.5466];

[output_FPCN_DAN_SN] = Bayesian_analysis_function(FPCN_DAN_SN_before,FPCN_DAN_SN_after,three_bw); 

FPCN_DAN_SN_Bayes_Rule = output_FPCN_DAN_SN.Bayes_Rule;
perfcurve_FPCN_DAN_SN_x = output_FPCN_DAN_SN.perfcurve_x;
perfcurve_FPCN_DAN_SN_y = output_FPCN_DAN_SN.perfcurve_y;
perfcurve_FPCN_DAN_SN_AUC = output_FPCN_DAN_SN.perfcurve_AUC;

save('FPCN_DAN_SN_Bayes_Rule.mat','FPCN_DAN_SN_Bayes_Rule');
save('perfcurve_FPCN_DAN_SN_x.mat','perfcurve_FPCN_DAN_SN_x');
save('perfcurve_FPCN_DAN_SN_y.mat','perfcurve_FPCN_DAN_SN_y');
save('perfcurve_FPCN_DAN_SN_AUC.mat','perfcurve_FPCN_DAN_SN_AUC');

%% DMN and FPCN 

DMN_FPCN_before = [avg_ROI_dynamic_DMN_before avg_ROI_dynamic_FPCN_before];
DMN_FPCN_after = [avg_ROI_dynamic_DMN_after avg_ROI_dynamic_FPCN_after];
two_bw = [4.4185 11.5466];

[output_DMN_FPCN] = Bayesian_analysis_function(DMN_FPCN_before,DMN_FPCN_after, two_bw); 

DMN_FPCN_Bayes_Rule = output_DMN_FPCN.Bayes_Rule;
perfcurve_DMN_FPCN_x = output_DMN_FPCN.perfcurve_x;
perfcurve_DMN_FPCN_y = output_DMN_FPCN.perfcurve_y;
perfcurve_DMN_FPCN_AUC = output_DMN_FPCN.perfcurve_AUC;

save('DMN_FPCN_Bayes_Rule.mat','DMN_FPCN_Bayes_Rule');
save('perfcurve_DMN_FPCN_x.mat','perfcurve_DMN_FPCN_x');
save('perfcurve_DMN_FPCN_y.mat','perfcurve_DMN_FPCN_y');
save('perfcurve_DMN_FPCN_AUC.mat','perfcurve_DMN_FPCN_AUC');

%% DMN and DAN

DMN_DAN_before = [avg_ROI_dynamic_DMN_before avg_ROI_dynamic_DAN_before];
DMN_DAN_after = [avg_ROI_dynamic_DMN_after avg_ROI_dynamic_DAN_after];
two_bw = [4.4185 11.5466];

[output_DMN_DAN] = Bayesian_analysis_function(DMN_DAN_before,DMN_DAN_after, two_bw); 

DMN_DAN_Bayes_Rule = output_DMN_DAN.Bayes_Rule;
perfcurve_DMN_DAN_x = output_DMN_DAN.perfcurve_x;
perfcurve_DMN_DAN_y = output_DMN_DAN.perfcurve_y;
perfcurve_DMN_DAN_AUC = output_DMN_DAN.perfcurve_AUC;

save('DMN_DAN_Bayes_Rule.mat','DMN_DAN_Bayes_Rule');
save('perfcurve_DMN_DAN_x.mat','perfcurve_DMN_DAN_x');
save('perfcurve_DMN_DAN_y.mat','perfcurve_DMN_DAN_y');
save('perfcurve_DMN_DAN_AUC.mat','perfcurve_DMN_DAN_AUC');

%% DMN and SN

DMN_SN_before = [avg_ROI_dynamic_DMN_before avg_ROI_dynamic_SN_before];
DMN_SN_after = [avg_ROI_dynamic_DMN_after avg_ROI_dynamic_SN_after];
two_bw = [4.4185 11.5466];

[output_DMN_SN] = Bayesian_analysis_function(DMN_SN_before,DMN_SN_after, two_bw); 

DMN_SN_Bayes_Rule = output_DMN_SN.Bayes_Rule;
perfcurve_DMN_SN_x = output_DMN_SN.perfcurve_x;
perfcurve_DMN_SN_y = output_DMN_SN.perfcurve_y;
perfcurve_DMN_SN_AUC = output_DMN_SN.perfcurve_AUC;

save('DMN_SN_Bayes_Rule.mat','DMN_SN_Bayes_Rule');
save('perfcurve_DMN_SN_x.mat','perfcurve_DMN_SN_x');
save('perfcurve_DMN_SN_y.mat','perfcurve_DMN_SN_y');
save('perfcurve_DMN_SN_AUC.mat','perfcurve_DMN_SN_AUC');

%% FPCN and SN

FPCN_SN_before = [avg_ROI_dynamic_FPCN_before avg_ROI_dynamic_SN_before];
FPCN_SN_after = [avg_ROI_dynamic_FPCN_after avg_ROI_dynamic_SN_after];
two_bw = [4.4185 11.5466];

[output_FPCN_SN] = Bayesian_analysis_function(FPCN_SN_before,FPCN_SN_after,two_bw); 

FPCN_SN_Bayes_Rule = output_FPCN_SN.Bayes_Rule;
perfcurve_FPCN_SN_x = output_FPCN_SN.perfcurve_x;
perfcurve_FPCN_SN_y = output_FPCN_SN.perfcurve_y;
perfcurve_FPCN_SN_AUC = output_FPCN_SN.perfcurve_AUC;

save('FPCN_SN_Bayes_Rule.mat','FPCN_SN_Bayes_Rule');
save('perfcurve_FPCN_SN_x.mat','perfcurve_FPCN_SN_x');
save('perfcurve_FPCN_SN_y.mat','perfcurve_FPCN_SN_y');
save('perfcurve_FPCN_SN_AUC.mat','perfcurve_FPCN_SN_AUC');

%% FPCN and DAN 

FPCN_DAN_before = [avg_ROI_dynamic_FPCN_before avg_ROI_dynamic_DAN_before];
FPCN_DAN_after = [avg_ROI_dynamic_FPCN_after avg_ROI_dynamic_DAN_after];
two_bw = [4.4185 11.5466];

[output_FPCN_DAN] = Bayesian_analysis_function(FPCN_DAN_before,FPCN_DAN_after,two_bw); 

FPCN_DAN_Bayes_Rule = output_FPCN_DAN.Bayes_Rule;
perfcurve_FPCN_DAN_x = output_FPCN_DAN.perfcurve_x;
perfcurve_FPCN_DAN_y = output_FPCN_DAN.perfcurve_y;
perfcurve_FPCN_DAN_AUC = output_FPCN_DAN.perfcurve_AUC;

save('FPCN_DAN_Bayes_Rule.mat','FPCN_DAN_Bayes_Rule');
save('perfcurve_FPCN_DAN_x.mat','perfcurve_FPCN_DAN_x');
save('perfcurve_FPCN_DAN_y.mat','perfcurve_FPCN_DAN_y');
save('perfcurve_FPCN_DAN_AUC.mat','perfcurve_FPCN_DAN_AUC');

%% DAN and SN

DAN_SN_before = [avg_ROI_dynamic_DAN_before avg_ROI_dynamic_SN_before];
DAN_SN_after = [avg_ROI_dynamic_DAN_after avg_ROI_dynamic_SN_after];
two_bw = [4.4185 11.5466];

[output_DAN_SN] = Bayesian_analysis_function(DAN_SN_before,DAN_SN_after, two_bw); 

DAN_SN_Bayes_Rule = output_DAN_SN.Bayes_Rule;
perfcurve_DAN_SN_x = output_DAN_SN.perfcurve_x;
perfcurve_DAN_SN_y = output_DAN_SN.perfcurve_y;
perfcurve_DAN_SN_AUC = output_DAN_SN.perfcurve_AUC;

save('DAN_SN_Bayes_Rule.mat','DAN_SN_Bayes_Rule');
save('perfcurve_DAN_SN_x.mat','perfcurve_DAN_SN_x');
save('perfcurve_DAN_SN_y.mat','perfcurve_DAN_SN_y');
save('perfcurve_DAN_SN_AUC.mat','perfcurve_DAN_SN_AUC');

%% DMN

labels = [zeros(90,1); ones(90,1)]; %Before squeezing is 0 and after squeezing is 1
[perfcurve_DMN_x, perfcurve_DMN_y, ~, perfcurve_DMN_AUC] = perfcurve(labels, [avg_ROI_dynamic_DMN_before;avg_ROI_dynamic_DMN_after], 1);

save('perfcurve_DMN_x.mat','perfcurve_DMN_x');
save('perfcurve_DMN_y.mat','perfcurve_DMN_y');
save('perfcurve_DMN_AUC.mat','perfcurve_DMN_AUC');

%% FPCN

labels = [zeros(90,1); ones(90,1)]; %Before squeezing is 0 and after squeezing is 1
[perfcurve_FPCN_x, perfcurve_FPCN_y, ~, perfcurve_FPCN_AUC] = perfcurve(labels, [avg_ROI_dynamic_FPCN_before;avg_ROI_dynamic_FPCN_after], 1);

save('perfcurve_FPCN_x.mat','perfcurve_FPCN_x');
save('perfcurve_FPCN_y.mat','perfcurve_FPCN_y');
save('perfcurve_FPCN_AUC.mat','perfcurve_FPCN_AUC');

%% DAN

labels = [zeros(90,1); ones(90,1)]; %Before squeezing is 0 and after squeezing is 1
[perfcurve_DAN_x, perfcurve_DAN_y, ~, perfcurve_DAN_AUC] = perfcurve(labels, [avg_ROI_dynamic_DAN_before;avg_ROI_dynamic_DAN_after], 1);

save('perfcurve_DAN_x.mat','perfcurve_DAN_x');
save('perfcurve_DAN_y.mat','perfcurve_DAN_y');
save('perfcurve_DAN_AUC.mat','perfcurve_DAN_AUC');

%% SN

labels = [zeros(90,1); ones(90,1)]; %Before squeezing is 0 and after squeezing is 1
[perfcurve_SN_x, perfcurve_SN_y, ~, perfcurve_SN_AUC] = perfcurve(labels, [avg_ROI_dynamic_SN_before;avg_ROI_dynamic_SN_after], 1);

save('perfcurve_SN_x.mat','perfcurve_SN_x');
save('perfcurve_SN_y.mat','perfcurve_SN_y');
save('perfcurve_SN_AUC.mat','perfcurve_SN_AUC');

%% DMN Non-Overlapping Windows

labels_nonoverlap = [zeros(13,1); ones(13,1)]; %Before squeezing is 0 and after squeezing is 1
[perfcurve_DMN_nonoverlap_x, perfcurve_DMN_nonoverlap_y, ~, perfcurve_DMN_nonoverlap_AUC] = perfcurve(labels_nonoverlap, [avg_ROI_dynamic_DMN_before_nonoverlap;avg_ROI_dynamic_DMN_after_nonoverlap], 1);

save('perfcurve_DMN_nonoverlap_x.mat','perfcurve_DMN_nonoverlap_x');
save('perfcurve_DMN_nonoverlap_y.mat','perfcurve_DMN_nonoverlap_y');
save('perfcurve_DMN_nonoverlap_AUC.mat','perfcurve_DMN_nonoverlap_AUC');

%% FPCN Non-Overlapping Windows

[perfcurve_FPCN_nonoverlap_x, perfcurve_FPCN_nonoverlap_y, ~, perfcurve_FPCN_nonoverlap_AUC] = perfcurve(labels_nonoverlap, [avg_ROI_dynamic_FPCN_before_nonoverlap;avg_ROI_dynamic_FPCN_after_nonoverlap], 1);

save('perfcurve_FPCN_nonoverlap_x.mat','perfcurve_FPCN_nonoverlap_x');
save('perfcurve_FPCN_nonoverlap_y.mat','perfcurve_FPCN_nonoverlap_y');
save('perfcurve_FPCN_nonoverlap_AUC.mat','perfcurve_FPCN_nonoverlap_AUC');

%% DAN Non-Overlapping Windows

[perfcurve_DAN_nonoverlap_x, perfcurve_DAN_nonoverlap_y, ~, perfcurve_DAN_nonoverlap_AUC] = perfcurve(labels_nonoverlap, [avg_ROI_dynamic_DAN_before_nonoverlap;avg_ROI_dynamic_DAN_after_nonoverlap], 1);

save('perfcurve_DAN_nonoverlap_x.mat','perfcurve_DAN_nonoverlap_x');
save('perfcurve_DAN_nonoverlap_y.mat','perfcurve_DAN_nonoverlap_y');
save('perfcurve_DAN_nonoverlap_AUC.mat','perfcurve_DAN_nonoverlap_AUC');

%% SN Non-Overlapping Windows

[perfcurve_SN_nonoverlap_x, perfcurve_SN_nonoverlap_y, ~, perfcurve_SN_nonoverlap_AUC] = perfcurve(labels_nonoverlap, [avg_ROI_dynamic_SN_before_nonoverlap;avg_ROI_dynamic_SN_after_nonoverlap], 1);

save('perfcurve_SN_nonoverlap_x.mat','perfcurve_SN_nonoverlap_x');
save('perfcurve_SN_nonoverlap_y.mat','perfcurve_SN_nonoverlap_y');
save('perfcurve_SN_nonoverlap_AUC.mat','perfcurve_SN_nonoverlap_AUC');

%% Right LC Non-Overlapping Windows
labels_LC = [zeros(150,1); ones(150,1)]; %Before squeezing is 0 and after squeezing is 1

[perfcurve_right_LC_x, perfcurve_right_LC_y, ~, perfcurve_right_LC_AUC] = perfcurve(labels_LC, [raw_LC_right_ts(1:150,:); raw_LC_right_ts(166:315,:)], 1);

save('perfcurve_right_LC_x.mat','perfcurve_right_LC_x');
save('perfcurve_right_LC_y.mat','perfcurve_right_LC_y');
save('perfcurve_right_LC_AUC.mat','perfcurve_right_LC_AUC');

%% Left LC Non-Overlapping Windows
labels_LC = [zeros(150,1); ones(150,1)]; %Before squeezing is 0 and after squeezing is 1

[perfcurve_left_LC_x, perfcurve_left_LC_y, ~, perfcurve_left_LC_AUC] = perfcurve(labels_LC, [raw_LC_left_ts(1:150,:); raw_LC_left_ts(166:315,:)], 1);

save('perfcurve_left_LC_x.mat','perfcurve_left_LC_x');
save('perfcurve_left_LC_y.mat','perfcurve_left_LC_y');
save('perfcurve_left_LC_AUC.mat','perfcurve_left_LC_AUC');
