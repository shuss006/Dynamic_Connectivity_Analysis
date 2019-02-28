%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
%                                                                         %
%                  Static and Dynamic Connectivity Analysis               %
%                           Author: Sana Hussain                          %
%                    University of California, Riverside                  %
%                                                                         %
%-------------------------------------------------------------------------%

%This script runs exactly the same analyses as those mentioned in
%dynamic_connectivity_analysis.m  I have split the code into functions to 
%make it more efficient, easier to understand, and reduce the number of 
%lines in a single script.  The functions called are: extract_bold_signal.m,
%static_correlation.m, dynamic_connectivity.m, LC_analysis.m.  The ideas,
%paradigm, and methods are exactly the same as those described in
%dynamic_connectivity_analysis.m, which is reiterated below.

%This code utilizes inputs of fMRI data and ROI masks created in FSL 
%in the form of nifti files to acquire the static and dynamic connectivity
%measurements between the ROIs in a network. These connectivity values
%are then averaged to ascertain the overall within-network connectivity of
%four networks: the default mode network (DMN), the fronto-parietal control
%network (FPCN), the dorsal attention network (DAN), and the salience
%network (SN).

%The fMRI data acquired is a pseudo-resting state design were subjects
%underwent a squeezing task followed by a subsequent "resting state" block.
%The paradigm is as follows: four 3 minute resting state epochs interspersed
%with three 30 second periods of squeezing a squeeze-ball. With a TR of 2s,
%each resting state block is 90 time points, and each squeezing period is
%15 time points. To investigate each resting state epoch separately, I
%divded the time series into each of these four blocks and analyzed each
%one separately.

%The locus coeruleus (LC) is also of interest in this experiment so the BOLD
%signal using an LC mask is acquired and eye tracking data (pupillometry data)
%was also acquired as a means of checking the reliability of the BOLD signal. 

%The NIFTI toolbox is needed to run this code


%Main Script where we call all the other functions

close all
clear all

test = load_untouch_nii('default_mode_network_ROIs.nii');
ROI_mask_DMN = test.img;

test = load_untouch_nii('fronto_parietal_control_network__ROIs.nii');
ROI_mask_FPCN = test.img;

test = load_untouch_nii('dorsal_attention_network_ROIs.nii');
ROI_mask_DAN = test.img;

test = load_untouch_nii('salience_network_ROIs.nii');
ROI_mask_SN = test.img;

test = load_untouch_nii('new_LC_atlas_transformed_bin_left.nii');
LC_mask_left = test.img;

test = load_untouch_nii('new_LC_atlas_transformed_bin_right.nii');
LC_mask_right = test.img;

test = load_untouch_nii('csf_mask_adjusted_bin.nii');
csf_mask = test.img;

test = load_untouch_nii('MNI_fourth_ventricle_mask_bin.nii');
fourth_ventricle_mask = test.img;

test = load_untouch_nii('data_in_MNI.nii');
timeseries = test.img;

%Log File Stuff
log_file = load('log_file.txt');

[nx, ny, nz, nt] = size(timeseries); %The time series for all 9 ROIs
ROI_DMN = max(ROI_mask_DMN(:));
ROI_FPCN = max(ROI_mask_FPCN(:));
ROI_DAN = max(ROI_mask_DAN(:));
ROI_SN = max(ROI_mask_SN(:));

%Divide up the time series into when the resting state epochs and squeezing
%sessions occurred
rs_epoch_0 = 1:1:120;
squeeze_session_1 = 121:1:135;
rs_epoch_1 = 136:1:255;
squeeze_session_2 = 256:1:270;
rs_epoch_2 = 271:1:390;
squeeze_session_3 = 391:1:405;
rs_epoch_3 = 406:1:525;
squeeze_session_4 = 526:1:540;
rs_epoch_4 = 541:1:660;

%Define other pertinent variables
total_blocks = 9;
rs_blocks_only = 5;
squeeze_sessions_only = 4;
all_time = {rs_epoch_0; squeeze_session_1; rs_epoch_1; squeeze_session_2; rs_epoch_2; squeeze_session_3; rs_epoch_3; squeeze_session_4; rs_epoch_4};
rs_time_only = {rs_epoch_0; rs_epoch_1; rs_epoch_2; rs_epoch_3; rs_epoch_4};
squeeze_time_only = {squeeze_session_1; squeeze_session_2; squeeze_session_3; squeeze_session_4};
rs_tr = length(rs_epoch_0);

%Find the % MVC during each squeeze session
for i = 1:squeeze_sessions_only
   mag_squeeze(i) = mean(log_file(squeeze_time_only{i,1},3));    
end

%Extract the BOLD Signal from each ROI within each network
DMN_BOLD = extract_BOLD_signal(ROI_DMN, ROI_mask_DMN, timeseries, total_blocks, all_time);
ROI_ts_DMN = DMN_BOLD.ROI_ts;
ROI_ts_DMN_avg = DMN_BOLD.ROI_ts_avg;
ROI_ts_DMN_mean = DMN_BOLD.ROI_ts_mean;
ROI_ts_DMN_std = DMN_BOLD.ROI_ts_std;

FPCN_BOLD = extract_BOLD_signal(ROI_FPCN, ROI_mask_FPCN, timeseries, total_blocks, all_time);
ROI_ts_FPCN = FPCN_BOLD.ROI_ts;
ROI_ts_FPCN_avg = FPCN_BOLD.ROI_ts_avg;
ROI_ts_FPCN_mean = FPCN_BOLD.ROI_ts_mean;
ROI_ts_FPCN_std = FPCN_BOLD.ROI_ts_std;

DAN_BOLD = extract_BOLD_signal(ROI_DAN, ROI_mask_DAN, timeseries, total_blocks, all_time);
ROI_ts_DAN = DAN_BOLD.ROI_ts;
ROI_ts_DAN_avg = DAN_BOLD.ROI_ts_avg;
ROI_ts_DAN_mean = DAN_BOLD.ROI_ts_mean;
ROI_ts_DAN_std = DAN_BOLD.ROI_ts_std;

SN_BOLD = extract_BOLD_signal(ROI_SN, ROI_mask_SN, timeseries, total_blocks, all_time);
ROI_ts_SN = SN_BOLD.ROI_ts;
ROI_ts_SN_avg = SN_BOLD.ROI_ts_avg;
ROI_ts_SN_mean = SN_BOLD.ROI_ts_mean;
ROI_ts_SN_std = SN_BOLD.ROI_ts_std;

LC_BOLD = LC_stuffz(LC_mask_left, LC_mask_right, timeseries, total_blocks, all_time, rs_blocks_only, rs_time_only);
LC_ts = LC_BOLD.LC_ts;
LC_ts_rs_blocks = LC_BOLD.LC_ts_rs_blocks;
LC_ts_mean = LC_BOLD.LC_ts_mean;
LC_ts_std = LC_BOLD.LC_ts_std;

ROI_ts_mean_all = [ROI_ts_DMN_mean, ROI_ts_FPCN_mean, ROI_ts_DAN_mean, ROI_ts_SN_mean];
ROI_ts_std_all = [ROI_ts_DMN_std, ROI_ts_FPCN_std, ROI_ts_DAN_std, ROI_ts_SN_std];

%Static Correlations between all ROIs within each network
DMN_static_con = static_correlation(ROI_DMN, ROI_ts_DMN, total_blocks, all_time, LC_ts, rs_blocks_only);
ROI_static_DMN_everything = DMN_static_con.ROI_static_everything;
ROI_static_DMN_everything_p = DMN_static_con.ROI_static_everything_p;
ROI_static_DMN_all = DMN_static_con.ROI_static_all;
ROI_static_DMN_all_p = DMN_static_con.ROI_static_all_p;
avg_ROI_static_DMN_all = DMN_static_con.avg_ROI_static_all;
DMN_LC_static_all = DMN_static_con.LC_static_all;
DMN_LC_static_all_p = DMN_static_con.LC_static_all_p;
spear_corr_DMN = DMN_static_con.spear_corr;
spear_corr_DMN_p = DMN_static_con.spear_corr_p;
avg_ROI_static_DMN_Fish_z = DMN_static_con.avg_ROI_static_Fish_z;

FPCN_static_con = static_correlation(ROI_FPCN, ROI_ts_FPCN, total_blocks, all_time, LC_ts, rs_blocks_only);
ROI_static_FPCN_everything = FPCN_static_con.ROI_static_everything;
ROI_static_FPCN_everything_p = FPCN_static_con.ROI_static_everything_p;
ROI_static_FPCN_all = FPCN_static_con.ROI_static_all;
ROI_static_FPCN_all_p = FPCN_static_con.ROI_static_all_p;
avg_ROI_static_FPCN_all = FPCN_static_con.avg_ROI_static_all;
FPCN_LC_static_all = FPCN_static_con.LC_static_all;
FPCN_LC_static_all_p = FPCN_static_con.LC_static_all_p;
spear_corr_FPCN = FPCN_static_con.spear_corr;
spear_corr_FPCN_p = FPCN_static_con.spear_corr_p;
avg_ROI_static_FPCN_Fish_z = FPCN_static_con.avg_ROI_static_Fish_z;

DAN_static_con = static_correlation(ROI_DAN, ROI_ts_DAN, total_blocks, all_time, LC_ts, rs_blocks_only);
ROI_static_DAN_everything = DAN_static_con.ROI_static_everything;
ROI_static_DAN_everything_p = DAN_static_con.ROI_static_everything_p;
ROI_static_DAN_all = DAN_static_con.ROI_static_all;
ROI_static_DAN_all_p = DAN_static_con.ROI_static_all_p;
avg_ROI_static_DAN_all = DAN_static_con.avg_ROI_static_all;
DAN_LC_static_all = DAN_static_con.LC_static_all;
DAN_LC_static_all_p = DAN_static_con.LC_static_all_p;
spear_corr_DAN = DAN_static_con.spear_corr;
spear_corr_DAN_p = DAN_static_con.spear_corr_p;
avg_ROI_static_DAN_Fish_z = DAN_static_con.avg_ROI_static_Fish_z;

SN_static_con = static_correlation(ROI_SN, ROI_ts_SN, total_blocks, all_time, LC_ts, rs_blocks_only);
ROI_static_SN_everything = SN_static_con.ROI_static_everything;
ROI_static_SN_everything_p = SN_static_con.ROI_static_everything_p;
ROI_static_SN_all = SN_static_con.ROI_static_all;
ROI_static_SN_all_p = SN_static_con.ROI_static_all_p;
avg_ROI_static_SN_all = SN_static_con.avg_ROI_static_all;
SN_LC_static_all = SN_static_con.LC_static_all;
SN_LC_static_all_p = SN_static_con.LC_static_all_p;
spear_corr_SN = SN_static_con.spear_corr;
spear_corr_SN_p = SN_static_con.spear_corr_p;
avg_ROI_static_SN_Fish_z = SN_static_con.avg_ROI_static_Fish_z;

mean_spear_all = [squeeze(mean(mean(spear_corr_DMN))), squeeze(mean(mean(spear_corr_FPCN))), squeeze(mean(mean(spear_corr_DAN))), squeeze(mean(mean(spear_corr_SN)))];
std_spear_all = [squeeze(mean(std(spear_corr_DMN))), squeeze(mean(std(spear_corr_FPCN))), squeeze(mean(std(spear_corr_DAN))), squeeze(mean(std(spear_corr_SN)))];

ROI_DMN = double(ROI_DMN);
ROI_FPCN = double(ROI_FPCN);
ROI_DAN = double(ROI_DAN);
ROI_SN = double(ROI_SN);
avg_ROI_static_all = transpose([avg_ROI_static_DMN_all; avg_ROI_static_FPCN_all; avg_ROI_static_DAN_all; avg_ROI_static_SN_all]);
avg_ROI_static_all_std_err = [repmat(std(avg_ROI_static_DMN_all)./sqrt(ROI_DMN),total_blocks,1), repmat(std(avg_ROI_static_FPCN_all)./sqrt(ROI_FPCN),total_blocks,1), repmat(std(avg_ROI_static_DAN_all)./sqrt(ROI_DAN),total_blocks,1), repmat(std(avg_ROI_static_SN_all)./sqrt(ROI_SN),total_blocks,1)];
avg_ROI_static_Fish_z = transpose([avg_ROI_static_DMN_Fish_z; avg_ROI_static_FPCN_Fish_z; avg_ROI_static_DAN_Fish_z; avg_ROI_static_SN_Fish_z]);
avg_ROI_static_Fish_z_std_err = [repmat(std(avg_ROI_static_DMN_Fish_z)./sqrt(ROI_DMN),total_blocks,1), repmat(std(avg_ROI_static_DMN_Fish_z)./sqrt(ROI_FPCN),total_blocks,1), repmat(std(avg_ROI_static_DMN_Fish_z)./sqrt(ROI_DAN),total_blocks,1), repmat(std(avg_ROI_static_DMN_Fish_z)./sqrt(ROI_SN),total_blocks,1)];


%Dynamic Connectivity analysis using a sliding time window approach where
%the time window is moved over 1 time point.  Also a dynamic connectivity
%analysis using nonoverlapping windows.
delta_t = 40;
delta_t_2 = 10;

%Create equal sized windows with the proper time points for each resting
%state epoch
nonoverlap_RS_0 = 1:delta_t_2:rs_epoch_0(1,end)+delta_t_2-1;
nonoverlap_RS_0(1,end) = nonoverlap_RS_0(1,end)-1;
nonoverlap_RS_1 = squeeze_session_1(1,end)+1:delta_t_2:rs_epoch_1(1,end)+delta_t_2-1;
nonoverlap_RS_1(1,end) = nonoverlap_RS_1(1,end)-1;
nonoverlap_RS_2 = squeeze_session_2(1,end)+1:delta_t_2:rs_epoch_2(1,end)+delta_t_2-1;
nonoverlap_RS_2(1,end) = nonoverlap_RS_2(1,end)-1;
nonoverlap_RS_3 = squeeze_session_3(1,end)+1:delta_t_2:rs_epoch_3(1,end)+delta_t_2-1;
nonoverlap_RS_3(1,end) = nonoverlap_RS_3(1,end)-1;
nonoverlap_RS_4 = squeeze_session_4(1,end)+1:delta_t_2:rs_epoch_4(1,end)+delta_t_2-1;
nonoverlap_RS_4(1,end) = nonoverlap_RS_4(1,end)-1;

nonoverlap_RS_all = cat(1,nonoverlap_RS_0, nonoverlap_RS_1, nonoverlap_RS_2, nonoverlap_RS_3, nonoverlap_RS_4)';

DMN_dynamic_con = dynamic_connectivity(ROI_DMN, ROI_ts_DMN, rs_blocks_only, rs_tr, delta_t, LC_ts_rs_blocks, nonoverlap_RS_all, rs_time_only);
ROI_dynamic_DMN_all_blocks = DMN_dynamic_con.ROI_dynamic_all_blocks;
avg_ROI_dynamic_DMN = DMN_dynamic_con.avg_ROI_dynamic;
DMN_LC_dynamic = DMN_dynamic_con.LC_dynamic;
ROI_dynamic_DMN_nonoverlap_all_blocks = DMN_dynamic_con.ROI_dynamic_nonoverlap_all_blocks;
ROI_dynamic_DMN_nonoverlap_all_blocks_p = DMN_dynamic_con.ROI_dynamic_nonoverlap_all_blocks_p;
avg_ROI_dynamic_DMN_nonoverlap = DMN_dynamic_con.avg_ROI_dynamic_nonoverlap;

FPCN_dynamic_con = dynamic_connectivity(ROI_FPCN, ROI_ts_FPCN, rs_blocks_only, rs_tr, delta_t, LC_ts_rs_blocks, nonoverlap_RS_all, rs_time_only);
ROI_dynamic_FPCN_all_blocks = FPCN_dynamic_con.ROI_dynamic_all_blocks;
avg_ROI_dynamic_FPCN = FPCN_dynamic_con.avg_ROI_dynamic;
FPCN_LC_dynamic = FPCN_dynamic_con.LC_dynamic;
ROI_dynamic_FPCN_nonoverlap_all_blocks = FPCN_dynamic_con.ROI_dynamic_nonoverlap_all_blocks;
ROI_dynamic_FPCN_nonoverlap_all_blocks_p = FPCN_dynamic_con.ROI_dynamic_nonoverlap_all_blocks_p;
avg_ROI_dynamic_FPCN_nonoverlap = FPCN_dynamic_con.avg_ROI_dynamic_nonoverlap;

DAN_dynamic_con = dynamic_connectivity(ROI_DAN, ROI_ts_DAN, rs_blocks_only, rs_tr, delta_t, LC_ts_rs_blocks, nonoverlap_RS_all, rs_time_only);
ROI_dynamic_DAN_all_blocks = DAN_dynamic_con.ROI_dynamic_all_blocks;
avg_ROI_dynamic_DAN = DAN_dynamic_con.avg_ROI_dynamic;
DAN_LC_dynamic = DAN_dynamic_con.LC_dynamic;
ROI_dynamic_DAN_nonoverlap_all_blocks = DAN_dynamic_con.ROI_dynamic_nonoverlap_all_blocks;
ROI_dynamic_DAN_nonoverlap_all_blocks_p = DAN_dynamic_con.ROI_dynamic_nonoverlap_all_blocks_p;
avg_ROI_dynamic_DAN_nonoverlap = DAN_dynamic_con.avg_ROI_dynamic_nonoverlap;

SN_dynamic_con = dynamic_connectivity(ROI_SN, ROI_ts_SN, rs_blocks_only, rs_tr, delta_t, LC_ts_rs_blocks, nonoverlap_RS_all, rs_time_only);
ROI_dynamic_SN_all_blocks = SN_dynamic_con.ROI_dynamic_all_blocks;
avg_ROI_dynamic_SN = SN_dynamic_con.avg_ROI_dynamic;
SN_LC_dynamic = SN_dynamic_con.LC_dynamic;
ROI_dynamic_SN_nonoverlap_all_blocks = SN_dynamic_con.ROI_dynamic_nonoverlap_all_blocks;
ROI_dynamic_SN_nonoverlap_all_blocks_p = SN_dynamic_con.ROI_dynamic_nonoverlap_all_blocks_p;
avg_ROI_dynamic_SN_nonoverlap = SN_dynamic_con.avg_ROI_dynamic_nonoverlap;

all_dynamic_var = var(avg_ROI_dynamic_DMN);
all_dynamic_var(:,:,2) = var(avg_ROI_dynamic_FPCN);
all_dynamic_var(:,:,3) = var(avg_ROI_dynamic_DAN);
all_dynamic_var(:,:,4) = var(avg_ROI_dynamic_SN);

all_dynamic_var_nonoverlap = var(avg_ROI_dynamic_DMN_nonoverlap);
all_dynamic_var_nonoverlap(:,:,2) = var(avg_ROI_dynamic_FPCN_nonoverlap);
all_dynamic_var_nonoverlap(:,:,3) = var(avg_ROI_dynamic_DAN_nonoverlap);
all_dynamic_var_nonoverlap(:,:,4) = var(avg_ROI_dynamic_SN_nonoverlap);

all_dynamic = avg_ROI_dynamic_DMN;
all_dynamic(:,:,2) = avg_ROI_dynamic_FPCN;
all_dynamic(:,:,3) = avg_ROI_dynamic_DAN;
all_dynamic(:,:,4) = avg_ROI_dynamic_SN;
abs_deriv_all_dynamic = abs(diff(all_dynamic));
length_abs = length(avg_ROI_dynamic_DMN);
abs_deriv_all_dynamic_std_err= std(abs_deriv_all_dynamic(:,:,1))./sqrt(length_abs);
abs_deriv_all_dynamic_std_err(:,:,2) = std(abs_deriv_all_dynamic(:,:,2))./sqrt(length_abs);
abs_deriv_all_dynamic_std_err(:,:,3) = std(abs_deriv_all_dynamic(:,:,3))./sqrt(length_abs);
abs_deriv_all_dynamic_std_err(:,:,4) = std(abs_deriv_all_dynamic(:,:,3))./sqrt(length_abs);


%% Plots

rs = [0 1 2 3 4];
all_blocks = [0 1 2 3 4 5 6 7 8];
num_net = 4;
ti_net = ['DMN', 'FPCN', 'DAN', 'SN'];
color_vec = [0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0; 0.25 0.25 0.25];

%Mean BOLD Signal for Each Block
figure
for i = 1:num_net
    subplot(1,num_net,i)
    scatter(all_blocks,ROI_ts_mean_all(:,i),100,'filled')
    hold on
    errorbar(all_blocks,ROI_ts_mean_all(:,i),ROI_ts_std_all(:,i),'LineWidth',1.5)
    xticks(all_blocks)
    xticklabels({'RS0', 'S1','RS1','S2','RS2','S3','RS3', 'S4','RS4'})  
end
suptitle('Mean BOLD Signal During Each Block')

%Static Correlation For Each Block
max_static = max(max(avg_ROI_static_all_std_err+avg_ROI_static_all));
min_static= min(min(avg_ROI_static_all-avg_ROI_static_all_std_err));

figure
for i = 1:num_net
    subplot(1,num_net,i)
    scatter(all_blocks,avg_ROI_static_all(:,i),100,'filled')
    hold on
    errorbar(all_blocks,avg_ROI_static_all(:,i),avg_ROI_static_all_std_err(:,i),'LineWidth',1.5)
    xticks(all_blocks)
    xticklabels({'RS0', 'S1','RS1','S2','RS2','S3','RS3', 'S4','RS4'})  
    ylim([min_static max_static])
end
suptitle('Static Correlations For Each Block')

%Static Correlation For Each Block After Fisher Z Transform
max_Fish = max(max(avg_ROI_static_Fish_z_std_err+avg_ROI_static_Fish_z));
min_Fish = min(min(avg_ROI_static_Fish_z-avg_ROI_static_Fish_z_std_err));

figure
for i = 1:num_net
    subplot(1,num_net,i)
    scatter(all_blocks,avg_ROI_static_Fish_z(:,i),100,'filled')
    hold on
    errorbar(all_blocks,avg_ROI_static_Fish_z(:,i),avg_ROI_static_Fish_z_std_err(:,i),'LineWidth',1.5)
    ylim([min_Fish max_Fish])
    xticks(all_blocks)
    xticklabels({'RS0', 'S1','RS1','S2','RS2','S3','RS3', 'S4','RS4'})  
end
suptitle('Static Correlations For Each Block After Fisher Z Transform')

%Static Correlation stuffz wrt RS 0

figure
for i = 1:rs_blocks_only
   subplot(2,3,i)
   ti = sprintf('RS Epoch %d', rs(i));
   max_DMN_spear = max(max(max(spear_corr_DMN)));
   min_DMN_spear = min(min(min(spear_corr_DMN)));
   imagesc(spear_corr_DMN(:,:,i), [min_DMN_spear max_DMN_spear])   
   title(ti)
end
suptitle('DMN Spearman Correlation with RS 0')

figure
for i = 1:rs_blocks_only
   subplot(2,3,i)
   ti = sprintf('RS Epoch %d', rs(i));
   max_FPCN_spear = max(max(max(spear_corr_FPCN)));
   min_FPCN_spear = min(min(min(spear_corr_FPCN)));
   imagesc(spear_corr_FPCN(:,:,i), [min_FPCN_spear max_FPCN_spear])   
   title(ti)
end
suptitle('FPCN Spearman Correlation with RS 0')

figure
for i = 1:rs_blocks_only
   subplot(2,3,i)
   ti = sprintf('RS Epoch %d', rs(i));
   max_DAN_spear = max(max(max(spear_corr_DAN)));
   min_DAN_spear = min(min(min(spear_corr_DAN)));
   imagesc(spear_corr_DAN(:,:,i), [min_DAN_spear max_DAN_spear])   
   title(ti)
end
suptitle('DAN Spearman Correlation with RS 0')

figure
for i = 1:rs_blocks_only
   subplot(2,3,i)
   ti = sprintf('RS Epoch %d', rs(i));
   max_SN_spear = max(max(max(spear_corr_SN)));
   min_SN_spear = min(min(min(spear_corr_SN)));
   imagesc(spear_corr_SN(:,:,i), [min_SN_spear max_SN_spear])   
   title(ti)
end
suptitle('SN Spearman Correlation with RS 0')

figure
for i = 1:num_net
    h(i) = scatter(rs,mean_spear_all(:,i,:),100,color_vec(i,:,:),'filled');
    hold on
    errorbar(rs,mean_spear_all(:,i),std_spear_all(:,i),'Color',color_vec(i,:,:))
end
legend([h(1), h(2), h(3), h(4)],'DMN','FPCN','DAN','SN','Location','best')
xticks(rs)
xticklabels({'RS0','RS1','RS2','RS3','RS4'})
suptitle('Mean Network Spearman Correlation with RS 0')

%Dynamic Correlation Scatter Plots Overlapping Windows
dynamic_length = length(DMN_LC_dynamic);
for i = 1:rs_blocks_only
    ti = sprintf('RS Epoch %d.', rs(i));
    
    figure(100)
    x_min = min(min(min(DMN_LC_dynamic)));
    x_max = max(max(max(DMN_LC_dynamic)));
    y_min = min(min(min(avg_ROI_dynamic_DMN)));
    y_max = max(max(max(avg_ROI_dynamic_DMN)));
    subplot(1,rs_blocks_only,i)
    scatter(DMN_LC_dynamic(:,i),avg_ROI_dynamic_DMN(1:dynamic_length,i),'filled','b')
    hold on 
    scatter(mean(DMN_LC_dynamic(:,i)),mean(avg_ROI_dynamic_DMN(1:dynamic_length,i)),200,'filled','s','r')
    xlim([x_min x_max])
    ylim([y_min y_max])
    xlabel('LC <-> DMN')
    ylabel('Within DMN')
    title(ti)
    
    figure(101)
    x_min = min(min(min(FPCN_LC_dynamic)));
    x_max = max(max(max(FPCN_LC_dynamic)));
    y_min = min(min(min(avg_ROI_dynamic_FPCN)));
    y_max = max(max(max(avg_ROI_dynamic_FPCN)));
    subplot(1,rs_blocks_only,i)
    scatter(FPCN_LC_dynamic(:,i),avg_ROI_dynamic_FPCN(1:dynamic_length,i),'filled','b')
    hold on 
    scatter(mean(FPCN_LC_dynamic(:,i)),mean(avg_ROI_dynamic_FPCN(1:dynamic_length,i)),200,'filled','s','r')
    xlim([x_min x_max])
    ylim([y_min y_max])
    xlabel('LC <-> FPCN')
    ylabel('Within FPCN')
    title(ti)
    
    figure(102)
    x_min = min(min(min(DAN_LC_dynamic)));
    x_max = max(max(max(DAN_LC_dynamic)));
    y_min = min(min(min(avg_ROI_dynamic_DAN)));
    y_max = max(max(max(avg_ROI_dynamic_DAN)));
    subplot(1,rs_blocks_only,i)
    scatter(DAN_LC_dynamic(:,i),avg_ROI_dynamic_DAN(1:dynamic_length,i),'filled','b')
    hold on 
    scatter(mean(DAN_LC_dynamic(:,i)),mean(avg_ROI_dynamic_DAN(1:dynamic_length,i)),200,'filled','s','r')
    xlim([x_min x_max])
    ylim([y_min y_max])
    xlabel('LC <-> DAN')
    ylabel('Within DAN')
    title(ti)
    
    figure(103)
    x_min = min(min(min(SN_LC_dynamic)));
    x_max = max(max(max(SN_LC_dynamic)));
    y_min = min(min(min(avg_ROI_dynamic_SN)));
    y_max = max(max(max(avg_ROI_dynamic_SN)));
    subplot(1,rs_blocks_only,i)
    scatter(SN_LC_dynamic(:,i),avg_ROI_dynamic_SN(1:dynamic_length,i),'filled','b')
    hold on 
    scatter(mean(SN_LC_dynamic(:,i)),mean(avg_ROI_dynamic_SN(1:dynamic_length,i)),200,'filled','s','r')
    xlim([x_min x_max])
    ylim([y_min y_max])
    xlabel('LC <-> SN')
    ylabel('Within SN')
    title(ti)
   
end

%Variance of Dynamic Connectivity Overlapping Windows
all_dynamic_var = var(avg_ROI_dynamic_DMN);
all_dynamic_var(:,:,2) = var(avg_ROI_dynamic_FPCN);
all_dynamic_var(:,:,3) = var(avg_ROI_dynamic_DAN);
all_dynamic_var(:,:,4) = var(avg_ROI_dynamic_SN);
max_var = max(max(max(all_dynamic_var)));
max_var_nonoverlap = max(max(max(all_dynamic_var_nonoverlap)));

figure
for i = 1:num_net
   subplot(1,num_net,i)
   bar(rs,all_dynamic_var(:,:,i),'FaceColor','b')
   xticks(rs)
   xticklabels({'RS0','RS1','RS2','RS3','RS4'})
   ylim([0 max_var])   
end
suptitle('Variance of Overlapping Window Dynamic Connectivity')

figure
for i = 1:num_net
   subplot(1,num_net,i)
   bar(rs,all_dynamic_var_nonoverlap(:,:,i),'FaceColor','b')
   xticks(rs)
   xticklabels({'RS0','RS1','RS2','RS3','RS4'})
   ylim([0 max_var_nonoverlap])   
end
suptitle('Variance of NonOverlapping Window Dynamic Connectivity')


%Raw LC BOLD Signal

yy = ones(length(squeeze_session_1),1)*1010;
figure
plot(squeeze(LC_ts))
hold on
plot(squeeze_session_1,yy)
hold on
plot(squeeze_session_2,yy)
hold on
plot(squeeze_session_3,yy)
hold on
plot(squeeze_session_4,yy)
xlim([0 660])
xlabel('TR')
ylabel('BOLD')
title('Raw LC BOLD Signal')


%Absolute Value Derivative
figure
for j = 1:num_net
    h(j) = scatter(rs,mean(abs_deriv_all_dynamic(:,:,j)),100,color_vec(j,:,:),'filled');
    hold on
    errorbar(rs,mean(abs_deriv_all_dynamic(:,:,j)),abs_deriv_all_dynamic_std_err(:,:,j),'LineWidth',2,'Color',color_vec(j,:,:))
end
legend([h(1), h(2), h(3), h(4)],'DMN','FPCN','DAN','SN','Location','best')
xticks(rs)
xticklabels({'RS0','RS1','RS2','RS3','RS4'})
suptitle('Abs. Val. of Deriv. Dyn. Conn.')



%Save everything
save('ROI_ts_DMN.mat','ROI_ts_DMN')
save('ROI_ts_FPCN.mat','ROI_ts_FPCN')
save('ROI_ts_DAN.mat','ROI_ts_DAN')
save('ROI_ts_SN.mat','ROI_ts_SN')
save('LC_ts.mat','LC_ts')

save('ROI_ts_DMN_avg.mat','ROI_ts_DMN_avg')
save('ROI_ts_FPCN_avg.mat','ROI_ts_FPCN_avg')
save('ROI_ts_DAN_avg.mat','ROI_ts_DAN_avg')
save('ROI_ts_SN_avg.mat','ROI_ts_SN_avg')

save('ROI_ts_DMN_mean.mat','ROI_ts_DMN_mean')
save('ROI_ts_FPCN_mean.mat','ROI_ts_FPCN_mean')
save('ROI_ts_DAN_mean.mat','ROI_ts_DAN_mean')
save('ROI_ts_SN_mean.mat','ROI_ts_SN_mean')
save('LC_ts_mean.mat','LC_ts_mean')

save('ROI_ts_DMN_std.mat','ROI_ts_DMN_std')
save('ROI_ts_FPCN_std.mat','ROI_ts_FPCN_std')
save('ROI_ts_DAN_std.mat','ROI_ts_DAN_std')
save('ROI_ts_SN_std.mat','ROI_ts_SN_std')
save('LC_ts_std.mat','LC_ts_std')

save('ROI_static_DMN_all.mat','ROI_static_DMN_all')
save('ROI_static_FPCN_all.mat','ROI_static_FPCN_all')
save('ROI_static_DAN_all.mat','ROI_static_DAN_all')
save('ROI_static_SN_all.mat','ROI_static_SN_all')

save('ROI_static_DMN_all_p.mat','ROI_static_DMN_all_p')
save('ROI_static_FPCN_all_p.mat','ROI_static_FPCN_all_p')
save('ROI_static_DAN_all_p.mat','ROI_static_DAN_all_p')
save('ROI_static_SN_all_p.mat','ROI_static_SN_all_p')

save('avg_ROI_static_DMN_all.mat','avg_ROI_static_DMN_all')
save('avg_ROI_static_FPCN_all.mat','avg_ROI_static_FPCN_all')
save('avg_ROI_static_DAN_all.mat','avg_ROI_static_DAN_all')
save('avg_ROI_static_SN_all.mat','avg_ROI_static_SN_all')

save('spear_corr_DMN.mat','spear_corr_DMN')
save('spear_corr_FPCN.mat','spear_corr_FPCN')
save('spear_corr_DAN.mat','spear_corr_DAN')
save('spear_corr_SN.mat','spear_corr_SN')

save('mean_spear_all.mat','mean_spear_all')
save('std_spear_all.mat','std_spear_all')

save('DMN_LC_static_all.mat','DMN_LC_static_all')
save('FPCN_LC_static_all.mat','FPCN_LC_static_all')
save('DAN_LC_static_all.mat','DAN_LC_static_all')
save('SN_LC_static_all.mat','SN_LC_static_all')

save('DMN_LC_static_all_p.mat','DMN_LC_static_all_p')
save('FPCN_LC_static_all_p.mat','FPCN_LC_static_all_p')
save('DAN_LC_static_all_p.mat','DAN_LC_static_all_p')
save('SN_LC_static_all_p.mat','SN_LC_static_all_p')

save('avg_ROI_dynamic_DMN.mat','avg_ROI_dynamic_DMN')
save('avg_ROI_dynamic_FPCN.mat','avg_ROI_dynamic_FPCN')
save('avg_ROI_dynamic_DAN.mat','avg_ROI_dynamic_DAN')
save('avg_ROI_dynamic_SN.mat','avg_ROI_dynamic_SN')

save('avg_ROI_dynamic_DMN_nonoverlap.mat','avg_ROI_dynamic_DMN_nonoverlap')
save('avg_ROI_dynamic_FPCN_nonoverlap.mat','avg_ROI_dynamic_FPCN_nonoverlap')
save('avg_ROI_dynamic_DAN_nonoverlap.mat','avg_ROI_dynamic_DAN_nonoverlap')
save('avg_ROI_dynamic_SN_nonoverlap.mat','avg_ROI_dynamic_SN_nonoverlap')

save('mag_squeeze.mat','mag_squeeze')

save('DMN_LC_dynamic.mat','DMN_LC_dynamic');
save('FPCN_LC_dynamic.mat','FPCN_LC_dynamic');
save('DAN_LC_dynamic.mat','DAN_LC_dynamic');
save('SN_LC_dynamic.mat','SN_LC_dynamic');

save('avg_ROI_static_DMN_all.mat','avg_ROI_static_DMN_all');
save('avg_ROI_static_FPCN_all.mat','avg_ROI_static_FPCN_all');
save('avg_ROI_static_DAN_all.mat','avg_ROI_static_DAN_all');
save('avg_ROI_static_SN_all.mat','avg_ROI_static_SN_all');

save('avg_ROI_static_all.mat','avg_ROI_static_all');
save('avg_ROI_static_all_std_err.mat','avg_ROI_static_all_std_err');

save('avg_ROI_static_Fish_z.mat','avg_ROI_static_Fish_z')
save('avg_ROI_static_Fish_z_std_err.mat','avg_ROI_static_Fish_z_std_err')

save('all_dynamic_var.mat','all_dynamic')
save('all_dynamic.mat','all_dynamic')
save('abs_deriv_all_dynamic.mat','abs_deriv_all_dynamic')
save('abs_deriv_all_dynamic_std_err.mat','abs_deriv_all_dynamic_std_err')








