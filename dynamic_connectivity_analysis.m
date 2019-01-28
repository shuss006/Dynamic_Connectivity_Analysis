%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
%                                                                         %
%                  Static and Dynamic Connectivity Analysis               %
%                           Author: Sana Hussain                          %
%                    University of California, Riverside                  %
%                                                                         %
%-------------------------------------------------------------------------%

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

%Commonly used phrases:
%"none" OR "no squeeze"      = first resting state epoch with no prior squeezing
%"during_201"                = first squeezing session
%"first_squeeze" OR "201"    = resting state epoch after the first squeeze
%"during_202"                = second squeezing session
%"second_squeeze" OR "202"   = resting state epoch after the second squeeze
%"during_203"                = third squeezing session
%"third_squeeze" OR "203"    = resting state epoch after the third squeeze

%The locus coeruleus (LC) is also of interest in this experiment so the BOLD
%signal using an LC mask is acquired and eye tracking data (pupillometry data)
%was also acquired as a means of checking the reliability of the BOLD signal. 

%The NIFTI toolbox is needed to run this code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

%Load the ROI masks and fMRI data

test = load_untouch_nii('default_mode_network_gopi_ROIs.nii');
ROI_mask_DMN = test.img;

test = load_untouch_nii('fronto_parietal_control_network_gopi_adjusted_ROIs.nii');
ROI_mask_FPCN = test.img;

test = load_untouch_nii('dorsal_attention_network_gopi_ROIs.nii');
ROI_mask_DAN = test.img;

test = load_untouch_nii('salience_network_ROIs.nii');
ROI_mask_SN = test.img;

test = load_untouch_nii('MNI_LC_mask_left_correct.nii');
LC_mask_left = test.img;

test = load_untouch_nii('MNI_LC_mask_right_correct.nii');
LC_mask_right = test.img;

test = load_untouch_nii('csf_mask_adjusted_bin.nii');
csf_mask = test.img;

test = load_untouch_nii('MNI_fourth_ventricle_mask_bin.nii');
fourth_ventricle_mask = test.img;

test = load_untouch_nii('data_in_MNI.nii');
timeseries = test.img;

[nx, ny, nz, nt] = size(timeseries); 
ROI_DMN = max(ROI_mask_DMN(:));
ROI_FPCN = max(ROI_mask_FPCN(:));
ROI_DAN = max(ROI_mask_DAN(:));
ROI_SN = max(ROI_mask_SN(:));

%%

%Extract the time series for each ROI from the masks loaded above.

%Start with the ROIs that make up DMN 
ROI_ts_DMN=zeros(nt,ROI_DMN);
for r = 1:ROI_DMN %Iterating over all ROIs
    count = 1;
    clear p;
    for x = 1:nx %Iterating over all values of x
        for y = 1:ny %Iterating over all values of y
            for z = 1:nz %Iterating over all values of z
                %Extract 1 ROI for the entire time series
                if(ROI_mask_DMN(x,y,z) == r) %Finding all of the ROIs in the time series
                    p(:,count) = squeeze(timeseries(x,y,z,:)); %Taking all the values in the time series for a each ROI
                    count = count+1;
                end
            end
        end
    end
    %Average the dummy variable, a
    for l = 1:nt
        ROI_ts_DMN(l,r) = mean(p(l,:));
    end
   clear p;
end

%Use the rest_filter function to bandwidth filter the DMN time series
for r = 1:ROI_DMN
    ROI_ts_DMN(:,r) = rest_filter(2,[0.01 0.1],ROI_ts_DMN(:,r));
end
save('ROI_ts_DMN.mat','ROI_ts_DMN');





%Repeat this analysis for the ROIs that make up FPCN
ROI_ts_FPCN=zeros(nt,ROI_FPCN);
for r = 1:ROI_FPCN %Iterating over all ROIs
    count = 1;
    clear p;
    for x = 1:nx %Iterating over all values of x
        for y = 1:ny %Iterating over all values of y
            for z = 1:nz %Iterating over all values of z
                %Extract 1 ROI for the entire time series
                if(ROI_mask_FPCN(x,y,z) == r) %Finding all of the ROIs in the time series
                    p(:,count) = squeeze(timeseries(x,y,z,:)); %Taking all the values in the time series for a each ROI
                    count = count+1;
                end
            end
        end
    end
    %Average the dummy variable, a
    for l = 1:nt
        ROI_ts_FPCN(l,r) = mean(p(l,:));
    end
   clear p;
end

for r = 1:ROI_FPCN
    ROI_ts_FPCN(:,r) = rest_filter(2,[0.01 0.1],ROI_ts_FPCN(:,r));
end
save('ROI_ts_FPCN.mat','ROI_ts_FPCN');




%Repeat this analysis for the ROIs that make up DAN
ROI_ts_DAN=zeros(nt,ROI_DAN);
for r = 1:ROI_DAN %Iterating over all ROIs
    count = 1;
    clear p;
    for x = 1:nx %Iterating over all values of x
        for y = 1:ny %Iterating over all values of y
            for z = 1:nz %Iterating over all values of z
                %Extract 1 ROI for the entire time series
                if(ROI_mask_DAN(x,y,z) == r) %Finding all of the ROIs in the time series
                    p(:,count) = squeeze(timeseries(x,y,z,:)); %Taking all the values in the time series for a each ROI
                    count = count+1;
                end
            end
        end
    end
    %Average the dummy variable, a
    for l = 1:nt
        ROI_ts_DAN(l,r) = mean(p(l,:));
    end
   clear p;
end

for r = 1:ROI_DAN
    ROI_ts_DAN(:,r) = rest_filter(2,[0.01 0.1],ROI_ts_DAN(:,r));
end
save('ROI_ts_DAN.mat','ROI_ts_DAN');





%Repeat this analysis for the ROIs that make up SN
ROI_ts_SN=zeros(nt,ROI_SN);
for r = 1:ROI_SN %Iterating over all ROIs
    count = 1;
    clear p;
    for x = 1:nx %Iterating over all values of x
        for y = 1:ny %Iterating over all values of y
            for z = 1:nz %Iterating over all values of z
                %Extract 1 ROI for the entire time series
                if(ROI_mask_SN(x,y,z) == r) %Finding all of the ROIs in the time series
                    p(:,count) = squeeze(timeseries(x,y,z,:)); %Taking all the values in the time series for a each ROI
                    count = count+1;
                end
            end
        end
    end
    %Average the dummy variable, a
    for l = 1:nt
        ROI_ts_SN(l,r) = mean(p(l,:));
    end
   clear p;
end

for r = 1:ROI_SN
    ROI_ts_SN(:,r) = rest_filter(2,[0.01 0.1],ROI_ts_SN(:,r));
end
save('ROI_ts_SN.mat','ROI_ts_SN');





%Repeat this analysis for the LC ROI mask
for r = 1:1
    clear b;
    clear f;
    for x = 1:nx
        count = 1;
        for y = 1:ny
            for z = 1:nz
                if(LC_mask_left(x,y,z) == 1)
                    b(:, count) = squeeze(timeseries(x,y,z,:));
                    count = count+1;
                end
                if(LC_mask_right(x,y,z) == 1)
                    f(:,count) = squeeze(timeseries(x,y,z,:));
                    count = count+1;
                end
            end
        end
    end
    %Average the time series
    for l = 1:nt
        LC_left_ts(l,r) = mean(b(l,:));
        LC_right_ts(l,r) = mean(f(l,:));
    end
end

LC_left_ts(:,1) = rest_filter(2,[0.01 0.1],LC_left_ts(:,1));
LC_right_ts(:,1) = rest_filter(2,[0.01 0.1],LC_right_ts(:,1));

LC_ts = (LC_left_ts+LC_right_ts)./2;
save('LC_ts.mat','LC_ts');

%Now divide each network time series into each of the 90 time point resting
%state epoch

LC_ts_none = LC_ts(1:90,1);
LC_ts_during_201 = LC_ts(91:105,1);
LC_ts_201 = LC_ts(106:195,1);
LC_ts_during_202 = LC_ts(196:210,1);
LC_ts_202 = LC_ts(211:300,1);
LC_ts_during_203 = LC_ts(301:315,1);
LC_ts_203 = LC_ts(316:end,1);

save('LC_ts_none.mat','LC_ts_none')
save('LC_ts_201.mat','LC_ts_201')
save('LC_ts_202.mat','LC_ts_202')
save('LC_ts_203.mat','LC_ts_203')
save('LC_ts.mat','LC_ts');
save('z_LC_ts.mat', 'z_LC_ts');

ROI_ts_size_DMN_none = ROI_ts_DMN(1:90,:);
ROI_ts_size_DMN_during201 = ROI_ts_DMN(91:105,:);
ROI_ts_size_DMN_201 = ROI_ts_DMN(106:195,:);
ROI_ts_size_DMN_during202 = ROI_ts_DMN(196:210,:);
ROI_ts_size_DMN_202 = ROI_ts_DMN(211:300,:);
ROI_ts_size_DMN_during203 = ROI_ts_DMN(301:315,:);
ROI_ts_size_DMN_203 = ROI_ts_DMN(316:end,:);


ROI_ts_size_FPCN_none = ROI_ts_FPCN(1:90,:);
ROI_ts_size_FPCN_during201 = ROI_ts_FPCN(91:105,:);
ROI_ts_size_FPCN_201 = ROI_ts_FPCN(106:195,:);
ROI_ts_size_FPCN_during202 = ROI_ts_FPCN(196:210,:);
ROI_ts_size_FPCN_202 = ROI_ts_FPCN(211:300,:);
ROI_ts_size_FPCN_during203 = ROI_ts_FPCN(301:315,:);
ROI_ts_size_FPCN_203 = ROI_ts_FPCN(316:end,:);

ROI_ts_size_DAN_none = ROI_ts_DAN(1:90,:);
ROI_ts_size_DAN_during201 = ROI_ts_DAN(91:105,:);
ROI_ts_size_DAN_201 = ROI_ts_DAN(106:195,:);
ROI_ts_size_DAN_during202 = ROI_ts_DAN(196:210,:);
ROI_ts_size_DAN_202 = ROI_ts_DAN(211:300,:);
ROI_ts_size_DAN_during203 = ROI_ts_DAN(301:315,:);
ROI_ts_size_DAN_203 = ROI_ts_DAN(316:end,:);


ROI_ts_size_SN_none = ROI_ts_SN(1:90,:);
ROI_ts_size_SN_during201 = ROI_ts_SN(91:105,:);
ROI_ts_size_SN_201 = ROI_ts_SN(106:195,:);
ROI_ts_size_SN_during202 = ROI_ts_SN(196:210,:);
ROI_ts_size_SN_202 = ROI_ts_SN(211:300,:);
ROI_ts_size_SN_during203 = ROI_ts_SN(301:315,:);
ROI_ts_size_SN_203 = ROI_ts_SN(316:end,:);

%%

%Average the ROI timeseries to create a single time series representing the
%overall network during the resting state epochs and squeezing sessions

%DMN Resting State Epochs
for k = 1:90
    count = 1;
    for r1 = 1:ROI_DMN
            output_none(count,1) = squeeze(ROI_ts_size_DMN_none(k,r1));
            output_201(count,1) = squeeze(ROI_ts_size_DMN_201(k,r1));
            output_202(count,1) = squeeze(ROI_ts_size_DMN_202(k,r1));
            output_203(count,1) = squeeze(ROI_ts_size_DMN_203(k,r1));
            count = count+1;
    end
    ROI_ts_DMN_none(k,1) = mean(output_none(:,1));
    ROI_ts_DMN_201(k,1) = mean(output_201(:,1));
    ROI_ts_DMN_202(k,1) = mean(output_202(:,1));
    ROI_ts_DMN_203(k,1) = mean(output_203(:,1));
    
    clear output_none
    clear output_201
    clear output_202
    clear output_203
    clear count
end

ROI_ts_DMN_all = [ROI_ts_DMN_none, ROI_ts_DMN_201, ROI_ts_DMN_202, ROI_ts_DMN_203];
save('ROI_ts_DMN_all.mat','ROI_ts_DMN_all')

%DMN Squeezing Sessions
for k = 1:15
    count = 1;
    for r1 = 1:ROI_DMN
            output_201(count,1) = squeeze(ROI_ts_size_DMN_during201(k,r1));
            output_202(count,1) = squeeze(ROI_ts_size_DMN_during202(k,r1));
            output_203(count,1) = squeeze(ROI_ts_size_DMN_during203(k,r1));
            count = count+1;
    end
    ROI_ts_DMN_during_201(k,1) = mean(output_201(:,1));
    ROI_ts_DMN_during_202(k,1) = mean(output_202(:,1));
    ROI_ts_DMN_during_203(k,1) = mean(output_203(:,1));
    clear output_201
    clear output_202
    clear output_203
    clear count
end

%FPCN Resting State Epochs
for k = 1:90
    count = 1;
    for r1 = 1:ROI_FPCN
            output_none(count,1) = squeeze(ROI_ts_size_FPCN_none(k,r1));
            output_201(count,1) = squeeze(ROI_ts_size_FPCN_201(k,r1));
            output_202(count,1) = squeeze(ROI_ts_size_FPCN_202(k,r1));
            output_203(count,1) = squeeze(ROI_ts_size_FPCN_203(k,r1));
            count = count+1;
    end
    ROI_ts_FPCN_none(k,1) = mean(output_none(:,1));
    ROI_ts_FPCN_201(k,1) = mean(output_201(:,1));
    ROI_ts_FPCN_202(k,1) = mean(output_202(:,1));
    ROI_ts_FPCN_203(k,1) = mean(output_203(:,1));
    
    clear output_none
    clear output_201
    clear output_202
    clear output_203
    clear count
end

%FPCN Squeezing Sessions
for k = 1:15
    count = 1;
    for r1 = 1:ROI_FPCN
            output_201(count,1) = squeeze(ROI_ts_size_FPCN_during201(k,r1));
            output_202(count,1) = squeeze(ROI_ts_size_FPCN_during202(k,r1));
            output_203(count,1) = squeeze(ROI_ts_size_FPCN_during203(k,r1));
            count = count+1;
    end
    ROI_ts_FPCN_during_201(k,1) = mean(output_201(:,1));
    ROI_ts_FPCN_during_202(k,1) = mean(output_202(:,1));
    ROI_ts_FPCN_during_203(k,1) = mean(output_203(:,1));
    clear output_201
    clear output_202
    clear output_203
    clear count
end

%DAN Resting State Epochs
for k = 1:90
    count = 1;
    for r1 = 1:ROI_DAN
            output_none(count,1) = squeeze(ROI_ts_size_DAN_none(k,r1));
            output_201(count,1) = squeeze(ROI_ts_size_DAN_201(k,r1));
            output_202(count,1) = squeeze(ROI_ts_size_DAN_202(k,r1));
            output_203(count,1) = squeeze(ROI_ts_size_DAN_203(k,r1));
            count = count+1;
    end
    ROI_ts_DAN_none(k,1) = mean(output_none(:,1));
    ROI_ts_DAN_201(k,1) = mean(output_201(:,1));
    ROI_ts_DAN_202(k,1) = mean(output_202(:,1));
    ROI_ts_DAN_203(k,1) = mean(output_203(:,1));
    
    clear output_none
    clear output_201
    clear output_202
    clear output_203
    clear count
end
ROI_ts_DAN_all = [ROI_ts_DAN_none, ROI_ts_DAN_201, ROI_ts_DAN_202, ROI_ts_DAN_203];
save('ROI_ts_DAN_all.mat','ROI_ts_DAN_all')

%DAN Squeezing Sessions
for k = 1:15
    count = 1;
    for r1 = 1:ROI_DAN
            output_201(count,1) = squeeze(ROI_ts_size_DAN_during201(k,r1));
            output_202(count,1) = squeeze(ROI_ts_size_DAN_during202(k,r1));
            output_203(count,1) = squeeze(ROI_ts_size_DAN_during203(k,r1));
            count = count+1;
    end
    ROI_ts_DAN_during_201(k,1) = mean(output_201(:,1));
    ROI_ts_DAN_during_202(k,1) = mean(output_202(:,1));
    ROI_ts_DAN_during_203(k,1) = mean(output_203(:,1));
    clear output_201
    clear output_202
    clear output_203
    clear count
end

%SN Resting State Epochs
for k = 1:90
    count = 1;
    for r1 = 1:ROI_SN
            output_none(count,1) = squeeze(ROI_ts_size_SN_none(k,r1));
            output_201(count,1) = squeeze(ROI_ts_size_SN_201(k,r1));
            output_202(count,1) = squeeze(ROI_ts_size_SN_202(k,r1));
            output_203(count,1) = squeeze(ROI_ts_size_SN_203(k,r1));
            count = count+1;
    end
    ROI_ts_SN_none(k,1) = mean(output_none(:,1));
    ROI_ts_SN_201(k,1) = mean(output_201(:,1));
    ROI_ts_SN_202(k,1) = mean(output_202(:,1));
    ROI_ts_SN_203(k,1) = mean(output_203(:,1));
    
    clear output_none
    clear output_201
    clear output_202
    clear output_203
    clear count
end
ROI_ts_SN_all = [ROI_ts_SN_none, ROI_ts_SN_201, ROI_ts_SN_202, ROI_ts_SN_203];
save('ROI_ts_SN_all.mat','ROI_ts_SN_all')

%SN Squeezing Sessions
for k = 1:15
    count = 1;
    for r1 = 1:ROI_SN
            output_201(count,1) = squeeze(ROI_ts_size_SN_during201(k,r1));
            output_202(count,1) = squeeze(ROI_ts_size_SN_during202(k,r1));
            output_203(count,1) = squeeze(ROI_ts_size_SN_during203(k,r1));
            count = count+1;
    end
    ROI_ts_SN_during_201(k,1) = mean(output_201(:,1));
    ROI_ts_SN_during_202(k,1) = mean(output_202(:,1));
    ROI_ts_SN_during_203(k,1) = mean(output_203(:,1));
    clear output_201
    clear output_202
    clear output_203
    clear count
end

%%

%Compute the static correlations and the corresponding p-values between
%each ROI in each network.

%DMN Static Correlations
for r1 = 1:ROI_DMN
    for r2 = 1:ROI_DMN
        [oo,p_oo] = corrcoef(squeeze(ROI_ts_DMN(1:90,r1)),squeeze(ROI_ts_DMN(1:90,r2)));
        ROI_static_DMN_none(r1,r2) = oo(2,1);  
        ROI_static_DMN_none_p(r1,r2) = p_oo(2,1);
        [vvv,p_vvv] = corrcoef(squeeze(ROI_ts_DMN(91:105,r1)),squeeze(ROI_ts_DMN(91:105,r2)));
        ROI_static_DMN_during_201(r1,r2) = vvv(2,1);
        ROI_static_DMN_during_201_p(r1,r2) = p_vvv(2,1);
        [vv, p_vv] = corrcoef(squeeze(ROI_ts_DMN(106:195,r1)),squeeze(ROI_ts_DMN(106:195,r2)));
        ROI_static_DMN_201(r1,r2) = vv(2,1);      
        ROI_static_DMN_201_p(r1,r2) = p_vv(2,1);
        [www, p_www] = corrcoef(squeeze(ROI_ts_DMN(196:210,r1)),squeeze(ROI_ts_DMN(196:210,r2)));
        ROI_static_DMN_during_202(r1,r2) = www(2,1);
        ROI_static_DMN_during_202_p(r1,r2) = p_www(2,1);
        [ww, p_ww] = corrcoef(squeeze(ROI_ts_DMN(211:300,r1)),squeeze(ROI_ts_DMN(211:300,r2)));
        ROI_static_DMN_202(r1,r2) = ww(2,1);  
        ROI_static_DMN_202_p(r1,r2)  = p_ww(2,1);
        [zzz, p_zzz] = corrcoef(squeeze(ROI_ts_DMN(301:315,r1)),squeeze(ROI_ts_DMN(301:315,r2)));
        ROI_static_DMN_during_203(r1,r2) = zzz(2,1);
        ROI_static_DMN_during_203_p(r1,r2) = p_zzz(2,1);
        [zz, p_zz] = corrcoef(squeeze(ROI_ts_DMN(316:end,r1)),squeeze(ROI_ts_DMN(316:end,r2)));
        ROI_static_DMN_203(r1,r2) = zz(2,1);
        ROI_static_DMN_203_p(r1,r2) = p_zz(2,1);
    end
end

%Placing all the DMN Static Correlations Values into a single array
count = 1;
for r1 = 1:ROI_DMN
    for r2 = 1:ROI_DMN
        if r2>r1
        all_ROI_static_DMN_none(count,1) = ROI_static_DMN_none(r1,r2);
        all_ROI_static_DMN_none_p(count,1) = ROI_static_DMN_none_p(r1,r2);
        all_ROI_static_DMN_during_201(count,1) = ROI_static_DMN_during_201(r1,r2);
        all_ROI_static_DMN_during_201_p(count,1) = ROI_static_DMN_during_201_p(r1,r2);
        all_ROI_static_DMN_201(count,1) = ROI_static_DMN_201(r1,r2);
        all_ROI_static_DMN_201_p(count,1) = ROI_static_DMN_201_p(r1,r2);
        all_ROI_static_DMN_during_202(count,1) = ROI_static_DMN_during_202(r1,r2);
        all_ROI_static_DMN_during_202_p(count,1) = ROI_static_DMN_during_202_p(r1,r2);
        all_ROI_static_DMN_202(count,1) = ROI_static_DMN_202(r1,r2);
        all_ROI_static_DMN_202_p(count,1) = ROI_static_DMN_202_p(r1,r2);
        all_ROI_static_DMN_during_203(count,1) = ROI_static_DMN_during_203(r1,r2);
        all_ROI_static_DMN_during_203_p(count,1) = ROI_static_DMN_during_203_p(r1,r2);
        all_ROI_static_DMN_203(count,1) = ROI_static_DMN_203(r1,r2);
        all_ROI_static_DMN_203_p(count,1) = ROI_static_DMN_203_p(r1,r2);
        
        count = count + 1;
        end
    end
end

%Correlating each static correlation values with those from the "no
%squeeze" condition (i.e. the first resting state epoch)
[DMN_R_none_none, DMN_p_none_none, DMN_rl_none_none, DMN_ru_none_none] = corrcoef(ROI_static_DMN_none,ROI_static_DMN_none);
[DMN_R_none_201, DMN_p_none_201, DMN_rl_none_201, DMN_ru_none_201] = corrcoef(ROI_static_DMN_none,ROI_static_DMN_201);
[DMN_R_none_202, DMN_p_none_202, DMN_rl_none_202, DMN_ru_none_202] = corrcoef(ROI_static_DMN_none,ROI_static_DMN_202);
[DMN_R_none_203, DMN_p_none_203, DMN_rl_none_203, DMN_ru_none_203] = corrcoef(ROI_static_DMN_none,ROI_static_DMN_203);


%FPCN Static Correlations
for r1 = 1:ROI_FPCN
    for r2 = 1:ROI_FPCN
        oo = corrcoef(squeeze(ROI_ts_FPCN(1:90,r1)),squeeze(ROI_ts_FPCN(1:90,r2)));
        ROI_static_FPCN_none(r1,r2) = oo(2,1);   
        vvv = corrcoef(squeeze(ROI_ts_FPCN(91:105,r1)),squeeze(ROI_ts_FPCN(91:105,r2)));
        ROI_static_FPCN_during_201(r1,r2) = vvv(2,1);
        vv = corrcoef(squeeze(ROI_ts_FPCN(106:195,r1)),squeeze(ROI_ts_FPCN(106:195,r2)));
        ROI_static_FPCN_201(r1,r2) = vv(2,1);      
        www = corrcoef(squeeze(ROI_ts_FPCN(196:210,r1)),squeeze(ROI_ts_FPCN(196:210,r2)));
        ROI_static_FPCN_during_202(r1,r2) = www(2,1);
        ww = corrcoef(squeeze(ROI_ts_FPCN(211:300,r1)),squeeze(ROI_ts_FPCN(211:300,r2)));
        ROI_static_FPCN_202(r1,r2) = ww(2,1);  
        zzz = corrcoef(squeeze(ROI_ts_FPCN(301:315,r1)),squeeze(ROI_ts_FPCN(301:315,r2)));
        ROI_static_FPCN_during_203(r1,r2) = zzz(2,1);
        zz = corrcoef(squeeze(ROI_ts_FPCN(316:end,r1)),squeeze(ROI_ts_FPCN(316:end,r2)));
        ROI_static_FPCN_203(r1,r2) = zz(2,1);
    end
end

%Placing all the FPCN Static Correlations Values into a single array
count = 1;
for r1 = 1:ROI_FPCN
    for r2 = 1:ROI_FPCN
        if r2>r1
        all_ROI_static_FPCN_none(count,1) = ROI_static_FPCN_none(r1,r2);
        all_ROI_static_FPCN_during_201(count,1) = ROI_static_FPCN_during_201(r1,r2);
        all_ROI_static_FPCN_201(count,1) = ROI_static_FPCN_201(r1,r2);
        all_ROI_static_FPCN_during_202(count,1) = ROI_static_FPCN_during_202(r1,r2);
        all_ROI_static_FPCN_202(count,1) = ROI_static_FPCN_202(r1,r2);
        all_ROI_static_FPCN_during_203(count,1) = ROI_static_FPCN_during_203(r1,r2);
        all_ROI_static_FPCN_203(count,1) = ROI_static_FPCN_203(r1,r2);
        
        count = count + 1;
        end
    end
end


%DAN Static Correlations
for r1 = 1:ROI_DAN
    for r2 = 1:ROI_DAN
        [oo, p_oo] = corrcoef(squeeze(ROI_ts_DAN(1:90,r1)),squeeze(ROI_ts_DAN(1:90,r2)));
        ROI_static_DAN_none(r1,r2) = oo(2,1);   
        ROI_static_DAN_none_p(r1,r2) = p_oo(2,1);
        [vvv, p_vvv] = corrcoef(squeeze(ROI_ts_DAN(91:105,r1)),squeeze(ROI_ts_DAN(91:105,r2)));
        ROI_static_DAN_during_201(r1,r2) = vvv(2,1);
        ROI_static_DAN_during_201_p(r1,r2) = p_vvv(2,1);
        [vv, p_vv] = corrcoef(squeeze(ROI_ts_DAN(106:195,r1)),squeeze(ROI_ts_DAN(106:195,r2)));
        ROI_static_DAN_201(r1,r2) = vv(2,1);    
        ROI_static_DAN_201_p(r1,r2) = p_vv(2,1);
        [www, p_www] = corrcoef(squeeze(ROI_ts_DAN(196:210,r1)),squeeze(ROI_ts_DAN(196:210,r2)));
        ROI_static_DAN_during_202(r1,r2) = www(2,1);
        ROI_static_DAN_during_202_p(r1,r2) = p_www(2,1);
        [ww, p_ww] = corrcoef(squeeze(ROI_ts_DAN(211:300,r1)),squeeze(ROI_ts_DAN(211:300,r2)));
        ROI_static_DAN_202(r1,r2) = ww(2,1);  
        ROI_static_DAN_202_p(r1,r2) = p_ww(2,1);
        [zzz, p_zzz] = corrcoef(squeeze(ROI_ts_DAN(301:315,r1)),squeeze(ROI_ts_DAN(301:315,r2)));
        ROI_static_DAN_during_203(r1,r2) = zzz(2,1);
        ROI_static_DAN_during_203_p(r1,r2) = p_zzz(2,1);
        [zz, p_zz] = corrcoef(squeeze(ROI_ts_DAN(316:end,r1)),squeeze(ROI_ts_DAN(316:end,r2)));
        ROI_static_DAN_203(r1,r2) = zz(2,1);
        ROI_static_DAN_203_p(r1,r2) = p_zz(2,1);
    end
end

%Placing all the DAN Static Correlations Values into a single array
count = 1;
for r1 = 1:ROI_DAN
    for r2 = 1:ROI_DAN
        if r2>r1
        all_ROI_static_DAN_none(count,1) = ROI_static_DAN_none(r1,r2);
        all_ROI_static_DAN_none_p(count,1) = ROI_static_DAN_none_p(r1,r2);
        all_ROI_static_DAN_during_201(count,1) = ROI_static_DAN_during_201(r1,r2);
        all_ROI_static_DAN_during_201_p(count,1) = ROI_static_DAN_during_201_p(r1,r2);
        all_ROI_static_DAN_201(count,1) = ROI_static_DAN_201(r1,r2);
        all_ROI_static_DAN_201_p(count,1) = ROI_static_DAN_201_p(r1,r2);
        all_ROI_static_DAN_during_202(count,1) = ROI_static_DAN_during_202(r1,r2);
        all_ROI_static_DAN_during_202_p(count,1) = ROI_static_DAN_during_202_p(r1,r2);
        all_ROI_static_DAN_202(count,1) = ROI_static_DAN_202(r1,r2);
        all_ROI_static_DAN_202_p(count,1) = ROI_static_DAN_202_p(r1,r2);
        all_ROI_static_DAN_during_203(count,1) = ROI_static_DAN_during_203(r1,r2);
        all_ROI_static_DAN_during_203_p(count,1) = ROI_static_DAN_during_203_p(r1,r2);
        all_ROI_static_DAN_203(count,1) = ROI_static_DAN_203(r1,r2);
        all_ROI_static_DAN_203_p(count,1) = ROI_static_DAN_203_p(r1,r2);
        count = count + 1;
        end
    end
end

%Correlating each static correlation values with those from the "no
%squeeze" condition (i.e. the first resting state epoch)
[DAN_R_none_none, DAN_p_none_none, DAN_rl_none_none, DAN_ru_none_none] = corrcoef(ROI_static_DAN_none,ROI_static_DAN_none);
[DAN_R_none_201, DAN_p_none_201, DAN_rl_none_201, DAN_ru_none_201] = corrcoef(ROI_static_DAN_none,ROI_static_DAN_201);
[DAN_R_none_202, DAN_p_none_202, DAN_rl_none_202, DAN_ru_none_202] = corrcoef(ROI_static_DAN_none,ROI_static_DAN_202);
[DAN_R_none_203, DAN_p_none_203, DAN_rl_none_203, DAN_ru_none_203] = corrcoef(ROI_static_DAN_none,ROI_static_DAN_203);

%SN Static Correlations
for r1 = 1:ROI_SN
    for r2 = 1:ROI_SN
        [oo, p_oo] = corrcoef(squeeze(ROI_ts_SN(1:90,r1)),squeeze(ROI_ts_SN(1:90,r2)));
        ROI_static_SN_none(r1,r2) = oo(2,1);   
        ROI_static_SN_none_p(r1,r2) = p_oo(2,1);
        [vvv, p_vvv] = corrcoef(squeeze(ROI_ts_SN(91:105,r1)),squeeze(ROI_ts_SN(91:105,r2)));
        ROI_static_SN_during_201(r1,r2) = vvv(2,1);
        ROI_static_SN_during_201_p(r1,r2) = p_vvv(2,1);
        [vv, p_vv] = corrcoef(squeeze(ROI_ts_SN(106:195,r1)),squeeze(ROI_ts_SN(106:195,r2)));
        ROI_static_SN_201(r1,r2) = vv(2,1);      
        ROI_static_SN_201_p(r1,r2) = p_vv(2,1);
        [www, p_www] = corrcoef(squeeze(ROI_ts_SN(196:210,r1)),squeeze(ROI_ts_SN(196:210,r2)));
        ROI_static_SN_during_202(r1,r2) = www(2,1);
        ROI_static_SN_during_202_p(r1,r2) = p_www(2,1);
        [ww, p_www] = corrcoef(squeeze(ROI_ts_SN(211:300,r1)),squeeze(ROI_ts_SN(211:300,r2)));
        ROI_static_SN_202(r1,r2) = ww(2,1);  
        ROI_static_SN_202_p(r1,r2) = p_www(2,1);
        [zzz, p_zzz] = corrcoef(squeeze(ROI_ts_SN(301:315,r1)),squeeze(ROI_ts_SN(301:315,r2)));
        ROI_static_SN_during_203(r1,r2) = zzz(2,1);
        ROI_static_SN_during_203_p(r1,r2) = p_zzz(2,1);
        [zz, p_zz] = corrcoef(squeeze(ROI_ts_SN(316:end,r1)),squeeze(ROI_ts_SN(316:end,r2)));
        ROI_static_SN_203(r1,r2) = zz(2,1);
        ROI_static_SN_203_p(r1,r2) = p_zz(2,1);
    end
end

%Placing all the SN Static Correlations Values into a single array
count = 1;
for r1 = 1:ROI_SN
    for r2 = 1:ROI_SN
        if r2>r1
        all_ROI_static_SN_none(count,1) = ROI_static_SN_none(r1,r2);
        all_ROI_static_SN_none_p(count,1) = ROI_static_SN_none_p(r1,r2);
        all_ROI_static_SN_during_201(count,1) = ROI_static_SN_during_201(r1,r2);
        all_ROI_static_SN_during_201_p(count,1) = ROI_static_SN_during_201_p(r1,r2);
        all_ROI_static_SN_201(count,1) = ROI_static_SN_201(r1,r2);
        all_ROI_static_SN_201_p(count,1) = ROI_static_SN_201_p(r1,r2);
        all_ROI_static_SN_during_202(count,1) = ROI_static_SN_during_202(r1,r2);
        all_ROI_static_SN_during_202_p(count,1) = ROI_static_SN_during_202_p(r1,r2);
        all_ROI_static_SN_202(count,1) = ROI_static_SN_202(r1,r2);
        all_ROI_static_SN_202_p(count,1) = ROI_static_SN_202_p(r1,r2);
        all_ROI_static_SN_during_203(count,1) = ROI_static_SN_during_203(r1,r2);
        all_ROI_static_SN_during_203_p(count,1) = ROI_static_SN_during_203_p(r1,r2);
        all_ROI_static_SN_203(count,1) = ROI_static_SN_203(r1,r2);
        all_ROI_static_SN_203_p(count,1) = ROI_static_SN_203_p(r1,r2);
        
        count = count + 1;
        end
    end
end

%Correlating each static correlation values with those from the "no
%squeeze" condition (i.e. the first resting state epoch)
[SN_R_none_none, SN_p_none_none, SN_rl_none_none, SN_ru_none_none] = corrcoef(ROI_static_SN_none,ROI_static_SN_none);
[SN_R_none_201, SN_p_none_201, SN_rl_none_201, SN_ru_none_201] = corrcoef(ROI_static_SN_none,ROI_static_SN_201);
[SN_R_none_202, SN_p_none_202, SN_rl_none, SN_ru_none_202] = corrcoef(ROI_static_SN_none,ROI_static_SN_202);
[SN_R_none_203, SN_p_none_203, SN_rl_none_203, SN_ru_none_203] = corrcoef(ROI_static_SN_none,ROI_static_SN_203);

%%

%Compute the static correlations between the LC BOLD signal and each of the
%networks' BOLD signals

%DMN statically correlated with LC during each resting state epoch
count = 1;
for r1 = 1:ROI_DMN
    yy = corrcoef(squeeze(ROI_ts_DMN(1:90,r1)),squeeze(LC_ts(1:90,1)));
    all_LC_static_DMN_none(count,1) = yy(2,1);    
    hh = corrcoef(squeeze(ROI_ts_DMN(106:195,r1)),squeeze(LC_ts(106:195,1)));
    all_LC_static_DMN_201(count,1) = hh(2,1);
    ee = corrcoef(squeeze(ROI_ts_DMN(211:300,r1)),squeeze(LC_ts(211:300,1)));
    all_LC_static_DMN_202(count,1) = ee(2,1);    
    ll = corrcoef(squeeze(ROI_ts_DMN(316:end,r1)),squeeze(LC_ts(316:end,1)));
    all_LC_static_DMN_203(count,1) = ll(2,1);    
    count = count + 1;
end

%FPCN statically correlated with LC during each resting state epoch
count = 1;
for r1 = 1:ROI_FPCN
    yy = corrcoef(squeeze(ROI_ts_FPCN(1:90,r1)),squeeze(LC_ts(1:90,1)));
    all_LC_static_FPCN_none(count,1) = yy(2,1);    
    hh = corrcoef(squeeze(ROI_ts_FPCN(106:195,r1)),squeeze(LC_ts(106:195,1)));
    all_LC_static_FPCN_201(count,1) = hh(2,1);
    ee = corrcoef(squeeze(ROI_ts_FPCN(211:300,r1)),squeeze(LC_ts(211:300,1)));
    all_LC_static_FPCN_202(count,1) = ee(2,1);    
    ll = corrcoef(squeeze(ROI_ts_FPCN(316:end,r1)),squeeze(LC_ts(316:end,1)));
    all_LC_static_FPCN_203(count,1) = ll(2,1);    
    count = count + 1;
end

%DAN statically correlated with LC during each resting state epoch
count = 1;
for r1 = 1:ROI_DAN
    yy = corrcoef(squeeze(ROI_ts_DAN(1:90,r1)),squeeze(LC_ts(1:90,1)));
    all_LC_static_DAN_none(count,1) = yy(2,1);    
    hh = corrcoef(squeeze(ROI_ts_DAN(106:195,r1)),squeeze(LC_ts(106:195,1)));
    all_LC_static_DAN_201(count,1) = hh(2,1);
    ee = corrcoef(squeeze(ROI_ts_DAN(211:300,r1)),squeeze(LC_ts(211:300,1)));
    all_LC_static_DAN_202(count,1) = ee(2,1);    
    ll = corrcoef(squeeze(ROI_ts_DAN(316:end,r1)),squeeze(LC_ts(316:end,1)));
    all_LC_static_DAN_203(count,1) = ll(2,1);    
    count = count + 1;
end

%SN statically correlated with LC during each resting state epoch
count = 1;
for r1 = 1:ROI_SN
    yy = corrcoef(squeeze(ROI_ts_SN(1:90,r1)),squeeze(LC_ts(1:90,1)));
    all_LC_static_SN_none(count,1) = yy(2,1);    
    hh = corrcoef(squeeze(ROI_ts_SN(106:195,r1)),squeeze(LC_ts(106:195,1)));
    all_LC_static_SN_201(count,1) = hh(2,1);
    ee = corrcoef(squeeze(ROI_ts_SN(211:300,r1)),squeeze(LC_ts(211:300,1)));
    all_LC_static_SN_202(count,1) = ee(2,1);    
    ll = corrcoef(squeeze(ROI_ts_SN(316:end,r1)),squeeze(LC_ts(316:end,1)));
    all_LC_static_SN_203(count,1) = ll(2,1);    
    count = count + 1;
end


%% 

%Sliding Dynamic Correlation Analysis Using Overlapping Time Windows

%In this section of code, I am performing a dynamic correlation analysis
%among the network ROIs' BOLD signal (to obtain the within-network
%connectivity) as well as between these networks' BOLD signals and that of
%the LC (to obtain the LC-network connectivity). I use a sliding time
%window approach with window size delta_t and shave off the first and last
%delta_t_2 time points of each time series.

%In this section of code, the sliding windows overlap

delta_t = 40;
delta_t_2 = 10;

%Within-network DMN first resting state epoch (no squeeze condition)
for t = 1:90-delta_t
    for r1 = 1:ROI_DMN
        for r2 = 1:ROI_DMN
            dd = corrcoef(squeeze(ROI_ts_DMN(t:t+delta_t,r1)),squeeze(ROI_ts_DMN(t:t+delta_t,r2)));
            ROI_dynamic_DMN_none(r1,r2,t) = dd(2,1);
        end
    end
end

kk = size(ROI_dynamic_DMN_none);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_DMN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_DMN_none(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_DMN_none(k,1) = mean(output(:,1));
    clear output
    clear count
end
save('avg_ROI_dynamic_DMN_none.mat','avg_ROI_dynamic_DMN_none')
average_avg_ROI_dynamic_DMN_none = mean(avg_ROI_dynamic_DMN_none);

%Within-network DMN second resting state epoch (after first squeeze)
for t = 106:195-delta_t
    for r1 = 1:ROI_DMN
        for r2 = 1:ROI_DMN
            dd = corrcoef(squeeze(ROI_ts_DMN(t:t+delta_t,r1)),squeeze(ROI_ts_DMN(t:t+delta_t,r2)));
            ROI_dynamic_DMN_201(r1,r2,t) = dd(2,1);
        end
    end
end

ROI_dynamic_DMN_201_2 = ROI_dynamic_DMN_201(:,:,106:end);
kk = size(ROI_dynamic_DMN_201_2);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_DMN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_DMN_201_2(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_DMN_201(k,1) = mean(output(:,1));
    clear output
    clear count
end
save('avg_ROI_dynamic_DMN_201.mat','avg_ROI_dynamic_DMN_201')
average_avg_ROI_dynamic_DMN_20 = mean(avg_ROI_dynamic_DMN_201);

%Within-network DMN third resting state epoch (after second squeeze)
for t = 211:300-delta_t
    for r1 = 1:ROI_DMN
        for r2 = 1:ROI_DMN
            dd = corrcoef(squeeze(ROI_ts_DMN(t:t+delta_t,r1)),squeeze(ROI_ts_DMN(t:t+delta_t,r2)));
            ROI_dynamic_DMN_202(r1,r2,t) = dd(2,1);
        end
    end
end

ROI_dynamic_DMN_202_2 = ROI_dynamic_DMN_202(:,:,211:end);
kk = size(ROI_dynamic_DMN_202_2);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_DMN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_DMN_202_2(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_DMN_202(k,1) = mean(output(:,1));
    clear output
    clear count
end
save('avg_ROI_dynamic_DMN_202.mat','avg_ROI_dynamic_DMN_202')
average_avg_ROI_dynamic_DMN_202 = mean(avg_ROI_dynamic_DMN_202);

%Within-network DMN fourth resting state epoch (after third squeeze)
for t = 316:nt-delta_t
    for r1 = 1:ROI_DMN
        for r2 = 1:ROI_DMN
            dd = corrcoef(squeeze(ROI_ts_DMN(t:t+delta_t,r1)),squeeze(ROI_ts_DMN(t:t+delta_t,r2)));
            ROI_dynamic_DMN_203(r1,r2,t) = dd(2,1);
        end
    end
end

ROI_dynamic_DMN_203_2 = ROI_dynamic_DMN_203(:,:,316:end);
kk = size(ROI_dynamic_DMN_203_2);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_DMN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_DMN_203_2(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_DMN_203(k,1) = mean(output(:,1));
    clear output
    clear count
end
save('avg_ROI_dynamic_DMN_203.mat','avg_ROI_dynamic_DMN_203')
average_avg_ROI_dynamic_DMN_203 = mean(avg_ROI_dynamic_DMN_203);

%Within-network FPCN first resting state epoch (no squeeze condition)
for t = 1:90-delta_t
    for r1 = 1:ROI_FPCN
        for r2 = 1:ROI_FPCN
            dd = corrcoef(squeeze(ROI_ts_FPCN(t:t+delta_t,r1)),squeeze(ROI_ts_FPCN(t:t+delta_t,r2)));
            ROI_dynamic_FPCN_none(r1,r2,t) = dd(2,1);
        end
    end
end

kk = size(ROI_dynamic_FPCN_none);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_FPCN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_FPCN_none(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_FPCN_none(k,1) = mean(output(:,1));
    clear output
    clear count
end
save('avg_ROI_dynamic_FPCN_none.mat','avg_ROI_dynamic_FPCN_none')
average_avg_ROI_dynamic_FPCN_none = mean(avg_ROI_dynamic_FPCN_none);

%Within-network FPCN second resting state epoch (after first squeeze)
for t = 106:195-delta_t
    for r1 = 1:ROI_FPCN
        for r2 = 1:ROI_FPCN
            dd = corrcoef(squeeze(ROI_ts_FPCN(t:t+delta_t,r1)),squeeze(ROI_ts_FPCN(t:t+delta_t,r2)));
            ROI_dynamic_FPCN_201(r1,r2,t) = dd(2,1);
        end
    end
end

ROI_dynamic_FPCN_201_2 = ROI_dynamic_FPCN_201(:,:,106:end);
kk = size(ROI_dynamic_FPCN_201_2);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_FPCN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_FPCN_201_2(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_FPCN_201(k,1) = mean(output(:,1));
    clear output
    clear count
end
save('avg_ROI_dynamic_FPCN_201.mat','avg_ROI_dynamic_FPCN_201')
average_avg_ROI_dynamic_FPCN_201 = mean(avg_ROI_dynamic_FPCN_201);

%Within-network FPCN third resting state epoch (after second squeeze)
for t = 211:300-delta_t
    for r1 = 1:ROI_FPCN
        for r2 = 1:ROI_FPCN
            dd = corrcoef(squeeze(ROI_ts_FPCN(t:t+delta_t,r1)),squeeze(ROI_ts_FPCN(t:t+delta_t,r2)));
            ROI_dynamic_FPCN_202(r1,r2,t) = dd(2,1);
        end
    end
end

ROI_dynamic_FPCN_202_2 = ROI_dynamic_FPCN_202(:,:,211:end);
kk = size(ROI_dynamic_FPCN_202_2);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_FPCN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_FPCN_202_2(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_FPCN_202(k,1) = mean(output(:,1));
    clear output
    clear count
end
save('avg_ROI_dynamic_FPCN_202.mat','avg_ROI_dynamic_FPCN_202')
average_avg_ROI_dynamic_FPCN_202 = mean(avg_ROI_dynamic_FPCN_202);

%Within-network DMN fourth resting state epoch (after third squeeze)
for t = 316:nt-delta_t
    for r1 = 1:ROI_FPCN
        for r2 = 1:ROI_FPCN
            dd = corrcoef(squeeze(ROI_ts_FPCN(t:t+delta_t,r1)),squeeze(ROI_ts_FPCN(t:t+delta_t,r2)));
            ROI_dynamic_FPCN_203(r1,r2,t) = dd(2,1);
        end
    end
end

ROI_dynamic_FPCN_203_2 = ROI_dynamic_FPCN_203(:,:,316:end);
kk = size(ROI_dynamic_FPCN_203_2);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_FPCN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_FPCN_203_2(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_FPCN_203(k,1) = mean(output(:,1));
    clear output
    clear count
end
save('avg_ROI_dynamic_FPCN_203.mat','avg_ROI_dynamic_FPCN_203')
average_avg_ROI_dynamic_FPCN_203 = mean(avg_ROI_dynamic_FPCN_203);


%Within-network DAN first resting state epoch (no squeeze condition)
for t = 1:90-delta_t
    for r1 = 1:ROI_DAN
        for r2 = 1:ROI_DAN
            dd = corrcoef(squeeze(ROI_ts_DAN(t:t+delta_t,r1)),squeeze(ROI_ts_DAN(t:t+delta_t,r2)));
            ROI_dynamic_DAN_none(r1,r2,t) = dd(2,1);
        end
    end
end

kk = size(ROI_dynamic_DAN_none);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_DAN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_DAN_none(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_DAN_none(k,1) = mean(output(:,1));
    clear output
    clear count
end
save('avg_ROI_dynamic_DAN_none.mat','avg_ROI_dynamic_DAN_none')
average_avg_ROI_dynamic_DAN_none = mean(avg_ROI_dynamic_DAN_none);

%Within-network DAN second resting state epoch (after first squeeze)
for t = 106:195-delta_t
    for r1 = 1:ROI_DAN
        for r2 = 1:ROI_DAN
            dd = corrcoef(squeeze(ROI_ts_DAN(t:t+delta_t,r1)),squeeze(ROI_ts_DAN(t:t+delta_t,r2)));
            ROI_dynamic_DAN_201(r1,r2,t) = dd(2,1);
        end
    end
end

ROI_dynamic_DAN_201_2 = ROI_dynamic_DAN_201(:,:,106:end);
kk = size(ROI_dynamic_DAN_201_2);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_DAN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_DAN_201_2(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_DAN_201(k,1) = mean(output(:,1));
    clear output
    clear count
end
save('avg_ROI_dynamic_DAN_201.mat','avg_ROI_dynamic_DAN_201')
average_avg_ROI_dynamic_DAN_201 = mean(avg_ROI_dynamic_DAN_201);

%Within-network DAN third resting state epoch (after second squeeze)
for t = 211:300-delta_t
    for r1 = 1:ROI_DAN
        for r2 = 1:ROI_DAN
            dd = corrcoef(squeeze(ROI_ts_DAN(t:t+delta_t,r1)),squeeze(ROI_ts_DAN(t:t+delta_t,r2)));
            ROI_dynamic_DAN_202(r1,r2,t) = dd(2,1);
        end
    end
end

ROI_dynamic_DAN_202_2 = ROI_dynamic_DAN_202(:,:,211:end);
kk = size(ROI_dynamic_DAN_202_2);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_DAN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_DAN_202_2(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_DAN_202(k,1) = mean(output(:,1));
    clear output
    clear count
end
save('avg_ROI_dynamic_DAN_202.mat','avg_ROI_dynamic_DAN_202')
average_avg_ROI_dynamic_DAN_202 = mean(avg_ROI_dynamic_DAN_202);

%Within-network DAN fourth resting state epoch (after third squeeze)
for t = 316:nt-delta_t
    for r1 = 1:ROI_DAN
        for r2 = 1:ROI_DAN
            dd = corrcoef(squeeze(ROI_ts_DAN(t:t+delta_t,r1)),squeeze(ROI_ts_DAN(t:t+delta_t,r2)));
            ROI_dynamic_DAN_203(r1,r2,t) = dd(2,1);
        end
    end
end

ROI_dynamic_DAN_203_2 = ROI_dynamic_DAN_203(:,:,316:end);
kk = size(ROI_dynamic_DAN_203_2);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_DAN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_DAN_203_2(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_DAN_203(k,1) = mean(output(:,1));
    clear output
    clear count
end
save('avg_ROI_dynamic_DAN_203.mat','avg_ROI_dynamic_DAN_203')
average_avg_ROI_dynamic_DAN_203 = mean(avg_ROI_dynamic_DAN_203);

%Within-network SN first resting state epoch (no squeeze condition)
for t = 1:90-delta_t
    for r1 = 1:ROI_SN
        for r2 = 1:ROI_SN
            dd = corrcoef(squeeze(ROI_ts_SN(t:t+delta_t,r1)),squeeze(ROI_ts_SN(t:t+delta_t,r2)));
            ROI_dynamic_SN_none(r1,r2,t) = dd(2,1);
        end
    end
end

kk = size(ROI_dynamic_SN_none);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_SN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_SN_none(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_SN_none(k,1) = mean(output(:,1));
    clear output
    clear count
end
save('avg_ROI_dynamic_SN_none.mat','avg_ROI_dynamic_SN_none')
average_avg_ROI_dynamic_DMN_none = mean(avg_ROI_dynamic_SN_none);

%Within-network SN second resting state epoch (after first squeeze)
for t = 106:195-delta_t
    for r1 = 1:ROI_SN
        for r2 = 1:ROI_SN
            dd = corrcoef(squeeze(ROI_ts_SN(t:t+delta_t,r1)),squeeze(ROI_ts_SN(t:t+delta_t,r2)));
            ROI_dynamic_SN_201(r1,r2,t) = dd(2,1);
        end
    end
end

ROI_dynamic_SN_201_2 = ROI_dynamic_SN_201(:,:,106:end);
kk = size(ROI_dynamic_SN_201_2);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_SN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_SN_201_2(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_SN_201(k,1) = mean(output(:,1));
    clear output
    clear count
end
save('avg_ROI_dynamic_SN_201.mat','avg_ROI_dynamic_SN_201')
average_avg_ROI_dynamic_SN_201 = mean(avg_ROI_dynamic_SN_201);

%Within-network SN third resting state epoch (after second squeeze)
for t = 211:300-delta_t
    for r1 = 1:ROI_SN
        for r2 = 1:ROI_SN
            dd = corrcoef(squeeze(ROI_ts_SN(t:t+delta_t,r1)),squeeze(ROI_ts_SN(t:t+delta_t,r2)));
            ROI_dynamic_SN_202(r1,r2,t) = dd(2,1);
        end
    end
end

ROI_dynamic_SN_202_2 = ROI_dynamic_SN_202(:,:,211:end);
kk = size(ROI_dynamic_SN_202_2);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_SN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_SN_202_2(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_SN_202(k,1) = mean(output(:,1));
    clear output
    clear count
end
save('avg_ROI_dynamic_SN_202.mat','avg_ROI_dynamic_SN_202')
average_avg_ROI_dynamic_SN_202 = mean(avg_ROI_dynamic_SN_202);

%Within-network SN fourth resting state epoch (after third squeeze)
for t = 316:nt-delta_t
    for r1 = 1:ROI_SN
        for r2 = 1:ROI_SN
            dd = corrcoef(squeeze(ROI_ts_SN(t:t+delta_t,r1)),squeeze(ROI_ts_SN(t:t+delta_t,r2)));
            ROI_dynamic_SN_203(r1,r2,t) = dd(2,1);
        end
    end
end

ROI_dynamic_SN_203_2 = ROI_dynamic_SN_203(:,:,316:end);
kk = size(ROI_dynamic_SN_203_2);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_SN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_SN_203_2(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_SN_203(k,1) = mean(output(:,1));
    clear output
    clear count
end
save('avg_ROI_dynamic_SN_203.mat','avg_ROI_dynamic_SN_203')
average_avg_ROI_dynamic_SN_203 = mean(avg_ROI_dynamic_SN_203);


%%

%Dynamic Sliding Window Correlation Analysis With Non-Overlapping Windows

%In this section of code, I use non-overlapping windows and perform a
%correlation within these windows

num_col = 4;
delta_t = 40;
delta_t_2 = 10;
length_resting_state = 90;
length_squeeze = 15+1;
%Create equal sized windows with the proper time points for each resting
%state epoch
nonoverlap_vector_none = 1:delta_t_2:90+delta_t_2-1; %First resting state epoch
nonoverlap_vector_201 = 106:delta_t_2:195+delta_t_2-1; %Second resting state epoch
nonoverlap_vector_202 = 211:delta_t_2:300+delta_t_2-1; %Third resting state epoch
nonoverlap_vector_203 = 316:delta_t_2:405+delta_t_2-1; %Fourth resting state epoch
nonoverlap_vector_all = cat(1,nonoverlap_vector_none, nonoverlap_vector_201, nonoverlap_vector_202, nonoverlap_vector_203)';

length_vector = size(nonoverlap_vector_all);

%Within-network DMN first resting state epoch (no squeeze condition)
for t = 1:length_vector(1)-1
    for r1 = 1:ROI_DMN
        for r2 = 1:ROI_DMN
            v = nonoverlap_vector_none(t);
            vv = nonoverlap_vector_none(t+1);
            [ll, p_ll] = corrcoef(squeeze(ROI_ts_DMN(v:vv,r1)),squeeze(ROI_ts_DMN(v:vv,r2)));
            ROI_dynamic_DMN_none_nonoverlap(r1,r2,t) = ll(2,1);
            ROI_dynamic_DMN_none_nonoverlap_p(r1,r2,t) = p_ll(2,1);
        end
    end
end

kk = size(ROI_dynamic_DMN_none_nonoverlap);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_DMN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_DMN_none_nonoverlap(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_DMN_none_nonoverlap(k,1) = mean(output(:,1));
    clear output
    clear count
end

%Within-network DMN second resting state epoch (after first squeeze)
for t = 1:length_vector(1)-1
    for r1 = 1:ROI_DMN
        for r2 = 1:ROI_DMN
            v = nonoverlap_vector_201(t);
            vv = nonoverlap_vector_201(t+1);
            [ll, p_ll] = corrcoef(squeeze(ROI_ts_DMN(v:vv,r1)),squeeze(ROI_ts_DMN(v:vv,r2)));
            ROI_dynamic_DMN_201_nonoverlap(r1,r2,t) = ll(2,1);
            ROI_dynamic_DMN_201_nonoverlap_p(r1,r2,t) = p_ll(2,1);
        end
    end
end

kk = size(ROI_dynamic_DMN_201_nonoverlap);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_DMN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_DMN_201_nonoverlap(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_DMN_201_nonoverlap(k,1) = mean(output(:,1));
    clear output
    clear count
end

%Within-network DMN third resting state epoch (after second squeeze)
for t = 1:length_vector(1)-1
    for r1 = 1:ROI_DMN
        for r2 = 1:ROI_DMN
            v = nonoverlap_vector_202(t);
            vv = nonoverlap_vector_202(t+1);
            [ll, p_ll] = corrcoef(squeeze(ROI_ts_DMN(v:vv,r1)),squeeze(ROI_ts_DMN(v:vv,r2)));
            ROI_dynamic_DMN_202_nonoverlap(r1,r2,t) = ll(2,1);
            ROI_dynamic_DMN_202_nonoverlap_p(r1,r2,t) = p_ll(2,1);
        end
    end
end

kk = size(ROI_dynamic_DMN_202_nonoverlap);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_DMN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_DMN_202_nonoverlap(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_DMN_202_nonoverlap(k,1) = mean(output(:,1));
    clear output
    clear count
end

%Within-network DMN fourth resting state epoch (after third squeeze)
for t = 1:length_vector(1)-1
    for r1 = 1:ROI_DMN
        for r2 = 1:ROI_DMN
            v = nonoverlap_vector_203(t);
            vv = nonoverlap_vector_203(t+1);
            [ll, p_ll] = corrcoef(squeeze(ROI_ts_DMN(v:vv,r1)),squeeze(ROI_ts_DMN(v:vv,r2)));
            ROI_dynamic_DMN_203_nonoverlap(r1,r2,t) = ll(2,1);
            ROI_dynamic_DMN_203_nonoverlap_p(r1,r2,t) = p_ll(2,1);
        end
    end
end

kk = size(ROI_dynamic_DMN_203_nonoverlap);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_DMN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_DMN_203_nonoverlap(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_DMN_203_nonoverlap(k,1) = mean(output(:,1));
    clear output
    clear count
end

avg_ROI_dynamic_DMN_all_nonoverlap = [avg_ROI_dynamic_DMN_none_nonoverlap,avg_ROI_dynamic_DMN_201_nonoverlap,avg_ROI_dynamic_DMN_202_nonoverlap,avg_ROI_dynamic_DMN_203_nonoverlap];
save('avg_ROI_dynamic_DMN_all_nonoverlap.mat','avg_ROI_dynamic_DMN_all_nonoverlap');

%Within-network DAN first resting state epoch (no squeeze condition)
for t = 1:length_vector(1)-1
    for r1 = 1:ROI_DAN
        for r2 = 1:ROI_DAN
            v = nonoverlap_vector_none(t);
            vv = nonoverlap_vector_none(t+1);
            [ll, p_ll] = corrcoef(squeeze(ROI_ts_DAN(v:vv,r1)),squeeze(ROI_ts_DAN(v:vv,r2)));
            ROI_dynamic_DAN_none_nonoverlap(r1,r2,t) = ll(2,1);
            ROI_dynamic_DAN_none_nonoverlap_p(r1,r2,t) = p_ll(2,1);
        end
    end
end

kk = size(ROI_dynamic_DAN_none_nonoverlap);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_DAN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_DAN_none_nonoverlap(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_DAN_none_nonoverlap(k,1) = mean(output(:,1));
    clear output
    clear count
end

%Within-network DAN second resting state epoch (after first squeeze)
for t = 1:length_vector(1)-1
    for r1 = 1:ROI_DAN
        for r2 = 1:ROI_DAN
            v = nonoverlap_vector_201(t);
            vv = nonoverlap_vector_201(t+1);
            [ll, p_ll] = corrcoef(squeeze(ROI_ts_DAN(v:vv,r1)),squeeze(ROI_ts_DAN(v:vv,r2)));
            ROI_dynamic_DAN_201_nonoverlap(r1,r2,t) = ll(2,1);
            ROI_dynamic_DAN_201_nonoverlap_p(r1,r2,t) = p_ll(2,1);
        end
    end
end

kk = size(ROI_dynamic_DAN_201_nonoverlap);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_DAN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_DAN_201_nonoverlap(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_DAN_201_nonoverlap(k,1) = mean(output(:,1));
    clear output
    clear count
end

%Within-network DAN third resting state epoch (after second squeeze)
for t = 1:length_vector(1)-1
    for r1 = 1:ROI_DAN
        for r2 = 1:ROI_DAN
            v = nonoverlap_vector_202(t);
            vv = nonoverlap_vector_202(t+1);
            [ll, p_ll] = corrcoef(squeeze(ROI_ts_DAN(v:vv,r1)),squeeze(ROI_ts_DAN(v:vv,r2)));
            ROI_dynamic_DAN_202_nonoverlap(r1,r2,t) = ll(2,1);
            ROI_dynamic_DAN_202_nonoverlap_p(r1,r2,t) = p_ll(2,1);
        end
    end
end

kk = size(ROI_dynamic_DAN_202_nonoverlap);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_DAN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_DAN_202_nonoverlap(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_DAN_202_nonoverlap(k,1) = mean(output(:,1));
    clear output
    clear count
end

%Within-network DAN fourth resting state epoch (after third squeeze)
for t = 1:length_vector(1)-1
    for r1 = 1:ROI_DAN
        for r2 = 1:ROI_DAN
            v = nonoverlap_vector_203(t);
            vv = nonoverlap_vector_203(t+1);
            [ll, p_ll] = corrcoef(squeeze(ROI_ts_DAN(v:vv,r1)),squeeze(ROI_ts_DAN(v:vv,r2)));
            ROI_dynamic_DAN_203_nonoverlap(r1,r2,t) = ll(2,1);
            ROI_dynamic_DAN_203_nonoverlap_p(r1,r2,t) = p_ll(2,1);
        end
    end
end

kk = size(ROI_dynamic_DAN_203_nonoverlap);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_DAN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_DAN_203_nonoverlap(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_DAN_203_nonoverlap(k,1) = mean(output(:,1));
    clear output
    clear count
end

avg_ROI_dynamic_DAN_all_nonoverlap = [avg_ROI_dynamic_DAN_none_nonoverlap,avg_ROI_dynamic_DAN_201_nonoverlap,avg_ROI_dynamic_DAN_202_nonoverlap,avg_ROI_dynamic_DAN_203_nonoverlap];
save('avg_ROI_dynamic_DAN_all_nonoverlap.mat','avg_ROI_dynamic_DAN_all_nonoverlap');

%Within-network SN first resting state epoch (no squeeze condition)
for t = 1:length_vector(1)-1
    for r1 = 1:ROI_SN
        for r2 = 1:ROI_SN
            v = nonoverlap_vector_none(t);
            vv = nonoverlap_vector_none(t+1);
            [ll, p_ll] = corrcoef(squeeze(ROI_ts_SN(v:vv,r1)),squeeze(ROI_ts_SN(v:vv,r2)));
            ROI_dynamic_SN_none_nonoverlap(r1,r2,t) = ll(2,1);
            ROI_dynamic_SN_none_nonoverlap_p(r1,r2,t) = p_ll(2,1);
        end
    end
end

kk = size(ROI_dynamic_SN_none_nonoverlap);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_SN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_SN_none_nonoverlap(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_SN_none_nonoverlap(k,1) = mean(output(:,1));
    clear output
    clear count
end

%Within-network SN second resting state epoch (after first squeeze)
for t = 1:length_vector(1)-1
    for r1 = 1:ROI_SN
        for r2 = 1:ROI_SN
            v = nonoverlap_vector_201(t);
            vv = nonoverlap_vector_201(t+1);
            [ll, p_ll] = corrcoef(squeeze(ROI_ts_SN(v:vv,r1)),squeeze(ROI_ts_SN(v:vv,r2)));
            ROI_dynamic_SN_201_nonoverlap(r1,r2,t) = ll(2,1);
            ROI_dynamic_SN_201_nonoverlap_p(r1,r2,t) = p_ll(2,1);
        end
    end
end

kk = size(ROI_dynamic_SN_201_nonoverlap);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_SN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_SN_201_nonoverlap(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_SN_201_nonoverlap(k,1) = mean(output(:,1));
    clear output
    clear count
end

%Within-network SN third resting state epoch (after second squeeze)
for t = 1:length_vector(1)-1
    for r1 = 1:ROI_SN
        for r2 = 1:ROI_SN
            v = nonoverlap_vector_202(t);
            vv = nonoverlap_vector_202(t+1);
            [ll, p_ll] = corrcoef(squeeze(ROI_ts_SN(v:vv,r1)),squeeze(ROI_ts_SN(v:vv,r2)));
            ROI_dynamic_SN_202_nonoverlap(r1,r2,t) = ll(2,1);
            ROI_dynamic_SN_202_nonoverlap_p(r1,r2,t) = p_ll(2,1);
        end
    end
end

kk = size(ROI_dynamic_SN_202_nonoverlap);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_SN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_SN_202_nonoverlap(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_SN_202_nonoverlap(k,1) = mean(output(:,1));
    clear output
    clear count
end

%Within-network SN fourth resting state epoch (after third squeeze)
for t = 1:length_vector(1)-1
    for r1 = 1:ROI_SN
        for r2 = 1:ROI_SN
            v = nonoverlap_vector_203(t);
            vv = nonoverlap_vector_203(t+1);
            [ll, p_ll] = corrcoef(squeeze(ROI_ts_SN(v:vv,r1)),squeeze(ROI_ts_SN(v:vv,r2)));
            ROI_dynamic_SN_203_nonoverlap(r1,r2,t) = ll(2,1);
            ROI_dynamic_SN_203_nonoverlap_p(r1,r2,t) = p_ll(2,1);
        end
    end
end

kk = size(ROI_dynamic_SN_203_nonoverlap);
for k = 1:kk(3)
    count = 1;
    for r1 = 2:ROI_SN
        for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_SN_203_nonoverlap(r1,r2,k));
            count = count+1;
        end
    end
    avg_ROI_dynamic_SN_203_nonoverlap(k,1) = mean(output(:,1));
    clear output
    clear count
end

avg_ROI_dynamic_SN_all_nonoverlap = [avg_ROI_dynamic_SN_none_nonoverlap,avg_ROI_dynamic_SN_201_nonoverlap,avg_ROI_dynamic_SN_202_nonoverlap,avg_ROI_dynamic_SN_203_nonoverlap];
save('avg_ROI_dynamic_SN_all_nonoverlap.mat','avg_ROI_dynamic_SN_all_nonoverlap');



%% 

%Dynamically correlate the LC with all the networks using the sliding time
%window approach and overlapping windows

ROI_DMN = double(ROI_DMN);
ROI_FPCN = double(ROI_FPCN);
ROI_DAN = double(ROI_DAN);
ROI_SN = double(ROI_SN);
% delta_t_new = 20;

h_size = size(avg_ROI_dynamic_DMN_none);

%DMN
for t = 1:h_size-delta_t       
        cc = corrcoef(squeeze(avg_ROI_dynamic_DMN_none(t:t+delta_t,1)),squeeze(LC_ts_none(t:t+delta_t,1)));
        LC_dynamic_DMN_none(t,1) = cc(2,1); 
        bb = corrcoef(squeeze(avg_ROI_dynamic_DMN_201(t:t+delta_t,1)),squeeze(LC_ts_201(t:t+delta_t,1)));
        LC_dynamic_DMN_201(t,1) = bb(2,1);
        aa = corrcoef(squeeze(avg_ROI_dynamic_DMN_202(t:t+delta_t,1)),squeeze(LC_ts_202(t:t+delta_t,1)));
        LC_dynamic_DMN_202(t,1) = aa(2,1);
        gg = corrcoef(squeeze(avg_ROI_dynamic_DMN_203(t:t+delta_t,1)),squeeze(LC_ts_203(t:t+delta_t,1)));
        LC_dynamic_DMN_203(t,1) = gg(2,1);
end
save('LC_dynamic_DMN_none.mat','LC_dynamic_DMN_none')
avg_LC_dynamic_DMN_none = mean(LC_dynamic_DMN_none);
save('LC_dynamic_DMN_201.mat','LC_dynamic_DMN_201')
avg_LC_dynamic_DMN_201 = mean(LC_dynamic_DMN_201);
save('LC_dynamic_DMN_202.mat','LC_dynamic_DMN_202')
avg_LC_dynamic_DMN_202 = mean(LC_dynamic_DMN_202);
save('LC_dynamic_DMN_203.mat','LC_dynamic_DMN_203')
avg_LC_dynamic_DMN_203 = mean(LC_dynamic_DMN_203);


%FPCN
for t = 1:h_size-delta_t       
        cc = corrcoef(squeeze(avg_ROI_dynamic_FPCN_none(t:t+delta_t,1)),squeeze(LC_ts_none(t:t+delta_t,1)));
        LC_dynamic_FPCN_none(t,1) = cc(2,1); 
        bb = corrcoef(squeeze(avg_ROI_dynamic_FPCN_201(t:t+delta_t,1)),squeeze(LC_ts_201(t:t+delta_t,1)));
        LC_dynamic_FPCN_201(t,1) = bb(2,1);
        aa = corrcoef(squeeze(avg_ROI_dynamic_FPCN_202(t:t+delta_t,1)),squeeze(LC_ts_202(t:t+delta_t,1)));
        LC_dynamic_FPCN_202(t,1) = aa(2,1);
        gg = corrcoef(squeeze(avg_ROI_dynamic_FPCN_203(t:t+delta_t,1)),squeeze(LC_ts_203(t:t+delta_t,1)));
        LC_dynamic_FPCN_203(t,1) = gg(2,1);
end
save('LC_dynamic_FPCN_none.mat','LC_dynamic_FPCN_none')
avg_LC_dynamic_FPCN_none = mean(LC_dynamic_FPCN_none);
save('LC_dynamic_FPCN_201.mat','LC_dynamic_FPCN_201')
avg_LC_dynamic_FPCN_201 = mean(LC_dynamic_FPCN_201);
save('LC_dynamic_FPCN_202.mat','LC_dynamic_FPCN_202')
avg_LC_dynamic_FPCN_202 = mean(LC_dynamic_FPCN_202);
save('LC_dynamic_FPCN_203.mat','LC_dynamic_FPCN_203')
avg_LC_dynamic_FPCN_203 = mean(LC_dynamic_FPCN_203);


%DAN
for t = 1:h_size-delta_t       
        cc = corrcoef(squeeze(avg_ROI_dynamic_DAN_none(t:t+delta_t,1)),squeeze(LC_ts_none(t:t+delta_t,1)));
        LC_dynamic_DAN_none(t,1) = cc(2,1); 
        bb = corrcoef(squeeze(avg_ROI_dynamic_DAN_201(t:t+delta_t,1)),squeeze(LC_ts_201(t:t+delta_t,1)));
        LC_dynamic_DAN_201(t,1) = bb(2,1);
        aa = corrcoef(squeeze(avg_ROI_dynamic_DAN_202(t:t+delta_t,1)),squeeze(LC_ts_202(t:t+delta_t,1)));
        LC_dynamic_DAN_202(t,1) = aa(2,1);
        gg = corrcoef(squeeze(avg_ROI_dynamic_DAN_203(t:t+delta_t,1)),squeeze(LC_ts_203(t:t+delta_t,1)));
        LC_dynamic_DAN_203(t,1) = gg(2,1);
end
save('LC_dynamic_DAN_none.mat','LC_dynamic_DAN_none')
avg_LC_dynamic_DAN_none = mean(LC_dynamic_DAN_none);
save('LC_dynamic_DAN_201.mat','LC_dynamic_DAN_201')
avg_LC_dynamic_DAN_201 = mean(LC_dynamic_DAN_201);
save('LC_dynamic_DAN_202.mat','LC_dynamic_DAN_202')
avg_LC_dynamic_DAN_202 = mean(LC_dynamic_DAN_202);
save('LC_dynamic_DAN_203.mat','LC_dynamic_DAN_203')
avg_LC_dynamic_DAN_203 = mean(LC_dynamic_DAN_203);


%SN
for t = 1:h_size-delta_t       
        cc = corrcoef(squeeze(avg_ROI_dynamic_SN_none(t:t+delta_t,1)),squeeze(LC_ts_none(t:t+delta_t,1)));
        LC_dynamic_SN_none(t,1) = cc(2,1); 
        bb = corrcoef(squeeze(avg_ROI_dynamic_SN_201(t:t+delta_t,1)),squeeze(LC_ts_201(t:t+delta_t,1)));
        LC_dynamic_SN_201(t,1) = bb(2,1);
        aa = corrcoef(squeeze(avg_ROI_dynamic_SN_202(t:t+delta_t,1)),squeeze(LC_ts_202(t:t+delta_t,1)));
        LC_dynamic_SN_202(t,1) = aa(2,1);
        gg = corrcoef(squeeze(avg_ROI_dynamic_SN_203(t:t+delta_t,1)),squeeze(LC_ts_203(t:t+delta_t,1)));
        LC_dynamic_SN_203(t,1) = gg(2,1);
end
save('LC_dynamic_SN_none.mat','LC_dynamic_SN_none')
avg_LC_dynamic_SN_none = mean(LC_dynamic_SN_none);
save('LC_dynamic_SN_201.mat','LC_dynamic_SN_201')
avg_LC_dynamic_SN_201 = mean(LC_dynamic_SN_201);
save('LC_dynamic_SN_202.mat','LC_dynamic_SN_202')
avg_LC_dynamic_SN_202 = mean(LC_dynamic_SN_202);
save('LC_dynamic_SN_203.mat','LC_dynamic_SN_203')
avg_LC_dynamic_SN_203 = mean(LC_dynamic_SN_203);
%%

%Eye Tracking and Pupillometry Data

%In this section of code, I am loading in the pupillometry data and
%separating it into the resting state epochs and squeezing sessions.
test = load('eye_tracking_data_no_blinks.mat');
eye_tracking_no_blinks = test.eye_tracking_data_no_blinks;
eye_tracking_none = eye_tracking_no_blinks(1:90,:);
eye_tracking_during_201 = eye_tracking_no_blinks(91:105,:);
eye_tracking_201 = eye_tracking_no_blinks(106:195,:);
eye_tracking_during_202 = eye_tracking_no_blinks(196:210,:);
eye_tracking_202 = eye_tracking_no_blinks(211:300,:);
eye_tracking_during_203 = eye_tracking_no_blinks(301:315,:);
eye_tracking_203 = eye_tracking_no_blinks(316:end,:);
pupil_size_leftpp_hori_all = eye_tracking_no_blinks(:,1);

%Extract the data for left pupil horizontal data, left pupil vertical data,
%right pupil horizontal data, right pupil vertical data
pupil_size_leftpp_hori = [eye_tracking_none(:,1), eye_tracking_201(:,1), eye_tracking_202(:,1), eye_tracking_203(1:90,1)];
pupil_size_leftpp_vert = [eye_tracking_none(:,2), eye_tracking_201(:,2), eye_tracking_202(:,2), eye_tracking_203(1:90,2)];
pupil_size_rightpp_hori = [eye_tracking_none(:,3), eye_tracking_201(:,3), eye_tracking_202(:,3), eye_tracking_203(1:90,3)];
pupil_size_rightpp_vert = [eye_tracking_none(:,4), eye_tracking_201(:,4), eye_tracking_202(:,4), eye_tracking_203(1:90,4)];

save('pupil_size_leftpp_hori.mat','pupil_size_leftpp_hori')
save('pupil_size_leftpp_hori_all.mat','pupil_size_leftpp_hori_all');

delta_t = 40;
delta_t_2 = 10;
length_resting_state = 90;
nonoverlap_vector_none = 1:delta_t_2:90+delta_t_2-1;

%Pupillometry data for all four measurements should be identical, so for
%now I will focus only on the left horizontal data for analysis. Now find
%the mean pupil size within each of the overlapping windows
for t = 1:90-delta_t
    for c = 1:num_col
        pupil_size_leftpp_hori_dynamics(t,c) = mean(pupil_size_leftpp_hori(t:t+delta_t,c));
    end
end
save('pupil_size_leftpp_hori_dynamics.mat','pupil_size_leftpp_hori_dynamics');

%Find the mean pupil size within each of the nonoverlapping windows
for c = 1:num_col
    for t = 1:delta_t_2-1
        v = nonoverlap_vector_none(t);
        vv = nonoverlap_vector_none(t+1);
        pupil_size_leftpp_hori_dynamics_nonoverlap(t,c) = mean(pupil_size_leftpp_hori(v:vv-1,c));
    end
end
save('pupil_size_leftpp_hori_dynamics_nonoverlap.mat','pupil_size_leftpp_hori_dynamics_nonoverlap')


%%

%To make plotting the data easier, I will concatenate the data into a
%single array/matrix to make looping over it easier. 
rs_epoch = [0 1 2 3];

avg_ROI_static_SN = [avg_ROI_static_SN_none, avg_ROI_static_SN_201, avg_ROI_static_SN_202, avg_ROI_static_SN_203];
avg_ROI_static_DMN = [avg_ROI_static_DMN_none, avg_ROI_static_DMN_201, avg_ROI_static_DMN_202, avg_ROI_static_DMN_203];
avg_ROI_static_DAN = [avg_ROI_static_DAN_none, avg_ROI_static_DAN_201, avg_ROI_static_DAN_202, avg_ROI_static_DAN_203];

save('avg_ROI_static_SN.mat','avg_ROI_static_SN')
save('avg_ROI_static_DMN.mat','avg_ROI_static_DMN')
save('avg_ROI_static_DAN.mat','avg_ROI_static_DAN')

avg_ROI_dynamic_DMN = [avg_ROI_dynamic_DMN_none,avg_ROI_dynamic_DMN_201,avg_ROI_dynamic_DMN_202, avg_ROI_dynamic_DMN_203(1:50,:)];
avg_ROI_dynamic_DAN = [avg_ROI_dynamic_DAN_none,avg_ROI_dynamic_DAN_201,avg_ROI_dynamic_DAN_202, avg_ROI_dynamic_DAN_203(1:50,:)];
avg_ROI_dynamic_SN = [avg_ROI_dynamic_SN_none,avg_ROI_dynamic_SN_201,avg_ROI_dynamic_SN_202, avg_ROI_dynamic_SN_203(1:50,:)];

var_ROI_dynamic_DMN_all = var(ROI_dynamic_DMN_none,0,3); var_ROI_dynamic_DMN_all(:,:,2) = var(ROI_dynamic_DMN_201,0,3); var_ROI_dynamic_DMN_all(:,:,3) = var(ROI_dynamic_DMN_202,0,3); var_ROI_dynamic_DMN_all(:,:,4) = var(ROI_dynamic_DMN_203,0,3);
var_ROI_dynamic_DAN_all = var(ROI_dynamic_DAN_none,0,3); var_ROI_dynamic_DAN_all(:,:,2) = var(ROI_dynamic_DAN_201,0,3); var_ROI_dynamic_DAN_all(:,:,3) = var(ROI_dynamic_DAN_202,0,3); var_ROI_dynamic_DAN_all(:,:,4) = var(ROI_dynamic_DAN_203,0,3);
var_ROI_dynamic_SN_all = var(ROI_dynamic_SN_none,0,3); var_ROI_dynamic_SN_all(:,:,2) = var(ROI_dynamic_SN_201,0,3); var_ROI_dynamic_SN_all(:,:,3) = var(ROI_dynamic_SN_202,0,3); var_ROI_dynamic_SN_all(:,:,4) = var(ROI_dynamic_SN_203,0,3);

ROI_static_DMN_all = ROI_static_DMN_none; ROI_static_DMN_all(:,:,2) = ROI_static_DMN_201; ROI_static_DMN_all(:,:,3) = ROI_static_DMN_202; ROI_static_DMN_all(:,:,4) = ROI_static_DMN_203;
ROI_static_DMN_all_p = ROI_static_DMN_none_p; ROI_static_DMN_all_p(:,:,2) = ROI_static_DMN_201_p; ROI_static_DMN_all_p(:,:,3) = ROI_static_DMN_202_p; ROI_static_DMN_all_p(:,:,4) = ROI_static_DMN_203_p;
ROI_static_DAN_all = ROI_static_DAN_none; ROI_static_DAN_all(:,:,2) = ROI_static_DAN_201; ROI_static_DAN_all(:,:,3) = ROI_static_DAN_202; ROI_static_DAN_all(:,:,4) = ROI_static_DAN_203;
ROI_static_DAN_all_p = ROI_static_DAN_none_p; ROI_static_DAN_all_p(:,:,2) = ROI_static_DAN_201_p; ROI_static_DAN_all_p(:,:,3) = ROI_static_DAN_202_p; ROI_static_DAN_all_p(:,:,4) = ROI_static_DAN_203_p;
ROI_static_SN_all = ROI_static_SN_none; ROI_static_SN_all(:,:,2) = ROI_static_SN_201; ROI_static_SN_all(:,:,3) = ROI_static_SN_202; ROI_static_SN_all(:,:,4) = ROI_static_SN_203; 
ROI_static_SN_all_p = ROI_static_SN_none_p; ROI_static_SN_all_p(:,:,2) = ROI_static_SN_201_p; ROI_static_SN_all_p(:,:,3) = ROI_static_SN_202_p; ROI_static_SN_all_p(:,:,4) = ROI_static_SN_203_p;

save('avg_ROI_dynamic_DMN.mat','avg_ROI_dynamic_DMN')
save('avg_ROI_dynamic_DAN.mat','avg_ROI_dynamic_DAN')
save('avg_ROI_dynamic_SN.mat','avg_ROI_dynamic_SN')

LC_ts_all = [LC_ts_none, LC_ts_201, LC_ts_202, LC_ts_203(1:90,1)];
save('LC_ts_all.mat','LC_ts_all')

%I will initialize some vectors to plot against the data.
color_vector = [0.12 0.95, 1; 0.2, 0.75, 1; 0 0 1; 0 0 0];
color_vector_during = [0.12, 0.95,1; 0.1,0.95,0.1; 0.2, 0.75, 1; 0.01,0.8,0.01; 0 0 1; 0.1,0.55,0.1; 0 0 0];  
squeeze_vector = {'No Squeeze','First Squeeze','Second Squeeze','Third Squeeze'};
squeeze_vector_during = {'No Squeeze','During','First Squeeze','During','Second Squeeze','During','Third Squeeze'};
time_vector = [1 90 105 195 210 300 315 405];
vector_during = [0 0.5 1 1.5 2 2.5 3];

LC_dynamic_DMN = [LC_dynamic_DMN_none, LC_dynamic_DMN_201, LC_dynamic_DMN_202, LC_dynamic_DMN_203];
LC_dynamic_DAN = [LC_dynamic_DAN_none, LC_dynamic_DAN_201, LC_dynamic_DAN_202, LC_dynamic_DAN_203];
LC_dynamic_SN = [LC_dynamic_SN_none, LC_dynamic_SN_201, LC_dynamic_SN_202, LC_dynamic_SN_203];
%%

%Now I create plots to visualize the data analyzed above.

%Scatter plot of mean pupil size during both the resting state epochs and the
%squeezing sessions
figure
for i = 1:7
    scatter(vector_during(1,i), mean(pupil_size_leftpp_hori_all(time_vector(i)+1:time_vector(i+1))),200,'filled','MarkerFaceColor',color_vector_during(i,:))
    hold on
    xlabel('Squeeze')
    ylabel('Mean Leftpp Hori')
    title('Mean Pupil Size for Each Squeeze')
    
end
set(gca,'Color',[0.95 0.95 0.95]) 
legend('No Squeeze','During','First Squeeze','During','Second Squeeze','During','Third Squeeze','Location','best')

%Scatter plot of mean LC activity vs. static correlations during each
%resting state epoch (for DMN, DAN, and SN)
figure
for i = 1:4
    subplot(1,3,1)
    set(gca,'Color',[0.95 0.95 0.95])
    scatter(mean(LC_ts_all(:,i)),avg_ROI_static_DMN(:,i),200,'filled','MarkerFaceColor',color_vector(i,:))
    xlabel('Right LC Activity')
    ylabel('Mean Connectivity')
    title('DMN')
    hold on
    subplot(1,3,2)
    set(gca,'Color',[0.95 0.95 0.95])
    scatter(mean(LC_ts_all(:,i)),avg_ROI_static_DAN(:,i),200,'filled','MarkerFaceColor',color_vector(i,:))
    xlabel('Right LC Activity')
    ylabel('Mean Connectivity')
    title('DAN')
    hold on
    subplot(1,3,3)
    set(gca,'Color',[0.95 0.95 0.95])
    scatter(mean(LC_ts_all(:,i)),avg_ROI_static_SN(:,i),200,'filled','MarkerFaceColor',color_vector(i,:))
    xlabel('Right LC Activity')
    ylabel('Mean Connectivity')
    title('SN')
    hold on
end
suptitle('Connectivity vs. LC Activity')

%Scatter plots of within-network dynamic connectivity vs. LC-network
%dynamic connectivity during each resting state epoch (for DMN, DAN, and
%SN)
for i = 1:4
   figure(1)
   subplot(2,2,i)
   ti = sprintf('Squeeze %d.', rs_epoch(i))
   scatter(LC_dynamic_DMN(:,i),avg_ROI_dynamic_DMN(1:10,i),50,'filled')
   hold on
   scatter(mean(LC_dynamic_DMN(:,i)),mean(avg_ROI_dynamic_DMN(1:10,i)),200,'s','filled')
   xlabel('LC and DMN Connectivity')
   ylabel('DMN Connectivity')
   title(ti)
   figure(2)
   subplot(2,2,i)
   ti = sprintf('Squeeze %d.', rs_epoch(i))
   scatter(LC_dynamic_DAN(:,i),avg_ROI_dynamic_DAN(1:10,i),50,'filled')
   hold on
   scatter(mean(LC_dynamic_DAN(:,i)),mean(avg_ROI_dynamic_DAN(1:10,i)),200,'s','filled')
   xlabel('LC and DMN Connectivity')
   ylabel('DAN Connectivity')
   title(ti)
   figure(3)
   subplot(2,2,i)
   ti = sprintf('Squeeze %d.', rs_epoch(i))
   scatter(LC_dynamic_SN(:,i),avg_ROI_dynamic_SN(1:10,i),50,'filled')
   hold on
   scatter(mean(LC_dynamic_SN(:,i)),mean(avg_ROI_dynamic_SN(1:10,i)),200,'s','filled')
   xlabel('LC and DMN Connectivity')
   ylabel('SN Connectivity')
   title(ti)    
end

%Bar graph of the variance of LC BOLD signal during each resting state
%epoch
figure
bar(rs_epoch, var(LC_ts_all))
ylabel('LC Time Series Variance')
title('TS Variance vs. Percent Squeeze')

%Scatter plot of each network static correlation (DMN, DAN, and SN) vs. 
%variance of LC BOLD signal
figure
for i = 1:length(rs_epoch)
    subplot(1,3,1)
    scatter(var(LC_ts_all(:,i)), avg_ROI_static_DMN(1,i),200,'filled','MarkerFaceColor',color_vector(i,:))
    set(gca,'Color',[0.95 0.95 0.95])
    xlabel('LC Time Series Variance')
    ylabel('Mean Connectivity')
    title('Default Mode Network')
    hold on
    subplot(1,3,2)
    scatter(var(LC_ts_all(:,i)), avg_ROI_static_DAN(1,i),200,'filled','MarkerFaceColor',color_vector(i,:))
    set(gca,'Color',[0.95 0.95 0.95])
    xlabel('LC Time Series Variance')
    ylabel('Mean Connectivity')
    title('Dorsal Attention Network')
    hold on
    subplot(1,3,3)
    scatter(var(LC_ts_all(:,i)), avg_ROI_static_SN(1,i),200,'filled','MarkerFaceColor',color_vector(i,:))
    set(gca,'Color',[0.95 0.95 0.95])
    xlabel('LC Time Series Variance')
    ylabel('Mean Connectivity')
    title('Salience Network')
    hold on
end
legend('No Squeeze','First Squeeze','Second Squeeze','Third Squeeze','Location','best')
suptitle('LC Time Series Variance vs. Static Connectivity')

%Bar graphs of the variance of the within-network dynamic connectivity for
%each resting state epoch (for DMN, DAN, and SN)
figure
subplot(1,3,1)
bar(rs_epoch, var(avg_ROI_dynamic_DMN))
ylabel('Dynamic Connectivity Variance')
title('DMN')
subplot(1,3,2)
bar(rs_epoch, var(avg_ROI_dynamic_DAN))
ylabel('Dynamic Connectivity Variance')
title('DAN')
subplot(1,3,3)
bar(rs_epoch, var(avg_ROI_dynamic_SN))
ylabel('Dynamic Connectivity Variance')
title('SN')
suptitle('Connectivity Variance vs. Percent Squeeze')

%Bar graph of the mean pupil size for each resting state epoch
figure
bar(rs_epoch, mean(pupil_size_leftpp_hori));
ylabel('Mean LeftPP Hori')
title('Average Pupil Size vs. Percent Squeeze')

%Scatter plot of LC BOLD signal vs. pupil size during each resting state
%epoch
figure
for i = 1:4
    scatter(mean(pupil_size_leftpp_hori(:,i)),mean(LC_ts_all(:,i)),200,'filled','MarkerFaceColor',color_vector(i,:))
    hold on
    set(gca,'Color',[0.95 0.95 0.95])
    xlabel('Mean LeftPP Hori')
    ylabel('Mean LC Activity')
    title('LC Activity vs. Mean Pupil Size')
    legend('No Squeeze','First Squeeze','Second Squeeze','Third Squeeze','Location','best')
end

%Scatter plot of static connectivity for each network (DMN, DAN, and SN)
%vs. the mean pupil size for each resting state epoch
figure
for i = 1:4
    subplot(1,3,1)
    set(gca,'Color',[0.95 0.95 0.95])
    scatter(mean(pupil_size_leftpp_hori(:,i)),avg_ROI_static_DMN(:,i),200,'filled','MarkerFaceColor',color_vector(i,:))
    xlabel('Mean LeftPP Hori')
    ylabel('Mean Network Connectivity')
    title('DMN')
    hold on
    subplot(1,3,2)
    set(gca,'Color',[0.95 0.95 0.95])
    scatter(mean(pupil_size_leftpp_hori(:,i)),avg_ROI_static_DAN(:,i),200,'filled','MarkerFaceColor',color_vector(i,:))
    xlabel('Mean LeftPP Hori')
    ylabel('Mean Network Connectivity')
    title('DAN')
    hold on
    subplot(1,3,3)
    set(gca,'Color',[0.95 0.95 0.95])
    scatter(mean(pupil_size_leftpp_hori(:,i)),avg_ROI_static_SN(:,i),200,'filled','MarkerFaceColor',color_vector(i,:))
    xlabel('Mean LeftPP Hori')
    ylabel('Mean Network Connectivity')
    title('SN') 
    hold on
end
legend('No Squeeze','First Squeeze','Second Squeeze','Third Squeeze','Location','best')
suptitle('Static Connectivity vs. Mean Pupil Size')

%Bar graph of the variance in pupil size for each resting state epoch
figure
bar(rs_epoch, var(pupil_size_leftpp_hori));
ylabel('Variance of LeftPP Hori')
title('Variance of LeftPP Hori vs. Percent Squeeze')

%Scatter plot of the static connectivity for each network (DMN, DAN, and
%SN) vs. the variance in pupil size for each resting state epoch
figure
for i = 1:4
    subplot(1,3,1)
    set(gca,'Color',[0.95 0.95 0.95])
    scatter(var(pupil_size_leftpp_hori(:,i)),avg_ROI_static_DMN(:,i),200,'filled','MarkerFaceColor',color_vector(i,:))
    xlabel('Variance LeftPP Hori')
    ylabel('Mean Network Connectivity')
    title('DMN')
    hold on
    subplot(1,3,2)
    set(gca,'Color',[0.95 0.95 0.95])
    scatter(var(pupil_size_leftpp_hori(:,i)),avg_ROI_static_DAN(:,i),200,'filled','MarkerFaceColor',color_vector(i,:))
    xlabel('Variance LeftPP Hori')
    ylabel('Mean Network Connectivity')
    title('DAN')
    hold on
    subplot(1,3,3)
    set(gca,'Color',[0.95 0.95 0.95])
    scatter(var(pupil_size_leftpp_hori(:,i)),avg_ROI_static_SN(:,i),200,'filled','MarkerFaceColor',color_vector(i,:))
    xlabel('Variance LeftPP Hori')
    ylabel('Mean Network Connectivity')
    title('SN') 
    hold on
end
legend('No Squeeze','First Squeeze','Second Squeeze','Third Squeeze','Location','best')
suptitle('Static Connectivity vs. Variance Pupil Size')

%Scatter plot of the variance in each networks' dynamic connectivity (DMN,
%DAN, and SN) vs. the variance in pupil size for each resting state network
figure
for i = 1:4
    subplot(1,3,1)
    set(gca,'Color',[0.95 0.95 0.95])
    scatter(var(pupil_size_leftpp_hori(:,i)),var(avg_ROI_dynamic_DMN(:,i)),200,'filled','MarkerFaceColor',color_vector(i,:))
    xlabel('Variance LeftPP Hori')
    ylabel('Variance of Dynamic Network Connectivity')
    title('DMN')
    hold on
    subplot(1,3,2)
    set(gca,'Color',[0.95 0.95 0.95])
    scatter(var(pupil_size_leftpp_hori(:,i)),var(avg_ROI_dynamic_DAN(:,i)),200,'filled','MarkerFaceColor',color_vector(i,:))
    xlabel('Variance LeftPP Hori')
    ylabel('Variance of Dynamic Network Connectivity')
    title('DAN')
    hold on
    subplot(1,3,3)
    set(gca,'Color',[0.95 0.95 0.95])
    scatter(var(pupil_size_leftpp_hori(:,i)),var(avg_ROI_dynamic_SN(:,i)),200,'filled','MarkerFaceColor',color_vector(i,:))
    xlabel('Variance LeftPP Hori')
    ylabel('Variance of Dynamic Network Connectivity')
    title('SN') 
    hold on
end
legend('No Squeeze','First Squeeze','Second Squeeze','Third Squeeze','Location','best')
suptitle('Variance of Dynamic Connectivity vs. Variance of Pupil Size')

%Scatter plot of the variance in each networks' dynamic connectivity (DMN,
%DAN, and SN) vs. the mean pupil size for each resting state epoch
figure
for i = 1:4
    subplot(1,3,1)
    set(gca,'Color',[0.95 0.95 0.95])
    scatter(mean(pupil_size_leftpp_hori(:,i)),var(avg_ROI_dynamic_DMN(:,i)),200,'filled','MarkerFaceColor',color_vector(i,:))
    xlabel('Mean LeftPP Hori')
    ylabel('Mean Network Connectivity')
    title('DMN')
    hold on
    subplot(1,3,2)
    set(gca,'Color',[0.95 0.95 0.95])
    scatter(mean(pupil_size_leftpp_hori(:,i)),var(avg_ROI_dynamic_DAN(:,i)),200,'filled','MarkerFaceColor',color_vector(i,:))
    xlabel('Mean LeftPP Hori')
    ylabel('Mean Network Connectivity')
    title('DAN')
    hold on
    subplot(1,3,3)
    set(gca,'Color',[0.95 0.95 0.95])
    scatter(mean(pupil_size_leftpp_hori(:,i)),var(avg_ROI_dynamic_SN(:,i)),200,'filled','MarkerFaceColor',color_vector(i,:))
    xlabel('Mean LeftPP Hori')
    ylabel('Mean Network Connectivity')
    title('SN') 
    hold on
end
legend('No Squeeze','First Squeeze','Second Squeeze','Third Squeeze','Location','best')
suptitle('Variance in Dynamic Connectivity vs. Mean Pupil Size')


%Imagesc matrix plots of static connectiviy for each resting state epoch.
min_DMN = min(min(min(ROI_static_DMN_all,[],3)));
figure
for k = 1:4
   subplot(2,2,k)
   imagesc(ROI_static_DMN_all(:,:,k),[min_DMN,1])        
   ti = sprintf('Resting State Epoch %d.', rs_epoch(k));
   hold on
   for i = 1:ROI_DMN
       for j = 1:ROI_DMN
           if i>j
               if ROI_static_DMN_all_p(i,j,k) < 0.05
                   plot(i,j,'r*','LineWidth',2,'MarkerSize',10)
                   hold on
               end
           end
       end
   end   
   title(ti)
end
suptitle('DMN Static Conn All ROI')

min_DAN = min(min(min(ROI_static_DAN_all,[],3)));
figure
for k = 1:4
   subplot(2,2,k)
   imagesc(ROI_static_DAN_all(:,:,k),[min_DAN,1])        
   ti = sprintf('Resting State Epoch %d.', rs_epoch(k));
   hold on
   for i = 1:ROI_DAN
       for j = 1:ROI_DAN
           if i>j
               if ROI_static_DAN_all_p(i,j,k) < 0.05
                   plot(i,j,'r*','LineWidth',2,'MarkerSize',10)
                   hold on
               end
           end
       end
   end   
   title(ti)
end
suptitle('DAN Static Conn All ROI')

min_SN = min(min(min(ROI_static_SN_all,[],3)));
figure
for k = 1:4
   subplot(2,2,k)
   imagesc(ROI_static_SN_all(:,:,k),[min_SN,1])        
   ti = sprintf('Resting State Epoch %d.', rs_epoch(k));
   hold on
   for i = 1:ROI_SN
       for j = 1:ROI_SN
           if i>j
               if ROI_static_SN_all_p(i,j,k) < 0.05
                   plot(i,j,'r*','LineWidth',2,'MarkerSize',10)
                   hold on
               end
           end
       end
   end   
   title(ti)
end
suptitle('SN Static Conn All ROI')

%Scatter plot of each networks' dynamic connectivity vs. the mean pupil
%size in each of the sliding windows during each resting state epoch (for
%DMN, DAN, and SN). And plot the mean of each axis.
for i = 1:4
   figure(1)
   subplot(2,2,i)
   ti = sprintf('Squeeze %d.', rs_epoch(i));
   scatter(pupil_size_leftpp_hori_dynamics(:,i),avg_ROI_dynamic_DMN(:,i),50,'filled')
   hold on
   scatter(mean(pupil_size_leftpp_hori_dynamics(:,i)),mean(avg_ROI_dynamic_DMN(:,i)),200,'s','filled')
   xlabel('LeftPP Hori')
   ylabel('DMN Connectivity')
   title(ti)
   figure(2)
   subplot(2,2,i)
   ti = sprintf('Squeeze %d.', rs_epoch(i));
   scatter(pupil_size_leftpp_hori_dynamics(:,i),avg_ROI_dynamic_DAN(:,i),50,'filled')
   hold on
   scatter(mean(pupil_size_leftpp_hori_dynamics(:,i)),mean(avg_ROI_dynamic_DAN(:,i)),200,'s','filled')
   xlabel('LeftPP Hori')
   ylabel('DAN Connectivity')
   title(ti)
   figure(3)
   subplot(2,2,i)
   ti = sprintf('Squeeze %d.', rs_epoch(i));
   scatter(pupil_size_leftpp_hori_dynamics(:,i),avg_ROI_dynamic_SN(:,i),50,'filled')
   hold on
   scatter(mean(pupil_size_leftpp_hori_dynamics(:,i)),mean(avg_ROI_dynamic_SN(:,i)),200,'s','filled')
   xlabel('LeftPP Hori')
   ylabel('SN Connectivity')
   title(ti)    
end
suptitle('Dynamic Connectiviy Overlapping Windows')

%Line plots of the absolute values of the derivative of the dynamic
%connectivity of each network (DMN, DAN, and SN)
for i = 1:4
   figure(7)
   subplot(2,2,i)
   ti = sprintf('Resting State Epoch %d.', rs_epoch(i));
   plot(squeeze(abs(diff(avg_ROI_dynamic_DMN(:,i)))))
   xlabel('TR')
   ylabel('Abs. Deriv. of Dyn. Conn.')
   title(ti)
   suptitle('DMN Abs. Deriv of Dynamics Conn.')
   figure(8)
   subplot(2,2,i)
   ti = sprintf('Resting State Epoch %d.', rs_epoch(i));
   plot(squeeze(abs(diff(avg_ROI_dynamic_DAN(:,i)))))
   xlabel('TR')
   ylabel('Abs. Deriv. of Dyn. Conn.')
   title(ti)
   suptitle('DAN Abs. Deriv of Dynamics Conn.')
   figure(9)
   subplot(2,2,i)
   ti = sprintf('Resting State Epoch %d.', rs_epoch(i));
   plot(squeeze(abs(diff(avg_ROI_dynamic_SN(:,i)))))
   xlabel('TR')
   ylabel('Abs. Deriv. of Dyn. Conn.')
   title(ti)
   suptitle('SN Abs. Deriv of Dynamics Conn.')  
    
end
