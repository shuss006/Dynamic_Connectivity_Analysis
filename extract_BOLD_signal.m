function [output] = extract_BOLD_signal(num_ROI, ROI_mask, timeseries, total_blocks, all_time)

[nx, ny, nz, nt] = size(timeseries);
ROI_ts=zeros(nt,num_ROI);

%Extract the BOLD signal from each voxel
for r = 1:num_ROI %Iterating over all ROIs
    count = 1;
    clear p;
    for x = 1:nx %Iterating over all values of x
        for y = 1:ny %Iterating over all values of y
            for z = 1:nz %Iterating over all values of z
                %Extract 1 ROI for the entire time series
                if(ROI_mask(x,y,z) == r) %Finding all of the ROIs in the time series
                    p(:,count) = squeeze(timeseries(x,y,z,:)); %Taking all the values in the time series for a each ROI
                    count = count+1;
                end
            end
        end
    end
    %Average the dummy variable, a
    for l = 1:nt
        ROI_ts(l,r) = mean(p(l,:));
    end
   clear p;
end

%Bandpass Filter
for r = 1:num_ROI
    ROI_ts(:,r) = rest_filter(2,[0.01 0.1],ROI_ts(:,r));
end


%Average the time series to get a single BOLD signal representing the
%entire network
for k = 1:nt
    count = 1;
    for r1 = 1:num_ROI
        output(count,1) = squeeze(ROI_ts(k,r1));
        count = count+1;
    end
    
    ROI_ts_avg(k,1) = mean(output(:,1));

    clear output
    clear count
end

%Find the mean and std of each block of the scan
for b = 1:total_blocks
    ROI_ts_mean(b,:) = mean(ROI_ts(all_time{b,1}));
    ROI_ts_std(b,:) = std(ROI_ts(all_time{b,1}));
end


output.ROI_ts = ROI_ts;
output.ROI_ts_avg = ROI_ts_avg;
output.ROI_ts_mean = ROI_ts_mean;
output.ROI_ts_std = ROI_ts_std;


end






