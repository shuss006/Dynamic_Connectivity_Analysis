function [output] = static_correlation(num_ROI,ROI_ts, total_blocks, all_time, LC_ts, rs_blocks_only)

%Correlate all ROIs within a network
for r1 = 1:num_ROI
    for r2 = 1:num_ROI
        for b = 1:total_blocks
            [corr_ouput, p_val_output] = corrcoef(squeeze(ROI_ts(all_time{b,1},r1)),squeeze(ROI_ts(all_time{b,1},r2)));
            ROI_static_everything(r1,r2,b) = corr_ouput(2,1);
            ROI_static_everything_p(r1,r2,b) = p_val_output(2,1);
            
        end
    end
end

for b = 1:total_blocks
    count = 1;
    for r1 = 1:num_ROI
        for r2 = 1:num_ROI
            if r2 > r1
                ROI_static_all(count,b) = ROI_static_everything(r1,r2,b);
                ROI_static_all_p(count,b) = ROI_static_everything_p(r1,r2,b);                
                count = count +1; 
            end
        end
    end
    clear count
end

avg_ROI_static_all = mean(ROI_static_all);
avg_ROI_static_Fish_z = 0.5*(log(1+avg_ROI_static_all)-log(1-avg_ROI_static_all));

%Static Correlations of Network with LC
for r1 = 1:num_ROI
    for t = 1:total_blocks
       [corr_val, stat_p_val] = corrcoef(squeeze(ROI_ts(all_time{t,1},r1)),squeeze(LC_ts(all_time{t,1},1)));    
       LC_static_all(r1,t) = corr_val(2,1);
       LC_static_all_p(r1,t) = stat_p_val(2,1);
    end
end

%Spearman Rank Correlation of all RS Epochs Compared with No Squeeze
%Condition
for i = 1:rs_blocks_only
   [spear_corr(:,:,i), spear_corr_p(:,:,i)] = corr(ROI_static_everything(:,:,1),ROI_static_everything(:,:,i),'Type','Spearman');
end


output.avg_ROI_static_all = avg_ROI_static_all;
output.ROI_static_everything = ROI_static_everything;
output.ROI_static_everything_p = ROI_static_everything_p;
output.ROI_static_all = ROI_static_all;
output.ROI_static_all_p = ROI_static_all_p;
output.avg_ROI_static_all = avg_ROI_static_all;
output.LC_static_all = LC_static_all;
output.LC_static_all_p = LC_static_all_p;
output.spear_corr = spear_corr;
output.spear_corr_p = spear_corr_p;
output.avg_ROI_static_Fish_z = avg_ROI_static_Fish_z;



end
