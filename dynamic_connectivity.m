function [output] = dynamic_connectivity(num_ROI, ROI_ts, rs_blocks_only, rs_tr, delta_t, LC_ts_rs_blocks, nonoverlap_RS_all, rs_time_only)

%Overlapping Windows

for b = 1:rs_blocks_only
    for t = 1:rs_tr-delta_t
        for r1 = 1:num_ROI
            for r2 = 1:num_ROI
                temp = rs_time_only{b,1};
                dyn_corr_val = corrcoef(squeeze(ROI_ts(temp(t:t+delta_t),r1)),squeeze(ROI_ts(temp(t:t+delta_t),r2)));
                ROI_dynamic_all_blocks(r1,r2,t,b) = dyn_corr_val(2,1);
            end
        end
    end
end

for b = 1:rs_blocks_only
    for k = 1:rs_tr-delta_t

        count = 1;
        for r1 = 2:num_ROI
            for r2 = r1-1
                output(count,1) = squeeze(ROI_dynamic_all_blocks(r1,r2,k,b));
                count = count + 1;
            end
        end
        avg_ROI_dynamic(k,b) = mean(output(:,1));
        clear output
        clear count
    end
end

%Correlate LC with the networks using overlapping windows
for b = 1:rs_blocks_only
    for t = 1:rs_tr-2*delta_t
        corr_val = corrcoef(squeeze(avg_ROI_dynamic(t:t+delta_t,b)),squeeze(LC_ts_rs_blocks(t:t+delta_t,b)));
        LC_dynamic(t,b) = corr_val(2,1);

    end
end

length_vector = size(nonoverlap_RS_all);
%Nonoverlapping Windows
for b = 1:rs_blocks_only
    for t = 1:length_vector(1)-1
        for r1 = 1:num_ROI
            for r2 = 1:num_ROI
                clear corr_val
                v = nonoverlap_RS_all(t,b);
                vv = nonoverlap_RS_all(t+1,b);
                [corr_val, p_val] = corrcoef(squeeze(ROI_ts(v:vv,r1)),squeeze(ROI_ts(v:vv,r2)));
                inter = corr_val(2,1);
                ROI_dynamic_nonoverlap_all_blocks(r1,r2,t,b) = inter;
                ROI_dynamic_nonoverlap_all_blocks_p(r1,r2,t,b) = p_val(2,1);

            end
        end
    end
end

for b = 1:rs_blocks_only
    count = 1;
    for r1 = 2:num_ROI
       for r2 = r1-1
            output(count,1) = squeeze(ROI_dynamic_nonoverlap_all_blocks(r1,r2,t,b));
            count = count + 1;
            end
        end
        avg_ROI_dynamic_nonoverlap(t,b) = mean(output(:,1));
        clear output
        clear count
end


output.ROI_dynamic_all_blocks = ROI_dynamic_all_blocks;
output.avg_ROI_dynamic = avg_ROI_dynamic;
output.LC_dynamic = LC_dynamic;
output.ROI_dynamic_nonoverlap_all_blocks = ROI_dynamic_nonoverlap_all_blocks;
output.ROI_dynamic_nonoverlap_all_blocks_p = ROI_dynamic_nonoverlap_all_blocks_p;
output.avg_ROI_dynamic_nonoverlap = avg_ROI_dynamic_nonoverlap;

end
