function [output] = LC_stuffz(LC_mask_left, LC_mask_right, timeseries, total_blocks, all_time, rs_blocks_only, rs_time_only)

%In this section of code, we are extracting the time series for the LC
%(that is, for the LC ROI)
clear LC_left_ts;
clear LC_right_ts;
[nx, ny, nz, nt] = size(timeseries);
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


for b = 1:total_blocks  
    LC_ts_mean(b,:) = mean(LC_ts(all_time{b,1}));
    LC_ts_std(b,:) = std(LC_ts(all_time{b,1}));
end

for i = 1:rs_blocks_only
    LC_ts_rs_blocks(:,i) = LC_ts(rs_time_only{i,1});
end




output.LC_ts = LC_ts;
output.LC_ts_rs_blocks = LC_ts_rs_blocks;
output.LC_ts_mean = LC_ts_mean;
output.LC_ts_std = LC_ts_std;

end
