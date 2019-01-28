%This code utilizes the connectivity matrices of networks (from
%Bayesian_analysis.m) to perform Bayes' rule and determine the probability
%of whether the subject squeezed a squeeze-ball.

%before = within-network connectivity values before squeezing
%after  = within-network connectivity values after squeezing
%bw     = arbitrary bandwith (values acquired online)
function [output] = Bayesian_analysis(before, after,bw)

size_before = size(before);
size_after = size(after);

rounding = [before after];
min_rounding = round(min(min(rounding)),2);
max_rounding = round(max(max(rounding)),2) + 0.01;

step_size = 0.1; 


gridx = min_rounding:step_size:max_rounding;


%Because I am investigating whether the connectivity pattern of a network
%in conjunction with one or more networks could be more informative than
%that of a single network alone, I use the switch/case operators to
%determine if I am examining 2 (case 2), 3 (case 3), or 4 (case 4) networks
%at once.

switch size_before(2)   
        
    case 2
        %To create a probability density function, I employ multivariable kernel
        %density estimation by creating a grid with as many dimensions as networks
        %being examined. The size of this grid is proportional to the spread of the
        %correlation coefficients in the networks. 
        [x1,x2] = ndgrid(gridx,gridx);
        x1 = x1(:,:)';
        x2 = x2(:,:)';
        xi = [x1(:) x2(:)];
        size_xi = size(xi);
        [f_before,~] = mvksdensity(before,xi,'Bandwidth',bw);%, 'Function','pdf');
        f_before = f_before./(sum(f_before));
        [f_after,~] = mvksdensity(after,xi,'Bandwidth',bw);% 'Kernal','normpdf');
        f_after = f_after./(sum(f_after));
        for k = 1:size_before(1)
            
            %The data points (the connectivity profile of the networks) fall within this grid, through not
            %exactly at the intersections. Thus, the Euclidean distances between the
            %correlation coefficients and the intersection of the grid are calculated
            %in effort to find the minimum distance between the two. Then
            %the height at the grid intersection closest to the data point
            %is observed because this height corresponds to the probability that the connectivity
            %pattern occurred before squeezing, and the probability that it
            %occurred after.
            
            % Find the probability of each "before" data GIVEN before AND after
            % squeezing - 90 data points
            euclidean_distance_before = sqrt((xi(:,1) - before(k,1)).^2 + (xi(:,2) - before(k,2)).^2);
            [m_before(k),ID_before(k)] = min(euclidean_distance_before);
            locations_before(k,:) = xi(ID_before(k),:); %these are the locations on the grid where we want to look at
            probability_before_before(k,:) = f_before(ID_before(k)); % p(d_before | before) Probability density of "before" function at the "before" data points
            probability_before_after(k,:) = f_after(ID_before(k)); % p(d_before | after) Probability densith of "after" function at the "before" data points
    
            euclidean_distance_after = sqrt((xi(:,1) - after(k,1)).^2 + (xi(:,2) - after(k,2)).^2);
            % Find the probability of each "after" data GIVEN before AND after
            % squeezing - 90 data points
            [m_after(k),ID_after(k)] = min(euclidean_distance_after);
            locations_after(k,:) = xi(ID_after(k),:); %these are the locations on the grid where we want to look at
            probability_after_before(k,:) = f_before(ID_after(k)); % p(d_after | before) Probability density of the "before" function at the "after" data points
            probability_after_after(k,:) = f_after(ID_after(k)); % p(d_after | after) Probability of the "after" function at the "after" data points

        end
        
        
        
     
    case 3
        [x1,x2,x3] = ndgrid(gridx,gridx,gridx);
        x1 = x1(:,:)';
        x2 = x2(:,:)';
        x3 = x3(:,:)';
        xi = [x1(:) x2(:) x3(:)];
        size_xi = size(xi);
        [f_before,~] = mvksdensity(before,xi,'Bandwidth',bw);%, 'Function','pdf');
        f_before = f_before./(sum(f_before));
        [f_after,~] = mvksdensity(after,xi,'Bandwidth',bw);% 'Kernal','normpdf');
        f_after = f_after./(sum(f_after));
        for k = 1:size_before(1)            
            % Find the probability of each "before" data GIVEN before AND after
            % squeezing - 90 data points
            euclidean_distance_before = sqrt((xi(:,1) - before(k,1)).^2 + (xi(:,2) - before(k,2)).^2 + (xi(:,3) - before(k,3)).^2);
            [m_before(k),ID_before(k)] = min(euclidean_distance_before);
            locations_before(k,:) = xi(ID_before(k),:); %these are the locations on the grid where we want to look at
            probability_before_before(k,:) = f_before(ID_before(k)); % p(d_before | before) Probability density of "before" function at the "before" data points
            probability_before_after(k,:) = f_after(ID_before(k)); % p(d_before | after) Probability densith of "after" function at the "before" data points
    
            euclidean_distance_after = sqrt((xi(:,1) - after(k,1)).^2 + (xi(:,2) - after(k,2)).^2 + (xi(:,3) - after(k,3)).^2);
            % Find the probability of each "after" data GIVEN before AND after
            % squeezing - 90 data points
            [m_after(k),ID_after(k)] = min(euclidean_distance_after);
            locations_after(k,:) = xi(ID_after(k),:); %these are the locations on the grid where we want to look at 
            probability_after_before(k,:) = f_before(ID_after(k)); % p(d_after | before) Probability density of the "before" function at the "after" data points
            probability_after_after(k,:) = f_after(ID_after(k)); % p(d_after | after) Probability of the "after" function at the "after" data points


        end
        
        
    case 4
        [x1,x2,x3,x4] = ndgrid(gridx,gridx,gridx,gridx);
        x1 = x1(:,:)';
        x2 = x2(:,:)';
        x3 = x3(:,:)';
        x4 = x4(:,:)';
        xi = [x1(:) x2(:) x3(:) x4(:)];
        size_xi = size(xi);
        [f_before,~] = mvksdensity(before,xi,'Bandwidth',bw);%, 'Function','pdf');
        f_before = f_before./(sum(f_before));
        [f_after,~] = mvksdensity(after,xi,'Bandwidth',bw);% 'Kernal','normpdf');
        f_after = f_after./(sum(f_after));
        
        for k = 1:size_before(1)
            % Find the probability of each "before" data GIVEN before AND after
            % squeezing - 90 data points
            euclidean_distance_before = sqrt((xi(:,1) - before(k,1)).^2 + (xi(:,2) - before(k,2)).^2 + (xi(:,3) - before(k,3)).^2 + (xi(:,4) - after(k,4)).^2);
            [m_before(k),ID_before(k)] = min(euclidean_distance_before);
            locations_before(k,:) = xi(ID_before(k),:); %these are the locations on the grid where we want to look at 
            probability_before_before(k,:) = f_before(ID_before(k)); % p(d_before | before) Probability density of "before" function at the "before" data points
            probability_before_after(k,:) = f_after(ID_before(k)); % p(d_before | after) Probability densith of "after" function at the "before" data points
    
            euclidean_distance_after = sqrt((xi(:,1) - after(k,1)).^2 + (xi(:,2) - after(k,2)).^2 + (xi(:,3) - after(k,3)).^2 + (xi(:,4) - after(k,4)).^2);
            % Find the probability of each "after" data GIVEN before AND after
            % squeezing - 90 data points
            [m_after(k),ID_after(k)] = min(euclidean_distance_after);
            locations_after(k,:) = xi(ID_after(k),:); %these are the locations on the grid where we want to look at
            probability_after_before(k,:) = f_before(ID_after(k)); % p(d_after | before) Probability density of the "before" function at the "after" data points
            probability_after_after(k,:) = f_after(ID_after(k)); % p(d_after | after) Probability of the "after" function at the "after" data points


        end
        
end

%Concatenate the before and after probabilities
p_allDots_before = [probability_before_before; probability_after_before];
p_allDots_after = [probability_before_after; probability_after_after];

%Now I perform Bayes' rule
%The prior probability is assumed to be 0.5 because there is an equal
%probability that the connectivity occurred before or after the squeeze.
Bayes_Rule = p_allDots_after ./ (p_allDots_before + p_allDots_after); % p(after | allDots)
labels = [zeros(90,1); ones(90,1)]; %Before squeezing is 0 and after squeezing is 1

%Acquire the x values, y values, and AUC value for an ROC curve
[perfcurve_x,perfcurve_y,~,perfcurve_AUC] = perfcurve(labels,Bayes_Rule,1);

%Create outputs for this function
output.Bayes_Rule = Bayes_Rule;
output.perfcurve_x = perfcurve_x;
output.perfcurve_y = perfcurve_y;
output.perfcurve_AUC = perfcurve_AUC;
output.size_before = size_before;
output.euclidean_distance_before = euclidean_distance_before;



end