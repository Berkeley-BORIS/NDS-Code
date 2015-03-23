function [] = calc_weighted_medians()

meridian = {'VM','HM'};

version = {'corrected'};

tasks = {'nature_walk_1','ordering_coffee','sandwich'};
weights = [0.2 0.51 0.29];

weighted_medians = zeros(1,207);
cis = zeros(2,207);

for m = 1:length(meridian)
    for v = 1:length(version)

        load(['/Users/natdispstats/Documents/cameras/disparity_stats/' meridian{m} '_dispdata_' version{v} '.mat']);
        dat = [nature_walk_1 ordering_coffee sandwich];
        weights = [repmat(weights(1)/size(nature_walk_1,2),size(nature_walk_1,1),size(nature_walk_1,2)) ...
            repmat(weights(2)/size(ordering_coffee,2),size(ordering_coffee,1),size(ordering_coffee,2)) ...
            repmat(weights(3)/size(sandwich,2),size(sandwich,1),size(sandwich,2))];
        
        for ecc = 1:size(dat,1)
        %for ecc = 100
            display(ecc)
            bools = ~isnan(dat(ecc,:));
            indices = find(~isnan(dat(ecc,:)));
            if sum(bools) >= 3000;
                weighted_medians(ecc) = weighted_median_func(dat(ecc,indices),weights(ecc,indices));
                
                %bootmedian = @(x)(weighted_median_func(x,weights(ecc,indices)));
                bootmedian = @(x)(weighted_median_func(dat(ecc,x),weights(ecc,x)));
                
                %ci = bootci(1000,bootmedian,dat(ecc,indices));
                ci = bootci(1000,bootmedian,indices);
                cis(:,ecc) = ci;
            else
                weighted_medians(ecc) = NaN;
                cis(:,ecc) = [NaN NaN];
            end
        end
        
        %save(['/Users/natdispstats/Documents/cameras/disparity_stats/' meridian{m} '_dispdata_' version{v} '_weighted_median_and_ci.mat'],'weighted_medians','cis');
        
    end
end

keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [wmed] = weighted_median_func(x,w)

[sortx,order] = sort(x);
sortw = w(order);

midpoint = sum(sortw)/2;
csumw = cumsum(sortw);
j = find(csumw<=midpoint,1,'last');
if csumw(j) == midpoint
     wmed = mean(sortx([j j+1]));
else
     wmed = sortx(j+1);
end