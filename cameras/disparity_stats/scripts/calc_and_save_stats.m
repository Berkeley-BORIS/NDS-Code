function [] = calc_and_save_stats(data_mat_path)

display('loading');
tic
data = load(data_mat_path);
toc
display('calc');
tic
for j = 1:207
    for k = 1:207
        quartile1_mat(j,k) = quantile(data.disparity(j,k,~isnan(data.disparity(j,k,:))),0.25,3);
        quartile3_mat(j,k) = quantile(data.disparity(j,k,~isnan(data.disparity(j,k,:))),0.75,3);
    end
end
%quartile1_mat = quantile(data.disparity(~isnan(data.disparity)),0.25,3);
median_mat = nanmedian(data.disparity,3);
%quartile3_mat = quantile(data.disparity(~isnan(data.disparity)),0.75,3);
var_mat = quartile3_mat - quartile1_mat;
count_mat = sum(~isnan(data.disparity),3);
percent_mat = sum(~isnan(data.disparity),3)./size(data.disparity,3);

mean_mat = nanmean(data.disparity,3);
std_mat = nanstd(data.disparity,0,3);
sterr_mat = std_mat./sqrt(count_mat);

toc
display('saving');
tic
clear data.disparity;
save(strrep(data_mat_path,'.mat','_stats.mat'),'quartile1_mat','median_mat','quartile3_mat','var_mat','count_mat','percent_mat','mean_mat','std_mat','sterr_mat');
toc