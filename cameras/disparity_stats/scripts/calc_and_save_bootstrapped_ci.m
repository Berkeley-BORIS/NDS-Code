function [] = calc_and_save_bootstrapped_ci(data_mat_path)

display('loading');
tic
data = load(data_mat_path);
toc
display('calc');

median_boot = zeros(1,1000);
sterrboot_mat = NaN*ones(207,207);

median_mat = nanmedian(data.disparity,3);
mat_size = size(data.disparity,3);
tic
for j = 1:207
    display(num2str(j));
    for k = 1:207
        if sum(~isnan(data.disparity(j,k,:))) ~= 0
            cnt = 1;
            for x = 1:1000
                sample = randsample(reshape(data.disparity(j,k,:),1,mat_size),mat_size,true);
                median_boot(cnt) = nanmedian(sample);
                cnt = cnt + 1;
            end
            sterrboot_mat(j,k) = std(median_boot);
        end
    end
end

toc
display('saving');
tic
clear data.disparity;
save(strrep(data_mat_path,'.mat','_bootstrapped_ci.mat'),'median_mat','sterrboot_mat');
toc