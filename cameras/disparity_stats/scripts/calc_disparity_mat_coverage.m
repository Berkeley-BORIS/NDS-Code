function [] = calc_disparity_mat_coverage(filename)

subjname = strtok(filename,'_');
loadpath = ['../data/' subjname '/' filename];
savepath = ['../data/' subjname '/' strtok(filename,'.') '_coverage.mat'];
disparity = load(loadpath);

disparity = disparity.disparity;

coverage_mat = zeros(size(disparity,3),5);

%total number of pixels in filled 2,4,6,8,10deg circle?
degperpixel = 2*atand(.5/583);

total = [];
masks = [];
cnt = 1;
for rad = ceil([2:2:10]./degperpixel)
    mat = zeros(size(disparity(:,:,1),1),size(disparity(:,:,1),1));
    for j = 1:length(mat)
        for k = 1:length(mat)
            if ((j-ceil(length(mat)/2)).^2 + (k-ceil(length(mat)/2)).^2) <= rad.^2
                mat(j,k) = 1;
            end
        end
    end
    total(cnt) = sum(sum(mat));
    masks(:,:,cnt) = mat;
    cnt = cnt+1;
end

masks = logical(masks);

for r = 1:size(masks,3)
    display(r)
    for m = 1:size(disparity,3)
        tmp = disparity(:,:,m);
        coverage_mat(m,r) = (sum(sum(~isnan(tmp((masks(:,:,r))))))/total(r)).*100;    
    end
end

save(savepath,'coverage_mat');
figure(); hold on; hist(coverage_mat);

