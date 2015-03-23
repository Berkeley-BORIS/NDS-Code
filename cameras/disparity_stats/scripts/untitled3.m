load('mean_vert_disp.mat')

eccentricity_deg = 0.0983*(-103:103);
deg_8 = find(eccentricity_deg <=8 & eccentricity_deg >=-8);

nat1 = load('../data/bwsnat5/bwsnat5_task_nature_walk_1_disparity.mat');
nat1.disparity_vert = nat1.disparity(deg_8,104,:);

nat1.disparity_vert_norm = zeros(163,size(nat1.disparity_vert,3));
for j = 1:size(nat1.disparity_vert,3)
    nat1.disparity_vert_norm(:,j) = nat1.disparity_vert(:,1,j) - nat1.disparity_vert(82,:,j);
end

nat1.disparity_vert_norm_median = nanmedian(nat1.disparity_vert_norm,2);

figure(); hold on;
scatter(reshape(nat1.disparity_vert_norm,1,163*3502),repmat(vd_ecc',3502,1),'ko')
plot(nat1.disparity_vert_norm_median,vd_ecc','r-')