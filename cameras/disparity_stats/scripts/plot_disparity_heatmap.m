bws_inside_walk_1 = load('../data/bwsdrdp6/bwsdrdp6_task_walk_1_disparity_stats.mat');
bws_inside_walk_2 = load('../data/bwsdrdp6/bwsdrdp6_task_walk_2_disparity_stats.mat');
bws_inside_walk_3 = load('../data/bwsdrdp6/bwsdrdp6_task_walk_3_disparity_stats.mat');
bws_nature_walk_1 = load('../data/bwsnat5/bwsnat5_task_nature_walk_1_disparity_stats.mat');
bws_nature_walk_2 = load('../data/bwsure1/bwsure1_task_nature_walk_2_disparity_stats.mat');
bws_campus_walk = load('../data/bwsnat5/bwsnat5_task_campus_walk_3_disparity_stats.mat');
bws_cafe = load('../data/bwscafe1/bwscafe1_task_ordering_coffee_disparity_stats.mat');
bws_sandwich = load('../data/bwssand1/bwssand1_task_sandwich_disparity_stats.mat');

% bws_inside_walk_1_ci = load('../data/bwsdrdp6/bwsdrdp6_task_walk_1_disparity_bootstrapped_ci.mat');
% bws_inside_walk_2_ci = load('../data/bwsdrdp6/bwsdrdp6_task_walk_2_disparity_bootstrapped_ci.mat');
% bws_inside_walk_3_ci = load('../data/bwsdrdp6/bwsdrdp6_task_walk_3_disparity_bootstrapped_ci.mat');
% bws_nature_walk_1_ci = load('../data/bwsnat5/bwsnat5_task_nature_walk_1_disparity_bootstrapped_ci.mat');
% bws_nature_walk_2_ci = load('../data/bwsure1/bwsure1_task_nature_walk_2_disparity_bootstrapped_ci.mat');
% bws_campus_walk_ci = load('../data/bwsnat5/bwsnat5_task_campus_walk_3_disparity_bootstrapped_ci.mat');
% bws_cafe_ci = load('../data/bwscafe1/bwscafe1_task_ordering_coffee_disparity_bootstrapped_ci.mat');
% bws_sandwich_ci = load('../data/bwssand1/bwssand1_task_sandwich_disparity_bootstrapped_ci.mat');

ges_inside_walk_1 = load('../data/gesdrdp1/gesdrdp1_task_walk_1_disparity_stats.mat');
ges_inside_walk_2 = load('../data/gesdrdp1/gesdrdp1_task_walk_2_disparity_stats.mat');
ges_inside_walk_3 = load('../data/gesdrdp1/gesdrdp1_task_walk_3_disparity_stats.mat');
ges_nature_walk_1 = load('../data/gesnat1/gesnat1_task_nature_walk_1_disparity_stats.mat');
ges_nature_walk_2 = load('../data/gesure1/gesure1_task_nature_walk_2_disparity_stats.mat');
ges_campus_walk = load('../data/gesnat1/gesnat1_task_campus_walk_3_disparity_stats.mat');
ges_cafe = load('../data/gescafe2/gescafe2_task_ordering_coffee_disparity_stats.mat');
ges_sandwich = load('../data/gessand1/gessand1_task_sandwich_disparity_stats.mat');

% ges_inside_walk_1_ci = load('../data/gesdrdp1/gesdrdp1_task_walk_1_disparity_bootstrapped_ci.mat');
% ges_inside_walk_2_ci = load('../data/gesdrdp1/gesdrdp1_task_walk_2_disparity_bootstrapped_ci.mat');
% ges_inside_walk_3_ci = load('../data/gesdrdp1/gesdrdp1_task_walk_3_disparity_bootstrapped_ci.mat');
% ges_nature_walk_1_ci = load('../data/gesnat1/gesnat1_task_nature_walk_1_disparity_bootstrapped_ci.mat');
% ges_nature_walk_2_ci = load('../data/gesure1/gesure1_task_nature_walk_2_disparity_bootstrapped_ci.mat');
% ges_campus_walk_ci = load('../data/gesnat1/gesnat1_task_campus_walk_3_disparity_bootstrapped_ci.mat');
% ges_cafe_ci = load('../data/gescafe2/gescafe2_task_ordering_coffee_disparity_bootstrapped_ci.mat');
% ges_sandwich_ci = load('../data/gessand1/gessand1_task_sandwich_disparity_bootstrapped_ci.mat');

hmf_inside_walk_1 = load('../data/hmfdrdp1/hmfdrdp1_task_walk_1_disparity_stats.mat');
hmf_inside_walk_2 = load('../data/hmfdrdp1/hmfdrdp1_task_walk_2_disparity_stats.mat');
hmf_inside_walk_3 = load('../data/hmfdrdp1/hmfdrdp1_task_walk_3_disparity_stats.mat');
hmf_nature_walk_1 = load('../data/hmfnat1/hmfnat1_task_nature_walk_1_disparity_stats.mat');
hmf_nature_walk_2 = load('../data/hmfure1/hmfure1_task_nature_walk_2_disparity_stats.mat');
hmf_campus_walk = load('../data/hmfnat1/hmfnat1_task_campus_walk_3_disparity_stats.mat');
hmf_cafe = load('../data/hmfcafe1/hmfcafe1_task_ordering_coffee_disparity_stats.mat');

% hmf_nature_walk_1_ci = load('../data/hmfnat1/hmfnat1_task_nature_walk_1_disparity_bootstrapped_ci.mat');
% hmf_campus_walk_ci = load('../data/hmfnat1/hmfnat1_task_campus_walk_3_disparity_bootstrapped_ci.mat');

%enforce minimum data count
min_count = 1000;
min_count_walks_all = 4;
min_count_other_all = 2;

bws_inside_walk_1.median_mat(bws_inside_walk_1.count_mat < min_count) = NaN;
bws_inside_walk_2.median_mat(bws_inside_walk_2.count_mat < min_count) = NaN;
bws_inside_walk_3.median_mat(bws_inside_walk_3.count_mat < min_count) = NaN;
bws_nature_walk_1.median_mat(bws_nature_walk_1.count_mat < min_count) = NaN;
bws_nature_walk_2.median_mat(bws_nature_walk_2.count_mat < min_count) = NaN;
bws_campus_walk.median_mat(bws_campus_walk.count_mat < min_count) = NaN;
bws_cafe.median_mat(bws_cafe.count_mat < min_count) = NaN;
bws_sandwich.median_mat(bws_sandwich.count_mat < min_count) = NaN;

ges_inside_walk_1.median_mat(ges_inside_walk_1.count_mat < min_count) = NaN;
ges_inside_walk_2.median_mat(ges_inside_walk_2.count_mat < min_count) = NaN;
ges_inside_walk_3.median_mat(ges_inside_walk_3.count_mat < min_count) = NaN;
ges_nature_walk_1.median_mat(ges_nature_walk_1.count_mat < min_count) = NaN;
ges_nature_walk_2.median_mat(ges_nature_walk_2.count_mat < min_count) = NaN;
ges_campus_walk.median_mat(ges_campus_walk.count_mat < min_count) = NaN;
ges_cafe.median_mat(ges_cafe.count_mat < min_count) = NaN;
ges_sandwich.median_mat(ges_sandwich.count_mat < min_count) = NaN;

hmf_inside_walk_1.median_mat(hmf_inside_walk_1.count_mat < min_count) = NaN;
hmf_inside_walk_2.median_mat(hmf_inside_walk_2.count_mat < min_count) = NaN;
hmf_inside_walk_3.median_mat(hmf_inside_walk_3.count_mat < min_count) = NaN;
hmf_nature_walk_1.median_mat(hmf_nature_walk_1.count_mat < min_count) = NaN;
hmf_nature_walk_2.median_mat(hmf_nature_walk_2.count_mat < min_count) = NaN;
hmf_campus_walk.median_mat(hmf_campus_walk.count_mat < min_count) = NaN;
hmf_cafe.median_mat(hmf_cafe.count_mat < min_count) = NaN;

all_inside_walks_count = nansum(cat(3,...
    bws_inside_walk_1.count_mat > min_count,bws_inside_walk_2.count_mat > min_count,bws_inside_walk_3.count_mat > min_count,...
    ges_inside_walk_1.count_mat > min_count,ges_inside_walk_2.count_mat > min_count,ges_inside_walk_3.count_mat > min_count,...
    hmf_inside_walk_1.count_mat > min_count,hmf_inside_walk_2.count_mat > min_count,hmf_inside_walk_3.count_mat > min_count),3);

all_outside_walks_count = nansum(cat(3,...
    bws_nature_walk_1.count_mat > min_count,bws_nature_walk_2.count_mat > min_count,bws_campus_walk.count_mat > min_count,...
    ges_nature_walk_1.count_mat > min_count,ges_nature_walk_2.count_mat > min_count,ges_campus_walk.count_mat > min_count,...
    hmf_nature_walk_1.count_mat > min_count,hmf_nature_walk_2.count_mat > min_count,hmf_campus_walk.count_mat > min_count),3);

all_other_tasks_count = nansum(cat(3,...
    bws_cafe.count_mat > min_count,bws_sandwich.count_mat > min_count,...
    ges_cafe.count_mat > min_count,ges_sandwich.count_mat > min_count,...
    hmf_cafe.count_mat > min_count),3);


%normalize
bws_inside_walk_1.median_mat = bws_inside_walk_1.median_mat - bws_inside_walk_1.median_mat(104,104);
bws_inside_walk_2.median_mat = bws_inside_walk_2.median_mat - bws_inside_walk_2.median_mat(104,104);
bws_inside_walk_3.median_mat = bws_inside_walk_3.median_mat - bws_inside_walk_3.median_mat(104,104);
bws_nature_walk_1.median_mat = bws_nature_walk_1.median_mat - bws_nature_walk_1.median_mat(104,104);
bws_nature_walk_2.median_mat = bws_nature_walk_2.median_mat - bws_nature_walk_2.median_mat(104,104);
bws_campus_walk.median_mat = bws_campus_walk.median_mat - bws_campus_walk.median_mat(104,104);
bws_cafe.median_mat = bws_cafe.median_mat - bws_cafe.median_mat(104,104);
bws_sandwich.median_mat = bws_sandwich.median_mat - bws_sandwich.median_mat(104,104);

ges_inside_walk_1.median_mat = ges_inside_walk_1.median_mat - ges_inside_walk_1.median_mat(104,104);
ges_inside_walk_2.median_mat = ges_inside_walk_2.median_mat - ges_inside_walk_2.median_mat(104,104);
ges_inside_walk_3.median_mat = ges_inside_walk_3.median_mat - ges_inside_walk_3.median_mat(104,104);
ges_nature_walk_1.median_mat = ges_nature_walk_1.median_mat - ges_nature_walk_1.median_mat(104,104);
ges_nature_walk_2.median_mat = ges_nature_walk_2.median_mat - ges_nature_walk_2.median_mat(104,104);
ges_campus_walk.median_mat = ges_campus_walk.median_mat - ges_campus_walk.median_mat(104,104);
ges_cafe.median_mat = ges_cafe.median_mat - ges_cafe.median_mat(104,104);
ges_sandwich.median_mat = ges_sandwich.median_mat - ges_sandwich.median_mat(104,104);

hmf_inside_walk_1.median_mat = hmf_inside_walk_1.median_mat - hmf_inside_walk_1.median_mat(104,104);
hmf_inside_walk_2.median_mat = hmf_inside_walk_2.median_mat - hmf_inside_walk_2.median_mat(104,104);
hmf_inside_walk_3.median_mat = hmf_inside_walk_3.median_mat - hmf_inside_walk_3.median_mat(104,104);
hmf_nature_walk_1.median_mat = hmf_nature_walk_1.median_mat - hmf_nature_walk_1.median_mat(104,104);
hmf_nature_walk_2.median_mat = hmf_nature_walk_2.median_mat - hmf_nature_walk_2.median_mat(104,104);
hmf_campus_walk.median_mat = hmf_campus_walk.median_mat - hmf_campus_walk.median_mat(104,104);
hmf_cafe.median_mat = hmf_cafe.median_mat - hmf_cafe.median_mat(104,104);


%average across subjects
all_inside_walks_mean = nanmean(cat(3,...
    bws_inside_walk_1.median_mat,bws_inside_walk_2.median_mat,bws_inside_walk_3.median_mat,...
    ges_inside_walk_1.median_mat,ges_inside_walk_2.median_mat,ges_inside_walk_3.median_mat,...
    hmf_inside_walk_1.median_mat,hmf_inside_walk_2.median_mat,hmf_inside_walk_3.median_mat),3);
all_outside_walks_mean = nanmean(cat(3,...
    bws_nature_walk_1.median_mat,bws_nature_walk_2.median_mat,bws_campus_walk.median_mat,...
    ges_nature_walk_1.median_mat,ges_nature_walk_2.median_mat,ges_campus_walk.median_mat,...
    hmf_nature_walk_1.median_mat,hmf_nature_walk_2.median_mat,hmf_campus_walk.median_mat),3);
all_other_mean = nanmean(cat(3,...
    bws_cafe.median_mat,bws_sandwich.median_mat,...
    ges_cafe.median_mat,ges_sandwich.median_mat,...
    hmf_cafe.median_mat),3);
all_mean = nanmean(cat(3,...
    bws_inside_walk_1.median_mat,bws_inside_walk_2.median_mat,bws_inside_walk_3.median_mat,...
    ges_inside_walk_1.median_mat,ges_inside_walk_2.median_mat,ges_inside_walk_3.median_mat,...
    hmf_inside_walk_1.median_mat,hmf_inside_walk_2.median_mat,hmf_inside_walk_3.median_mat,...
    bws_nature_walk_1.median_mat,bws_nature_walk_2.median_mat,bws_campus_walk.median_mat,...
    ges_nature_walk_1.median_mat,ges_nature_walk_2.median_mat,ges_campus_walk.median_mat,...
    hmf_nature_walk_1.median_mat,hmf_nature_walk_2.median_mat,hmf_campus_walk.median_mat,...
    bws_cafe.median_mat,bws_sandwich.median_mat,...
    ges_cafe.median_mat,ges_sandwich.median_mat,...
    hmf_cafe.median_mat),3);

all_inside_walks_mean(all_inside_walks_count < min_count_walks_all) = NaN;
all_outside_walks_mean(all_outside_walks_count < min_count_walks_all) = NaN;
all_other_mean(all_other_tasks_count < min_count_other_all) = NaN;




eccentricity_deg = 0.0983*(-103:103);
deg_8 = find(eccentricity_deg <=8 & eccentricity_deg >=-8);

degperpixel = 2*atand(.5/583);
rad = ceil(8./degperpixel);

%crop to 8 deg radius
all_outside_walks_mean = all_outside_walks_mean(104-82:104+82,104-82:104+82,:);

y = 1:size(all_outside_walks_mean,1);
x = 1:size(all_outside_walks_mean,2);
[x y] = meshgrid(x,y);

[row col] = find((y-83).^2 + (x-83).^2 > rad.^2);

all_outside_walks_mean(row,col,:) = NaN;



fig1 = figure(); hold on;
max_abs_disparity = max(max(abs(all_outside_walks_mean)));
sc(flipud(all_outside_walks_mean),'diff',[-max_abs_disparity max_abs_disparity],[1 1 1],isnan(all_outside_walks_mean));
cbar = colorbar;
ylabel(cbar,'median disparity (deg)'); title('outside walks');

fig1 = figure(); hold on;
max_abs_disparity = max(max(abs(all_inside_walks_mean)));
sc(flipud(all_inside_walks_mean),'diff',[-max_abs_disparity max_abs_disparity],[1 1 1],isnan(all_inside_walks_mean));
cbar = colorbar;
ylabel(cbar,'median disparity (deg)'); title('inside walks');

fig1 = figure(); hold on;
max_abs_disparity = max(max(abs(all_other_mean)));
sc(flipud(all_other_mean),'diff',[-max_abs_disparity max_abs_disparity],[1 1 1],isnan(all_other_mean));
cbar = colorbar;
ylabel(cbar,'median disparity (deg)'); title('other tasks');

keyboard