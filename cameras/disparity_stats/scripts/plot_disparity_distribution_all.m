
%bws_inside_walk_2 = load('../data/bwsdrdp6/bwsdrdp6_task_walk_2_disparity.mat');
%bws_inside_walk_3 = load('../data/bwsdrdp6/bwsdrdp6_task_walk_3_disparity.mat');

figure(); hold on; hold all;
bws_nature_walk_1 = load('../data/bwsnat5/bwsnat5_task_nature_walk_1_disparity.mat');
plot_disparity_distribution(bws_nature_walk_1.disparity,'r')
clear all;

bws_inside_walk_1 = load('../data/bwsdrdp6/bwsdrdp6_task_walk_1_disparity.mat');
plot_disparity_distribution(bws_inside_walk_1.disparity,'b')
clear all;

bws_sandwich = load('../data/bwssand1/bwssand1_task_sandwich_disparity.mat');
plot_disparity_distribution(bws_sandwich.disparity,'c')
clear all;

legend('bws nature walk 1','bws inside walk 1','bws sandwich')

figure(); hold on; hold all;
ges_nature_walk_1 = load('../data/gesnat1/gesnat1_task_nature_walk_1_disparity.mat');
plot_disparity_distribution(ges_nature_walk_1.disparity,'r')
clear all;

ges_inside_walk_1 = load('../data/gesdrdp1/gesdrdp1_task_walk_1_disparity.mat');
plot_disparity_distribution(ges_inside_walk_1.disparity,'b')
clear all;

ges_sandwich = load('../data/gessand1/gessand1_task_sandwich_disparity.mat');
plot_disparity_distribution(ges_sandwich.disparity,'c')
clear all;

legend('ges nature walk 1','ges inside walk 1','ges sandwich')

keyboard

figure(); hold on; hold all;
bws_nature_walk_1 = load('../data/bwsnat5/bwsnat5_task_nature_walk_1_disparity.mat');
plot_disparity_distribution(bws_nature_walk_1.disparity,'r')
clear all;

bws_nature_walk_2 = load('../data/bwsure1/bwsure1_task_nature_walk_2_disparity.mat');
plot_disparity_distribution(bws_nature_walk_2.disparity,'r:')
clear all;

bws_campus_walk = load('../data/bwsnat5/bwsnat5_task_campus_walk_3_disparity.mat');
plot_disparity_distribution(bws_campus_walk.disparity,'r--')
clear all;

ges_nature_walk_1 = load('../data/gesnat1/gesnat1_task_nature_walk_1_disparity.mat');
plot_disparity_distribution(ges_nature_walk_1.disparity,'b')
clear all;

ges_nature_walk_2 = load('../data/gesure1/gesure1_task_nature_walk_2_disparity.mat');
plot_disparity_distribution(ges_nature_walk_2.disparity,'b:')
clear all;

ges_campus_walk = load('../data/gesnat1/gesnat1_task_campus_walk_3_disparity.mat');
plot_disparity_distribution(ges_campus_walk.disparity,'b--')
clear all;

hmf_nature_walk_1 = load('../data/hmfnat1/hmfnat1_task_nature_walk_1_disparity.mat');
plot_disparity_distribution(hmf_nature_walk_1.disparity,'g')
clear all;

hmf_nature_walk_2 = load('../data/hmfure1/hmfure1_task_nature_walk_2_disparity.mat');
plot_disparity_distribution(hmf_nature_walk_2.disparity,'g:')
clear all;

hmf_campus_walk = load('../data/hmfnat1/hmfnat1_task_campus_walk_3_disparity.mat');
plot_disparity_distribution(hmf_campus_walk.disparity,'g--')
clear all;

legend('bws nature 1','bws nature 2','bws campus','ges nature 1','ges nature 2','ges campus','hmf nature 1','hmf nature 2','hmf campus');

figure(); hold on;
bws_inside_walk_1 = load('../data/bwsdrdp6/bwsdrdp6_task_walk_1_disparity.mat');
plot_disparity_distribution(bws_inside_walk_1.disparity,'r')
clear all;

bws_inside_walk_2 = load('../data/bwsdrdp6/bwsdrdp6_task_walk_2_disparity.mat');
plot_disparity_distribution(bws_inside_walk_2.disparity,'r:')
clear all;

bws_inside_walk_3 = load('../data/bwsdrdp6/bwsdrdp6_task_walk_3_disparity.mat');
plot_disparity_distribution(bws_inside_walk_3.disparity,'r--')
clear all;

ges_inside_walk_1 = load('../data/gesdrdp1/gesdrdp1_task_walk_1_disparity.mat');
plot_disparity_distribution(ges_inside_walk_1.disparity,'b')
clear all;

ges_inside_walk_2 = load('../data/gesdrdp1/gesdrdp1_task_walk_2_disparity.mat');
plot_disparity_distribution(ges_inside_walk_2.disparity,'b:')
clear all;

ges_inside_walk_3 = load('../data/gesdrdp1/gesdrdp1_task_walk_3_disparity.mat');
plot_disparity_distribution(ges_inside_walk_3.disparity,'b--')
clear all;

hmf_inside_walk_1 = load('../data/hmfdrdp1/hmfdrdp1_task_walk_1_disparity.mat');
plot_disparity_distribution(hmf_inside_walk_1.disparity,'g')
clear all;

hmf_inside_walk_2 = load('../data/hmfdrdp1/hmfdrdp1_task_walk_2_disparity.mat');
plot_disparity_distribution(hmf_inside_walk_2.disparity,'g:')
clear all;

hmf_inside_walk_3 = load('../data/hmfdrdp1/hmfdrdp1_task_walk_3_disparity.mat');
plot_disparity_distribution(hmf_inside_walk_3.disparity,'g--')
clear all;

legend('bws walk 1','bws reverse','bws walk 2','ges walk 1','ges reverse','ges walk 2','hmf walk 1','hmf reverse','hmf walk 2');

figure(); hold on;
bws_cafe = load('../data/bwscafe1/bwscafe1_task_ordering_coffee_disparity.mat');
plot_disparity_distribution(bws_cafe.disparity,'r')
clear all;

bws_sandwich = load('../data/bwssand1/bwssand1_task_sandwich_disparity.mat');
plot_disparity_distribution(bws_sandwich.disparity,'r:')
clear all;

ges_cafe = load('../data/gescafe2/gescafe2_task_ordering_coffee_disparity.mat');
plot_disparity_distribution(ges_cafe.disparity,'b')
clear all;

ges_sandwich = load('../data/gessand1/gessand1_task_sandwich_disparity.mat');
plot_disparity_distribution(ges_sandwich.disparity,'b:')
clear all;

legend('bws cafe','bws sandwich','ges cafe','ges sandwich');

