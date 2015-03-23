function [] = mean_maps(subj)

%normalize data so disparity at fixation is zero?
norm = 0;

%mask out data below a percentage-of-valid-pixels threshhold
thresh =-1;

file_flag = '';
if norm ==1
    file_flag = [file_flag '_norm'];
end
if thresh > -1
    file_flag = [file_flag '_thresh'];
end

%load in average empirical horopter disparities
load('average_horizontal_horopter.mat');
load('average_vertical_horopter.mat');

%load in data
if strcmp(subj,'bws')
    nw1 = load('../data/bwsnat5/bwsnat5_task_nature_walk_1_disparity_downsampled.mat');
    nw2 = load('../data/bwsure1/bwsure1_task_nature_walk_2_disparity_downsampled.mat');
    iw = load('../data/bwsdrdp6/bwsdrdp6_task_walk_1_disparity_downsampled.mat');
    ca = load('../data/bwscafe1/bwscafe1_task_ordering_coffee_disparity_downsampled.mat');
    sa = load('../data/bwssand1/bwssand1_task_sandwich_disparity_downsampled.mat');
elseif strcmp(subj,'ges')
    nw1 = load('../data/gesnat1/gesnat1_task_nature_walk_1_disparity_downsampled.mat');
    nw2 = load('../data/gesure1/gesure1_task_nature_walk_2_disparity_downsampled.mat');
    iw = load('../data/gesdrdp1/gesdrdp1_task_walk_1_disparity_downsampled.mat');
    ca = load('../data/gescafe2/gescafe2_task_ordering_coffee_disparity_downsampled.mat');
    sa = load('../data/gessand1/gessand1_task_sandwich_disparity_downsampled.mat');
end

disparity.task{1}       = 'nature walk 1';
disparity.alldata{1}    = nw1.disparity_mat;

disparity.task{2}       = 'nature walk 2';
disparity.alldata{2}    = nw2.disparity_mat;

disparity.task{3}       = 'nature walk both';
disparity.alldata{3}    = cat(3,nw1.disparity_mat,nw2.disparity_mat);

disparity.task{4}       = 'inside walk';
disparity.alldata{4}    = iw.disparity_mat;

disparity.task{5}       = 'cafe';
disparity.alldata{5}    = ca.disparity_mat;

disparity.task{6}       = 'making sandwich';
disparity.alldata{6}    = sa.disparity_mat;

disparity.task{7}       = 'all';
disparity.alldata{7}    = cat(3,nw1.disparity_mat,nw2.disparity_mat,iw.disparity_mat,ca.disparity_mat,sa.disparity_mat);

%calculate mean disparity, total non NaN values, onesided 95% ci, and max abs mean for each pixel
for x = 1:length(disparity.task)
    disparity.mean{x}   = nanmean(disparity.alldata{x},3);
    disparity.count{x}  = sum(~isnan(disparity.alldata{x}),3);
    disparity.maxcount{x}  = max(max(disparity.count{x}));
    disparity.mincount{x}  = min(min(disparity.count{x}));
    disparity.percent{x}  = sum(~isnan(disparity.alldata{x}),3)./size(disparity.alldata{x},3);
    disparity.maxpercent{x}  = max(max(disparity.percent{x}));
    disparity.minpercent{x}  = min(min(disparity.percent{x}));
    disparity.ci{x}     = 1.96.*(nanstd(disparity.alldata{x},0,3)./sqrt(disparity.count{x}));
    disparity.maxci{x}  = max(max(disparity.ci{x}));
    disparity.minci{x}  = min(min(disparity.ci{x}));
    disparity.maxabsmean{x} = max(abs([max(max(max(disparity.mean{x}))) min(min(min(disparity.mean{x})))]));
    
    %mask out pixels below a certain threshold of measurements
    if thresh > -1
        disparity.mean{x}(disparity.percent{x} <= thresh) = NaN;
        disparity.count{x}(disparity.percent{x} <= thresh) = NaN;
        disparity.maxcount{x}  = max(max(disparity.count{x}));
        disparity.mincount{x}  = min(min(disparity.count{x}));
        disparity.percent{x}(disparity.percent{x} <= thresh) = NaN;
        disparity.maxpercent{x}  = max(max(disparity.percent{x}));
        disparity.minpercent{x}  = min(min(disparity.percent{x}));
        disparity.ci{x}(disparity.percent{x} <= thresh) = NaN;
        disparity.maxci{x}  = max(max(disparity.ci{x}));
        disparity.minci{x}  = min(min(disparity.ci{x}));
        disparity.maxabsmean{x} = max(abs([max(max(max(disparity.mean{x}))) min(min(min(disparity.mean{x})))]));
    end
end

%normalize disparity so that disparity at fixation is zero
if norm == 1
    
    for x = 1:length(disparity.task)
        disparity.mean{x} = disparity.mean{x} - disparity.mean{x}(ceil(length(disparity.mean{x})/2),ceil(length(disparity.mean{x})/2));
        disparity.maxabsmean{x} = max(abs([max(max(max(disparity.mean{x}))) min(min(min(disparity.mean{x})))]));
    end
    
end    
    
%put all the mean plots together in a nice image
disparityplot_row1 = [ NaN*ones(69,5) flipud(disparity.mean{1}) ...
                       NaN*ones(69,5) flipud(disparity.mean{2}) ...
                       NaN*ones(69,5) flipud(disparity.mean{3}) ...
                       NaN*ones(69,5)];
disparityplot_row2 = [ NaN*ones(69,5) flipud(disparity.mean{4}) ...
                       NaN*ones(69,5) flipud(disparity.mean{5}) ...
                       NaN*ones(69,5) flipud(disparity.mean{6}) ...
                       NaN*ones(69,5)];
disparityplot_row3 = [ NaN*ones(69,5) NaN*ones(69,69) ...
                       NaN*ones(69,5) flipud(disparity.mean{7}) ...
                       NaN*ones(69,5) NaN*ones(69,69) ...
                       NaN*ones(69,5)];

%pad top and bottom
disparityplot_row1 = [ NaN*ones(10,size(disparityplot_row1,2)) ; disparityplot_row1 ; NaN*ones(10,size(disparityplot_row1,2))];
disparityplot_row2 = [ NaN*ones(10,size(disparityplot_row2,2)) ; disparityplot_row2 ; NaN*ones(10,size(disparityplot_row2,2))];
disparityplot_row3 = [ NaN*ones(10,size(disparityplot_row3,2)) ; disparityplot_row3 ; NaN*ones(10,size(disparityplot_row3,2))];
disparityplot = cat(1,NaN*ones(10,size(disparityplot_row2,2)),disparityplot_row3,disparityplot_row2,disparityplot_row1,NaN*ones(10,size(disparityplot_row2,2)));

[disparityplot_mask] = make_circles(disparityplot);

fig1 = figure(); hold on;
max_abs_disparity = max([disparity.maxabsmean{:}]);
sc(disparityplot,'diff',[-max_abs_disparity max_abs_disparity],[1 1 1],isnan(disparityplot),[0 0 0],isnan(disparityplot_mask));
cbar = colorbar;
ylabel(cbar,'mean disparity (deg)');

text_place = [ .17 .935 ; .45 .935 ; .73 .935 ; .17 .625 ; .45 .625 ; .73 .625 ; .45 .325 ];

for x = 1:length(disparity.task)
    text('horizontalalignment','center','units','normalized','position',[text_place(x,1),text_place(x,2)],'fontsize',14,'string',disparity.task{x});
end

print(fig1,'-depsc','-r300',['../data/plots/' subj '_mean_maps' file_flag '.eps']);


%plot mean disparity as a function of vertical and horziontal eccentrivity

%angular subtense of a pixel * pixel index gives eccentricity in degrees
eccentricity_deg = -0.2948*((1:length(disparity.mean{1})) - round(length(disparity.mean{1})/2));

% colors
disparity.color = {'r','m','b','c','g','y','k'};

fig2 = figure(); hold on;
plot(vh_ecc,vh_disparity,'k--','LineWidth',1.5);
for x = 1:length(disparity.task)
    
    disparity.verticalslice_mean{x} = disparity.mean{x}(:,35);
    disparity.verticalslice_ci{x} = disparity.ci{x}(:,35);
    
    h(x) = plot(eccentricity_deg,disparity.verticalslice_mean{x},disparity.color{x});
    boundedline(eccentricity_deg(~isnan(disparity.verticalslice_mean{x})),disparity.verticalslice_mean{x}(~isnan(disparity.verticalslice_mean{x})),disparity.verticalslice_ci{x}(~isnan(disparity.verticalslice_mean{x})),disparity.color{x},'alpha');
    
    if x == 7
        plot(eccentricity_deg,disparity.verticalslice_mean{x},disparity.color{x},'LineWidth',1.5);
    end
    
end

legend(h,disparity.task,'Location','EastOutside');
xlabel('vertical eccentricity deg (down/up)');
ylabel('disparity deg (uncrossed/crossed)');


print(fig2,'-depsc','-r300',['../data/plots/' subj '_vertical' file_flag '.eps']);


fig3 = figure(); hold on;
plot(hh_ecc,hh_disparity,'k--','LineWidth',1.5);
    
for x = 1:length(disparity.task)
    
    disparity.horizontalslice_mean{x} = disparity.mean{x}(35,:);
    disparity.horizontalslice_ci{x} = disparity.ci{x}(35,:);
    
    h(x) = plot(eccentricity_deg,disparity.horizontalslice_mean{x},disparity.color{x});
    boundedline(eccentricity_deg(~isnan(disparity.horizontalslice_mean{x})),disparity.horizontalslice_mean{x}(~isnan(disparity.horizontalslice_mean{x})),disparity.horizontalslice_ci{x}(~isnan(disparity.horizontalslice_mean{x})),disparity.color{x},'alpha');
    
    if x == 7
        plot(eccentricity_deg,disparity.horizontalslice_mean{x},disparity.color{x},'LineWidth',1.5);
    end
end
    
legend(h,disparity.task,'Location','EastOutside');
xlabel('horizontal eccentricity deg (left/right)');
ylabel('disparity deg (uncrossed/crossed)');


print(fig3,'-depsc','-r300',['../data/plots/' subj '_horizontal' file_flag '.eps']);



%put all the count plots together in a nice image
disparityplot_row1 = [ NaN*ones(69,5) flipud(disparity.count{1}) ...
                       NaN*ones(69,5) flipud(disparity.count{2}) ...
                       NaN*ones(69,5) flipud(disparity.count{3}) ...
                       NaN*ones(69,5)];
disparityplot_row2 = [ NaN*ones(69,5) flipud(disparity.count{4}) ...
                       NaN*ones(69,5) flipud(disparity.count{5}) ...
                       NaN*ones(69,5) flipud(disparity.count{6}) ...
                       NaN*ones(69,5)];
disparityplot_row3 = [ NaN*ones(69,5) NaN*ones(69,69) ...
                       NaN*ones(69,5) flipud(disparity.count{7}) ...
                       NaN*ones(69,5) NaN*ones(69,69) ...
                       NaN*ones(69,5)];

%pad top and bottom
disparityplot_row1 = [ NaN*ones(10,size(disparityplot_row1,2)) ; disparityplot_row1 ; NaN*ones(10,size(disparityplot_row1,2))];
disparityplot_row2 = [ NaN*ones(10,size(disparityplot_row2,2)) ; disparityplot_row2 ; NaN*ones(10,size(disparityplot_row2,2))];
disparityplot_row3 = [ NaN*ones(10,size(disparityplot_row3,2)) ; disparityplot_row3 ; NaN*ones(10,size(disparityplot_row3,2))];
disparityplot = cat(1,NaN*ones(10,size(disparityplot_row2,2)),disparityplot_row3,disparityplot_row2,disparityplot_row1,NaN*ones(10,size(disparityplot_row2,2)));

[disparityplot_mask] = make_circles(disparityplot);

fig4 = figure(); hold on;
max_count = max([disparity.maxcount{:}]);
min_count = min([disparity.mincount{:}]);
sc(disparityplot,'jet',[min_count max_count],[1 1 1],isnan(disparityplot),[0 0 0],isnan(disparityplot_mask));
cbar = colorbar;
ylabel(cbar,'data points');

text_place = [ .17 .935 ; .45 .935 ; .73 .935 ; .17 .625 ; .45 .625 ; .73 .625 ; .45 .325 ];

for x = 1:length(disparity.task)
    text('horizontalalignment','center','units','normalized','position',[text_place(x,1),text_place(x,2)],'fontsize',14,'string',disparity.task{x});
end

print(fig4,'-depsc','-r300',['../data/plots/' subj '_count_maps' file_flag '.eps']);


%put all the ci plots together in a nice image
disparityplot_row1 = [ NaN*ones(69,5) flipud(2.*disparity.ci{1}) ...
                       NaN*ones(69,5) flipud(2.*disparity.ci{2}) ...
                       NaN*ones(69,5) flipud(2.*disparity.ci{3}) ...
                       NaN*ones(69,5)];
disparityplot_row2 = [ NaN*ones(69,5) flipud(2.*disparity.ci{4}) ...
                       NaN*ones(69,5) flipud(2.*disparity.ci{5}) ...
                       NaN*ones(69,5) flipud(2.*disparity.ci{6}) ...
                       NaN*ones(69,5)];
disparityplot_row3 = [ NaN*ones(69,5) NaN*ones(69,69) ...
                       NaN*ones(69,5) flipud(2.*disparity.ci{7}) ...
                       NaN*ones(69,5) NaN*ones(69,69) ...
                       NaN*ones(69,5)];

%pad top and bottom
disparityplot_row1 = [ NaN*ones(10,size(disparityplot_row1,2)) ; disparityplot_row1 ; NaN*ones(10,size(disparityplot_row1,2))];
disparityplot_row2 = [ NaN*ones(10,size(disparityplot_row2,2)) ; disparityplot_row2 ; NaN*ones(10,size(disparityplot_row2,2))];
disparityplot_row3 = [ NaN*ones(10,size(disparityplot_row3,2)) ; disparityplot_row3 ; NaN*ones(10,size(disparityplot_row3,2))];
disparityplot = cat(1,NaN*ones(10,size(disparityplot_row2,2)),disparityplot_row3,disparityplot_row2,disparityplot_row1,NaN*ones(10,size(disparityplot_row2,2)));

[disparityplot_mask] = make_circles(disparityplot);

fig5 = figure(); hold on;
max_ci = max(2.*[disparity.maxci{:}]);
min_ci = min(2.*[disparity.minci{:}]);
sc(disparityplot,'jet',[min_ci max_ci],[1 1 1],isnan(disparityplot),[0 0 0],isnan(disparityplot_mask));
cbar = colorbar;
ylabel(cbar,'95% confidence interval deg');

text_place = [ .17 .935 ; .45 .935 ; .73 .935 ; .17 .625 ; .45 .625 ; .73 .625 ; .45 .325 ];

for x = 1:length(disparity.task)
    text('horizontalalignment','center','units','normalized','position',[text_place(x,1),text_place(x,2)],'fontsize',14,'string',disparity.task{x});
end

print(fig5,'-depsc','-r300',['../data/plots/' subj '_ci_maps' file_flag '.eps']);

%put all the ci plots together in a nice image
disparityplot_row1 = [ NaN*ones(69,5) flipud(disparity.percent{1}) ...
                       NaN*ones(69,5) flipud(disparity.percent{2}) ...
                       NaN*ones(69,5) flipud(disparity.percent{3}) ...
                       NaN*ones(69,5)];
disparityplot_row2 = [ NaN*ones(69,5) flipud(disparity.percent{4}) ...
                       NaN*ones(69,5) flipud(disparity.percent{5}) ...
                       NaN*ones(69,5) flipud(disparity.percent{6}) ...
                       NaN*ones(69,5)];
disparityplot_row3 = [ NaN*ones(69,5) NaN*ones(69,69) ...
                       NaN*ones(69,5) flipud(disparity.percent{7}) ...
                       NaN*ones(69,5) NaN*ones(69,69) ...
                       NaN*ones(69,5)];

%pad top and bottom
disparityplot_row1 = [ NaN*ones(10,size(disparityplot_row1,2)) ; disparityplot_row1 ; NaN*ones(10,size(disparityplot_row1,2))];
disparityplot_row2 = [ NaN*ones(10,size(disparityplot_row2,2)) ; disparityplot_row2 ; NaN*ones(10,size(disparityplot_row2,2))];
disparityplot_row3 = [ NaN*ones(10,size(disparityplot_row3,2)) ; disparityplot_row3 ; NaN*ones(10,size(disparityplot_row3,2))];
disparityplot = cat(1,NaN*ones(10,size(disparityplot_row2,2)),disparityplot_row3,disparityplot_row2,disparityplot_row1,NaN*ones(10,size(disparityplot_row2,2)));

[disparityplot_mask] = make_circles(disparityplot);

fig6 = figure(); hold on;
max_percent = max([disparity.maxpercent{:}]);
min_percent = min([disparity.minpercent{:}]);
sc(disparityplot,'jet',[min_percent max_percent],[1 1 1],isnan(disparityplot),[0 0 0],isnan(disparityplot_mask));
cbar = colorbar;
ylabel(cbar,'percent valid pixels');

text_place = [ .17 .935 ; .45 .935 ; .73 .935 ; .17 .625 ; .45 .625 ; .73 .625 ; .45 .325 ];

for x = 1:length(disparity.task)
    text('horizontalalignment','center','units','normalized','position',[text_place(x,1),text_place(x,2)],'fontsize',14,'string',disparity.task{x});
end

print(fig6,'-depsc','-r300',['../data/plots/' subj '_percent_maps' file_flag '.eps']);

keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create color map
function cmap = create_color_map()

cmap = [];
for x = 0:0.05:1
    cmap(end+1,:) = [x x 1];
end
for x = 0.95:-0.05:0.05
    cmap(end+1,:) = [1 x x];
end
cmap(end+1,:) = [1 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%replace nans with mean disparity of 9 surrounding pixels, so long as at
%least 4 surrounding pixels are not nans
function im2 = replace_nans(im)

N = size(im,1);
R = (floor(N/2)-1)^2;
    
im2 = im;
[Y,X] = find( isnan(im) );
for k = 1 : length(X)
    if( ((X(k)-N/2).^2 + (Y(k)-N/2).^2) <= R )
        blk = im(Y(k)-1:Y(k)+1,X(k)-1:X(k)+1);
        ind = find( ~isnan(blk) );
        if( length(ind) >= 3 )
            im2(Y(k),X(k)) = sum( blk(ind) )/length(ind);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [disparityplot_mask] = make_circles(disparityplot);
disparityplot_mask = zeros(size(disparityplot,1),size(disparityplot,2));
for j = 2:size(disparityplot,1)-1
    for k = 2:size(disparityplot,2)-1
        if isnan(disparityplot(j,k))
            if sum(sum(~isnan(disparityplot(j-1:j+1,k-1:k+1)))) >= 1
                disparityplot_mask(j,k) = NaN;
            end
        end
    end
end