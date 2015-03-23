function [] = plotting(subj)

norm = 1;

if strcmp(subj,'bws')
    nw1 = load('../data/bwsnat5/bwsnat5_task_nature_walk_1_disparity_downsampled.mat');
    nw2 = load('../data/bwsure1/bwsure1_task_nature_walk_2_disparity_downsampled.mat');
    iw = load('../data/bwsdrdp6/bwsdrdp6_task_walk_1_disparity_downsampled.mat');
    ca = load('../data/bwscafe1/bwscafe1_task_ordering_coffee_disparity_downsampled.mat');
    sa = load('../data/bwssand1/bwssand1_task_sandwich_disparity_downsampled.mat');
end

nw1 = nw1.disparity_mat;
nw2 = nw2.disparity_mat;
nwall = cat(3,nw1,nw2);
iw = iw.disparity_mat;
ca = ca.disparity_mat;
sa = sa.disparity_mat;
all = cat(3,nw1,nw2,iw,ca,sa);

%calculate gradients in x,y,z 
%pixel spacing in degrees
xyspacing = 2*atand(1.5/583);
%frame spacing in ms
zspacing = 33.3;
%spacing in x and y is in degrees

[x_all,y_all,z_all] = gradient(all,xyspacing,xyspacing,zspacing);

keyboard
xgradient_mean_all = nanmean(abs(x_all),3);
ygradient_mean_all = nanmean(abs(y_all),3);
zgradient_mean_all = nanmean(abs(z_all),3);

x_thresh = sum(x_all >= 1,3);
y_thresh = sum(y_all >= 1,3);
z_thresh = sum(z_all >= 1,3);

xcount = sqrt(sum(~isnan(x_all),3));
ycount = sqrt(sum(~isnan(y_all),3));
zcount = sqrt(sum(~isnan(z_all),3));

x_percent = x_thresh./xcount;
y_percent = y_thresh./ycount;
z_percent = z_thresh./zcount;

max_percent = max(max(max([x_percent ; y_percent])));
min_percent = min(min(min([x_percent ; y_percent])));

xgradient_ci_all = 1.96.*(nanstd(x_all,0,3)./xcount);
ygradient_ci_all = 1.96.*(nanstd(y_all,0,3)./ycount);
zgradient_ci_all = 1.96.*(nanstd(z_all,0,3)./zcount);

max_abs_xgradient_all = max(max(max(xgradient_mean_all)));
max_abs_ygradient_all = max(max(max(ygradient_mean_all)));
max_abs_zgradient_all = max(max(max(zgradient_mean_all)));
min_abs_xgradient_all = min(min(min(xgradient_mean_all)));
min_abs_ygradient_all = min(min(min(ygradient_mean_all)));
min_abs_zgradient_all = min(min(min(zgradient_mean_all)));


disparityplot_1 = [ NaN*ones(69,5) flipud(x_percent) NaN*ones(69,5) flipud(y_percent) NaN*ones(69,5) flipud(y_percent) NaN*ones(69,5)];

disparityplot_1 = cat(1,NaN*ones(10,size(disparityplot_1,2)),disparityplot_1,NaN*ones(10,size(disparityplot_1,2)));

[disparityplot_1_mask] = make_circles(disparityplot_1);


fig1 = figure(); hold on;
max_abs_disparity = max([max_abs_xgradient_all max_abs_ygradient_all max_abs_zgradient_all]);
min_abs_disparity = min([min_abs_xgradient_all min_abs_ygradient_all min_abs_zgradient_all]);

sc(disparityplot_1,'gray',[min_percent max_percent],[1 1 1],isnan(disparityplot_1),[0 0 0],isnan(disparityplot_1_mask));

cbar = colorbar;
ylabel(cbar,'mean disparity gradient (deg)');

text('horizontalalignment','center','units','normalized','position',[.17,.935],'fontsize',14,'string','horizontal');
text('horizontalalignment','center','units','normalized','position',[.45,.935],'fontsize',14,'string','vertical');
text('horizontalalignment','center','units','normalized','position',[.73,.935],'fontsize',14,'string','temporal');

keyboard
% if norm == 1
%     print(fig1,'-depsc','-r300',['../data/plots/' subj 'mean_maps_norm.eps']);
% else
%     print(fig1,'-depsc','-r300',['../data/plots/' subj 'mean_maps.eps']);
% end

eccentricity_deg = -0.0983*2*((1:length(disparity_mean_all)) - (length(disparity_mean_all)/2));

vertical_slice_nw1 = disparity_mean_nw1(:,35);
vertical_slice_nw1_ci = disparity_ci_nw1(:,35);

horizontal_slice_nw1 = disparity_mean_nw1(35,:);
horizontal_slice_nw1_ci = disparity_ci_nw1(35,:);

vertical_slice_nw2 = disparity_mean_nw2(:,35);
vertical_slice_nw2_ci = disparity_ci_nw2(:,35);

horizontal_slice_nw2 = disparity_mean_nw2(35,:);
horizontal_slice_nw2_ci = disparity_ci_nw2(35,:);

vertical_slice_nwall = disparity_mean_nwall(:,35);
vertical_slice_nwall_ci = disparity_ci_nwall(:,35);

horizontal_slice_nwall = disparity_mean_nwall(35,:);
horizontal_slice_nwall_ci = disparity_ci_nwall(35,:);

vertical_slice_iw = disparity_mean_iw(:,35);
vertical_slice_iw_ci = disparity_ci_iw(:,35);

horizontal_slice_iw = disparity_mean_iw(35,:);
horizontal_slice_iw_ci = disparity_ci_iw(35,:);

vertical_slice_ca = disparity_mean_ca(:,35);
vertical_slice_ca_ci = disparity_ci_ca(:,35);

horizontal_slice_ca = disparity_mean_ca(35,:);
horizontal_slice_ca_ci = disparity_ci_ca(35,:);

vertical_slice_sa = disparity_mean_sa(:,35);
vertical_slice_sa_ci = disparity_ci_sa(:,35);

horizontal_slice_sa = disparity_mean_sa(35,:);
horizontal_slice_sa_ci = disparity_ci_sa(35,:);

vertical_slice_all = disparity_mean_all(:,35);
vertical_slice_all_ci = disparity_ci_all(:,35);

horizontal_slice_all = disparity_mean_all(35,:);
horizontal_slice_all_ci = disparity_ci_all(35,:);

%vertical_std = disparity_std(:,51:53);
%vertical_std = nanmean(vertical_std,2);

fig2 = figure(); hold on;
h(1) = plot(eccentricity_deg,vertical_slice_nw1,'r');
boundedline(eccentricity_deg,vertical_slice_nw1,vertical_slice_nw1_ci,'r','alpha');
h(2) = plot(eccentricity_deg,vertical_slice_nw2,'m');
boundedline(eccentricity_deg,vertical_slice_nw2,vertical_slice_nw2_ci,'m','alpha');
h(3) = plot(eccentricity_deg,vertical_slice_nwall,'b');
boundedline(eccentricity_deg,vertical_slice_nwall,vertical_slice_nwall_ci,'b','alpha');
h(4) = plot(eccentricity_deg,vertical_slice_iw,'c');
boundedline(eccentricity_deg,vertical_slice_iw,vertical_slice_iw_ci,'c','alpha');
h(5) = plot(eccentricity_deg,vertical_slice_ca,'g');
boundedline(eccentricity_deg,vertical_slice_ca,vertical_slice_ca_ci,'g','alpha');
h(6) = plot(eccentricity_deg,vertical_slice_sa,'y');
boundedline(eccentricity_deg,vertical_slice_sa,vertical_slice_sa_ci,'y','alpha');
h(7) = plot(eccentricity_deg,vertical_slice_all,'k');
boundedline(eccentricity_deg,vertical_slice_all,vertical_slice_all_ci,'k','alpha');

legend(h,'nature walk 1','nature walk 2','nature walk both','inside walk','cafe','sandwich','all');
xlabel('vertical eccentricity deg (down/up)');
ylabel('disparity deg (uncrossed/crossed)');

% if norm == 1
%     print(fig2,'-depsc','-r300',['../data/plots/' subj 'vertical_norm.eps']);
% else
%     print(fig2,'-depsc','-r300',['../data/plots/' subj 'vertical.eps']);
% end

fig3 = figure(); hold on;
h(1) = plot(eccentricity_deg,horizontal_slice_nw1,'r');
boundedline(eccentricity_deg,horizontal_slice_nw1,horizontal_slice_nw1_ci,'r','alpha');
h(2) = plot(eccentricity_deg,horizontal_slice_nw2,'m');
boundedline(eccentricity_deg,horizontal_slice_nw2,horizontal_slice_nw2_ci,'m','alpha');
h(3) = plot(eccentricity_deg,horizontal_slice_nwall,'b');
boundedline(eccentricity_deg,horizontal_slice_nwall,horizontal_slice_nwall_ci,'b','alpha');
h(4) = plot(eccentricity_deg,horizontal_slice_iw,'c');
boundedline(eccentricity_deg,horizontal_slice_iw,horizontal_slice_iw_ci,'c','alpha');
h(5) = plot(eccentricity_deg,horizontal_slice_ca,'g');
boundedline(eccentricity_deg,horizontal_slice_ca,horizontal_slice_ca_ci,'g','alpha');
h(6) = plot(eccentricity_deg,horizontal_slice_sa,'y');
boundedline(eccentricity_deg,horizontal_slice_sa,horizontal_slice_sa_ci,'y','alpha');
h(7) = plot(eccentricity_deg,horizontal_slice_all,'k');
boundedline(eccentricity_deg,horizontal_slice_all,horizontal_slice_all_ci,'k','alpha');

legend(h,'nature walk 1','nature walk 2','nature walk both','inside walk','cafe','sandwich','all');
xlabel('horziontal eccentricity deg (left/right)');
ylabel('disparity deg (uncrossed/crossed)');

% if norm == 1
%     print(fig3,'-depsc','-r300',['../data/plots/' subj 'horizontal_norm.eps']);
% else
%     print(fig3,'-depsc','-r300',['../data/plots/' subj 'horizontal.eps']);
% end

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