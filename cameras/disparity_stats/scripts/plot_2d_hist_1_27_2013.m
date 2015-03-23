clear all; close all;

subjects = {'bwsdrdp6', 'bwsnat5', 'bwsure1', 'bwsnat5', 'bwssand2', ...
    'gesdrdp1', 'gesnat1', 'gesure1', 'gesnat1', 'gessand1',...
    'tkidrdp1', 'tkinat1', 'tkiure2', 'tkinat1', 'tkisand2',...
    'all','all','all','all','all'};

activities = repmat({'walk_1',...
    'nature_walk_1',...
    'nature_walk_2',...
    'campus_walk_3',...
    'sandwich'},1,3);

flags = {'', '_random', '_shuffle'};

%subj = 'bwsdrdp6';
%flag = '_random';
%task = 'walk_1';

for s = 2
    for f = 1
        
        subj = subjects{s};
        flag = flags{f};
        task = activities{s};
        
        %min number of measurements per bin to plot in hist
        bin_min = 50;
        
        %load in horopter disparities
        %load('mean_vert_disp.mat')
        load('average_horizontal_horopter.mat')
        load('average_vertical_horopter.mat')
        
        %make vector of pixel eccentricity in helmholtz
        eccentricity_deg = [];
        for k = 1:207
            eccentricity_deg(k) = sign(k-104)*acosd(dot([583 0]./norm([583 0]),[583 k-104]./norm([583 k-104])));
        end
        deg_8 = find(eccentricity_deg <=8 & eccentricity_deg >=-8);
        eccentricity_deg = eccentricity_deg(deg_8);
        
        %load data
        dmat = load(['../data/' subj flag '/' subj '_task_' task '_disparity.mat']);
        dmat.disparity_vert = dmat.disparity(deg_8,104,:);
        dmat.disparity_horz = dmat.disparity(104,deg_8,:);
        dmat.disparity_horz = reshape(dmat.disparity_horz,size(dmat.disparity_vert,1),1,size(dmat.disparity_vert,3));
        
        clear dmat.disparity
        
        %normalize to fixation disparity
        dmat.disparity_vert_norm = zeros(163,size(dmat.disparity_vert,3));
        for j = 1:size(dmat.disparity_vert,3)
            dmat.disparity_vert_norm(:,j) = dmat.disparity_vert(:,1,j) - dmat.disparity_vert(82,:,j);
        end
        dmat.disparity_horz_norm = zeros(163,size(dmat.disparity_horz,3));
        for j = 1:size(dmat.disparity_horz,3)
            dmat.disparity_horz_norm(:,j) = dmat.disparity_horz(:,1,j) - dmat.disparity_horz(82,:,j);
        end
        
        
        %medians
        dmat.disparity_vert_median = nanmedian(dmat.disparity_vert,3);
        dmat.disparity_vert_norm_median = nanmedian(dmat.disparity_vert_norm,2);
        dmat.disparity_horz_median = nanmedian(dmat.disparity_horz,3);
        dmat.disparity_horz_norm_median = nanmedian(dmat.disparity_horz_norm,2);
        
        %eccentricity indices
        vd_ecc_all = repmat(eccentricity_deg',size(dmat.disparity_vert,3),1);
        
        dmat.disparity_vert_norm = reshape(dmat.disparity_vert_norm,1,163*size(dmat.disparity_vert,3));
        dmat.disparity_vert = reshape(dmat.disparity_vert,1,163*size(dmat.disparity_vert,3));
        dmat.disparity_horz_norm = reshape(dmat.disparity_horz_norm,1,163*size(dmat.disparity_horz,3));
        dmat.disparity_horz = reshape(dmat.disparity_horz,1,163*size(dmat.disparity_horz,3));
        
        %VERTICAL MERIDIAN
        %normed hist
        step1 = 0.2;
        step2 = 0.03;
        centers{1} = [linspace(-8,0,8/step1)];
        centers{1} = [centers{1} fliplr(-centers{1}(1:end-1))];
        centers{2} = [linspace(-15,0,15/step2)];
        centers{2} = [centers{2} fliplr(-centers{2}(1:end-1))];
        
        %find bounds of disparity bins > bin_min
        [N] = hist3([-vd_ecc_all , dmat.disparity_vert_norm'],centers);
        [row,col] = find(N > bin_min);
        %stop hist3 form summing up all measurements beyond the bins
        dmat.disparity_vert_norm(dmat.disparity_vert_norm < centers{2}(min(col))) = NaN;
        dmat.disparity_vert_norm(dmat.disparity_vert_norm > centers{2}(max(col))) = NaN;
        centers{2} = centers{2}(min(col)-2:max(col)+2);
        
        fig1 = figure(); hold on; title([subj ' ' flag ' ' task ' vertical normed']);
        hist3([-vd_ecc_all , dmat.disparity_vert_norm'],centers);
        set(get(gca,'child'),'FaceColor','interp','CDataMode','auto','Edgecolor','none');
        zData = get(get(gca,'child'),'ZData');
        
        %colorbar
        zData(zData < bin_min) = 1;
        d = log10(zData);
        d(d == 0) = NaN;
        mn = min(d(:));
        rng = max(d(:))-mn;
        d = 1+63*(d-mn)/rng;
        hC = colorbar;
        L = round(logspace(log10(round(min(zData(zData>1))*.1)*10),log10(round(max(zData(:))*.1)*10),6)*.1)*10;
        l = round(1+63*(log10(L)-mn)/rng);
        set(get(gca,'child'),'ZData',d)
        
        figure(1); hold on;
        set(hC,'YTick',l);
        set(hC,'YTicklabel',L);
        set(get(hC,'YLabel'),'string','count');
        
        plot3(vh_ecc,-vh_disparity,100*ones(length(vh_ecc)),'k','LineWidth',2)
        plot3(-eccentricity_deg',dmat.disparity_vert_norm_median,100*ones(length(eccentricity_deg)),'r','LineWidth',2)
        plot3([-10 10],[0 0],[100 100],'k:');
        plot3([0 0],[.2 -.2],[100 100],'k:');
        xlabel('eccentricity(deg)'); ylabel('disparity(deg)');
        print(fig1,'-depsc','-r150',['plots_1_27_2013/' subj '_' task flag '_vertical_normed_hist.eps']);
        
        if strcmp(flag,'_random') == 0
            %unnormed hist
            step1 = 0.2;
            step2 = 0.03;
            centers{1} = [linspace(-8,0,8/step1)];
            centers{1} = [centers{1} fliplr(-centers{1}(1:end-1))];
            centers{2} = [linspace(-15,0,15/step2)];
            centers{2} = [centers{2} fliplr(-centers{2}(1:end-1))];
            
            %find bounds of disparity bins > bin_min
            [N] = hist3([-vd_ecc_all , dmat.disparity_vert'],centers);
            [row,col] = find(N > bin_min/2);
            %stop hist3 form summing up all measurements beyond the bins
            dmat.disparity_vert(dmat.disparity_vert < centers{2}(min(col))) = NaN;
            dmat.disparity_vert(dmat.disparity_vert > centers{2}(max(col))) = NaN;
            centers{2} = centers{2}(min(col):max(col));
            
            
            fig2 = figure(); hold on; title([subj ' ' flag ' ' task ' vertical unnormed']);
            hist3([-vd_ecc_all , dmat.disparity_vert'],centers);
            set(get(gca,'child'),'FaceColor','interp','CDataMode','auto','Edgecolor','none');
            zData = get(get(gca,'child'),'ZData');
            
            %colorbar
            zData(zData < bin_min/2) = 1;
            d = log10(zData);
            d(d == 0) = NaN;
            mn = min(d(:));
            rng = max(d(:))-mn;
            d = 1+63*(d-mn)/rng;
            hC = colorbar;
            L = round(logspace(log10(round(min(zData(zData>1)))),log10(round(max(zData(:)))),6));
            %L = round(logspace(log10(round(min(zData(zData>1))*.1)*10),log10(round(max(zData(:))*.1)*10),6)*.1)*10;
            l = round(1+63*(log10(L)-mn)/rng);
            set(get(gca,'child'),'ZData',d)
            
            figure(2); hold on;
            set(hC,'YTick',l);
            set(hC,'YTicklabel',L);
            set(get(hC,'YLabel'),'string','count');
            
            plot3(-eccentricity_deg',dmat.disparity_vert_median,100*ones(length(eccentricity_deg)),'r','LineWidth',3)
            plot3([-10 10],[0 0],[100 100],'k:');
            plot3([0 0],[.2 -.2],[100 100],'k:');
            
            xlabel('eccentricity(deg)'); ylabel('disparity(deg)');
            print(fig2,'-depsc','-r150',['plots_1_27_2013/' subj '_' task flag '_vertical_unnormed_hist.eps']);
            
        end
        
        
        %HORIZONTAL MERIDIAN
        %normed hist
        step1 = 0.2;
        step2 = 0.03;
        centers{1} = [linspace(-8,0,8/step1)];
        centers{1} = [centers{1} fliplr(-centers{1}(1:end-1))];
        centers{2} = [linspace(-15,0,15/step2)];
        centers{2} = [centers{2} fliplr(-centers{2}(1:end-1))];
        
        %find bounds of disparity bins > bin_min
        [N] = hist3([vd_ecc_all , dmat.disparity_horz_norm'],centers);
        [row,col] = find(N > bin_min);
        %stop hist3 form summing up all measurements beyond the bins
        dmat.disparity_horz_norm(dmat.disparity_horz_norm < centers{2}(min(col))) = NaN;
        dmat.disparity_horz_norm(dmat.disparity_horz_norm > centers{2}(max(col))) = NaN;
        centers{2} = centers{2}(min(col)-2:max(col)+2);
        
        fig3 = figure(); hold on; title([subj ' ' flag ' ' task ' horizontal normed']);
        hist3([vd_ecc_all , dmat.disparity_horz_norm'],centers);
        set(get(gca,'child'),'FaceColor','interp','CDataMode','auto','Edgecolor','none');
        zData = get(get(gca,'child'),'ZData');
        
        %colorbar
        zData(zData < bin_min) = 1;
        d = log10(zData);
        d(d == 0) = NaN;
        mn = min(d(:));
        rng = max(d(:))-mn;
        d = 1+63*(d-mn)/rng;
        hC = colorbar;
        L = round(logspace(log10(round(min(zData(zData>1))*.1)*10),log10(round(max(zData(:))*.1)*10),6)*.1)*10;
        l = round(1+63*(log10(L)-mn)/rng);
        set(get(gca,'child'),'ZData',d)
        
        if strcmp(flag,'_random') == 0
            figure(3); hold on;
        else
            figure(2); hold on;
        end
        set(hC,'YTick',l);
        set(hC,'YTicklabel',L);
        set(get(hC,'YLabel'),'string','count');
        
        plot3(hh_ecc(3:19),-hh_disparity(3:19),100*ones(length(hh_ecc(3:19))),'k','LineWidth',2)
        plot3(eccentricity_deg',dmat.disparity_horz_norm_median,100*ones(length(eccentricity_deg)),'r','LineWidth',2)
        plot3([-10 10],[0 0],[100 100],'k:');
        plot3([0 0],[.2 -.2],[100 100],'k:');
        
        xlabel('eccentricity(deg)'); ylabel('disparity(deg)');
        print(fig3,'-depsc','-r150',['plots_1_27_2013/' subj '_' task flag '_horizontal_normed_hist.eps']);
        
        %plot(vh_ecc,-vh_disparity,'k')
        %plot(-eccentricity_deg',dmat.disparity_horz_norm_median,'r-')
        
        if strcmp(flag,'_random') == 0
            
            %unnormed hist
            step1 = 0.2;
            step2 = 0.03;
            centers{1} = [linspace(-8,0,8/step1)];
            centers{1} = [centers{1} fliplr(-centers{1}(1:end-1))];
            centers{2} = [linspace(-15,0,15/step2)];
            centers{2} = [centers{2} fliplr(-centers{2}(1:end-1))];
            
            %find bounds of disparity bins > bin_min
            [N] = hist3([vd_ecc_all , dmat.disparity_horz'],centers);
            [row,col] = find(N > bin_min/2);
            %stop hist3 form summing up all measurements beyond the bins
            dmat.disparity_horz(dmat.disparity_horz < centers{2}(min(col))) = NaN;
            dmat.disparity_horz(dmat.disparity_horz > centers{2}(max(col))) = NaN;
            centers{2} = centers{2}(max([1 min(col)-2]):min([length(centers{2}) max(col)+2]));
            
            
            fig4 = figure(); hold on; title([subj ' ' flag ' ' task ' horizontal unnormed']);
            hist3([vd_ecc_all , dmat.disparity_horz'],centers);
            set(get(gca,'child'),'FaceColor','interp','CDataMode','auto','Edgecolor','none');
            zData = get(get(gca,'child'),'ZData');
            
            %colorbar
            zData(zData < bin_min/2) = 1;
            d = log10(zData);
            d(d == 0) = NaN;
            mn = min(d(:));
            rng = max(d(:))-mn;
            d = 1+63*(d-mn)/rng;
            hC = colorbar;
            L = round(logspace(log10(round(min(zData(zData>1)))),log10(round(max(zData(:)))),6));
            
            %L = round(logspace(log10(round(min(zData(zData>1))*.1)*10),log10(round(max(zData(:))*.1)*10),6)*.1)*10;
            l = round(1+63*(log10(L)-mn)/rng);
            set(get(gca,'child'),'ZData',d)
            
            figure(4); hold on;
            set(hC,'YTick',l);
            set(hC,'YTicklabel',L);
            set(get(hC,'YLabel'),'string','count');
            
            plot3(eccentricity_deg',dmat.disparity_horz_median,100*ones(length(eccentricity_deg)),'r','LineWidth',3)
            plot3([-10 10],[0 0],[100 100],'k:');
            plot3([0 0],[.2 -.2],[100 100],'k:');
            
            xlabel('eccentricity(deg)'); ylabel('disparity(deg)');
            print(fig4,'-depsc','-r150',['plots_1_27_2013/' subj '_' task flag '_horizontal_unnormed_hist.eps']);
            
        end
        
        close all;
        
    end
end