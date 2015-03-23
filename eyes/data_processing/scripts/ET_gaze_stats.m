function [] = ET_gaze_stats(subj,taskstr)
%subj = 'bwsnat5'
%task = 'task_nature_walk_1'

dst_dir = ['../data/' subj '/plots/'];
mkdir(dst_dir);

%load subj and tasks
if strcmp('bws',subj) && isempty(taskstr)
    tasks = [{'bwsnat5','task_nature_walk_1'},{'bwsnat5','task_campus_walk_3'},...
        {'bwsure1','task_nature_walk_2'},{'bwsdrdp6','task_walk_1'},...
        {'bwscafe1','task_ordering_coffee'},{'bwssand1','task_sandwich'}];
elseif strcmp('ges',subj) && isempty(taskstr)
    tasks = [{'gesnat1','task_nature_walk_1'},{'gesnat1','task_campus_walk_3'},...
        {'gesure1','task_nature_walk_2'},{'gesdrdp1','task_walk_1'},...
        {'gescafe2','task_ordering_coffee'},{'gessand1','task_sandwich'}];
elseif strcmp('hmf',subj) && isempty(taskstr)
    tasks = [{'hmfnat1','task_nature_walk_1'},{'hmfnat1','task_campus_walk_3'},...
        {'hmfure1','task_nature_walk_2'},{'hmfdrdp1','task_walk_1'},...
        {'hmfcafe1','task_ordering_coffee'}];    
else
    tasks = [{subj,taskstr}];
end

mat = [];
subj_label = [];
task_label = [];
verg_label = []; verg_vel_label = [];
hvers_label = []; hvers_vel_label = [];
vvers_label = []; vvers_vel_label = [];

for t = 1:2:length(tasks)-1
    
    %load data
    subj1 = tasks{t};
    taskstr1 = tasks{t+1};
    load(['../data/' subj1 '/' subj1 '_eye_movements.mat']);
    vvv = [];
    vvv_label = []; vvv_meas_label = [];
    
    load(['/users/natdispstats/Documents/cameras/disparity_stats/data/' subj1 '/' subj1 '_' taskstr1 '_disparity_coverage.mat']);
    
    [x] = textread(['/users/natdispstats/Documents/cameras/disparity_stats/data/' subj1 '/' subj1 '_' taskstr1 '_disparity_params.txt'], ...
        '%*s %*s %s',3);
    
    start_frame = str2num(cell2mat(x(2)));
    end_frame = str2num(cell2mat(x(3)));
    
    if strfind(subj1,'drdp')
        mat_temp = walk_1(start_frame+1:end_frame+1,:);
    else
        mat_temp = task(start_frame+1:end_frame+1,:);
    end
    
    %add visual field coverage info to data matrix (6deg vf radius)
    cnt = 1;
    for k = 1:size(mat_temp,1)
        if mat_temp(k,14) == 1 || mat_temp(k,14) == 3
            mat_temp(k,15) = coverage_mat(cnt,3);
            cnt = cnt + 1;
        elseif mat_temp(k,14) == 4
            mat_temp(k,15) = 0;
        elseif mat_temp(k,14) == 2
            mat_temp(k,15) = NaN;    
        end
    end
    
    mat = [mat ; mat_temp];

    %calculate percentiles and outliers of v/v/v, make into labels
    verg_label = [verg_label ; repmat({[num2str(quantile(mat_temp(:,11),[0.25 0.50]),'%.1f/') num2str(quantile(mat_temp(:,11),[0.75]),'%.1f') ...
        ' otlrs:' num2str(sum(mat_temp(:,11) < quantile(mat_temp(:,11),0.25) - 1.5*(quantile(mat_temp(:,11),0.75) - quantile(mat_temp(:,11),0.25))) + sum(mat_temp(:,11) > quantile(mat_temp(:,11),0.75) + 1.5*(quantile(mat_temp(:,11),0.75) - quantile(mat_temp(:,11),0.25)))) ...
        '/' num2str(sum(~isnan(mat_temp(:,11))))]},...
        length(mat_temp(:,11)),1)];
    hvers_label = [hvers_label ; repmat({[num2str(quantile(mat_temp(:,12),[0.25 0.50]),'%.1f/') num2str(quantile(mat_temp(:,12),[0.75]),'%.1f') ...
        ' otlrs:' num2str(sum(mat_temp(:,12) < quantile(mat_temp(:,12),0.25) - 1.5*(quantile(mat_temp(:,12),0.75) - quantile(mat_temp(:,12),0.25))) + sum(mat_temp(:,12) > quantile(mat_temp(:,12),0.75) + 1.5*(quantile(mat_temp(:,12),0.75) - quantile(mat_temp(:,12),0.25)))) ...
        '/' num2str(sum(~isnan(mat_temp(:,12))))]},...
        length(mat_temp(:,12)),1)];
    vvers_label = [vvers_label ; repmat({[num2str(quantile(mat_temp(:,13),[0.25 0.50]),'%.1f/') num2str(quantile(mat_temp(:,13),[0.75]),'%.1f') ...
        ' otlrs:' num2str(sum(mat_temp(:,13) < quantile(mat_temp(:,13),0.25) - 1.5*(quantile(mat_temp(:,13),0.75) - quantile(mat_temp(:,13),0.25))) + sum(mat_temp(:,13) > quantile(mat_temp(:,13),0.75) + 1.5*(quantile(mat_temp(:,13),0.75) - quantile(mat_temp(:,13),0.25)))) ...
        '/' num2str(sum(~isnan(mat_temp(:,13))))]},...
        length(mat_temp(:,13)),1)];
    
    
    subj_label = [subj_label ; repmat({subj1},length(mat_temp(:,11)),1)];
    task_label = [task_label ; repmat({taskstr1},length(mat_temp(:,11)),1)];
    
    %concate for individual plots
    vvv = [mat_temp(:,11) ; mat_temp(:,12) ; mat_temp(:,13)];
    vvv_label = {[repmat({'vergence'},length(mat_temp(:,11)),1) ; repmat({'h version'},length(mat_temp(:,12)),1) ; repmat({'v version'},length(mat_temp(:,13)),1)]};
    
    vvv_meas_label = {[repmat({[num2str(quantile(mat_temp(:,11),[0.25 0.50]),'%.1f/') num2str(quantile(mat_temp(:,11),[0.75]),'%.1f') ...
        ' otlrs:' num2str(sum(mat_temp(:,11) < quantile(mat_temp(:,11),0.25) - 1.5*(quantile(mat_temp(:,11),0.75) - quantile(mat_temp(:,11),0.25))) + sum(mat_temp(:,11) > quantile(mat_temp(:,11),0.75) + 1.5*(quantile(mat_temp(:,11),0.75) - quantile(mat_temp(:,11),0.25)))) ...
        '/' num2str(sum(~isnan(mat_temp(:,11))))]},...
        length(mat_temp(:,11)),1) ; ...
        repmat({[num2str(quantile(mat_temp(:,12),[0.25 0.50]),'%.1f/') num2str(quantile(mat_temp(:,12),[0.75]),'%.1f') ...
        ' otlrs:' num2str(sum(mat_temp(:,12) < quantile(mat_temp(:,12),0.25) - 1.5*(quantile(mat_temp(:,12),0.75) - quantile(mat_temp(:,12),0.25))) + sum(mat_temp(:,12) > quantile(mat_temp(:,12),0.75) + 1.5*(quantile(mat_temp(:,12),0.75) - quantile(mat_temp(:,12),0.25)))) ...
        '/' num2str(sum(~isnan(mat_temp(:,12))))]},...
        length(mat_temp(:,12)),1) ; ...
        repmat({[num2str(quantile(mat_temp(:,13),[0.25 0.50]),'%.1f/') num2str(quantile(mat_temp(:,13),[0.75]),'%.1f') ...
        ' otlrs:' num2str(sum(mat_temp(:,13) < quantile(mat_temp(:,13),0.25) - 1.5*(quantile(mat_temp(:,13),0.75) - quantile(mat_temp(:,13),0.25))) + sum(mat_temp(:,13) > quantile(mat_temp(:,13),0.75) + 1.5*(quantile(mat_temp(:,13),0.75) - quantile(mat_temp(:,13),0.25)))) ...
        '/' num2str(sum(~isnan(mat_temp(:,13))))]},...
        length(mat_temp(:,13)),1) ]};
    
    %how much missing data?
    blinks = sum(mat_temp(:,14) == 2)/size(mat_temp,1);
    missing = sum(mat_temp(:,14) == 4)/size(mat_temp,1);
    

    
    
    %vergence, version box plots
    fig1 = figure(); hold on;  
    scrsz = get(0,'ScreenSize');
    set(fig1,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
    subplot(1,6,1:3); hold on; title(regexprep([subj1 ' : ' taskstr1],'_',' '));
    h = boxplot(vvv,{vvv_label vvv_meas_label});
    delete(h(end,:)); axis auto;
    plot([.5 3.5],[0 0],'k--');
    ylabel('deg');
    

    subplot(1,6,4:6); hold on;
    bins = 0:10:90;
    n = [];
    for j = 1:length(bins)
        n(j) = length(find(mat_temp(:,15) >= bins(j) & mat_temp(:,15) < bins(j)+10));
    end
    %[n,x] = hist(mat_temp(:,15));
    bins(end+1) = -20;
    n(end+1) = sum(mat_temp(:,14) == 2);
    bins(end+1) = -30;
    n(end+1) = sum(mat_temp(:,14) == 4);
    h = bar(bins,(n./size(mat_temp,1))*100);
    labels = {'missing','blink','0-9','10-19','20-29','30-39','40-49',...
        '50-59','60-69','70-79','80-89','90-99'};
    set(gca,'XTick',sort(bins));
    set(gca,'XTickLabel',labels);
    xlabel('visual field coverage distribution'); ylabel('% of data');
    %x_loc = get(h,'Xdata');
    %y_height = get(h,'Ydata');
    arrayfun(@(x,y) text(x,y+0.2,num2str(y),'Color','k'),bins,2+round((n./size(mat_temp,1))*100));
    
    fig2 = figure(); hold on;  
    scrsz = get(0,'ScreenSize');
    set(fig2,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
    subplot(1,3,1); hold on;
    scatter(mat_temp(:,12),mat_temp(:,13),5,mat_temp(:,15));
    plot([-90 90],[0 0],'k--');
    plot([-90 90],[quantile(mat_temp(:,13),0.50) quantile(mat_temp(:,13),0.50)],'r');
    plot([-90 90],[quantile(mat_temp(:,13),0.25) quantile(mat_temp(:,13),0.25)],'b');
    plot([-90 90],[quantile(mat_temp(:,13),0.75) quantile(mat_temp(:,13),0.75)],'b');
    plot([-90 90],[quantile(mat_temp(:,13),0.10) quantile(mat_temp(:,13),0.10)],'k');
    plot([-90 90],[quantile(mat_temp(:,13),0.90) quantile(mat_temp(:,13),0.90)],'k');
    %plot([-90 90],[quantile(mat_temp(:,13),0.75)+1.5*(quantile(mat_temp(:,13),0.75) - quantile(mat_temp(:,13),0.25)) quantile(mat_temp(:,13),0.75)+1.5*(quantile(mat_temp(:,13),0.75) - quantile(mat_temp(:,13),0.25))],'k');
    %plot([-90 90],[quantile(mat_temp(:,13),0.25)-1.5*(quantile(mat_temp(:,13),0.75) - quantile(mat_temp(:,13),0.25)) quantile(mat_temp(:,13),0.25)-1.5*(quantile(mat_temp(:,13),0.75) - quantile(mat_temp(:,13),0.25))],'k');
    
    plot([quantile(mat_temp(:,12),0.50) quantile(mat_temp(:,12),0.50)],[-90 90],'r');
    plot([quantile(mat_temp(:,12),0.25) quantile(mat_temp(:,12),0.25)],[-90 90],'b');
    plot([quantile(mat_temp(:,12),0.75) quantile(mat_temp(:,12),0.75)],[-90 90],'b');
    plot([quantile(mat_temp(:,12),0.10) quantile(mat_temp(:,12),0.10)],[-90 90],'k');
    plot([quantile(mat_temp(:,12),0.90) quantile(mat_temp(:,12),0.90)],[-90 90],'k');
    %plot([quantile(mat_temp(:,12),0.75)+1.5*(quantile(mat_temp(:,12),0.75) - quantile(mat_temp(:,12),0.25)) quantile(mat_temp(:,12),0.75)+1.5*(quantile(mat_temp(:,12),0.75) - quantile(mat_temp(:,12),0.25))],[-90 90],'k');
    %plot([quantile(mat_temp(:,12),0.25)-1.5*(quantile(mat_temp(:,12),0.75) - quantile(mat_temp(:,12),0.25)) quantile(mat_temp(:,12),0.25)-1.5*(quantile(mat_temp(:,12),0.75) - quantile(mat_temp(:,12),0.25))],[-90 90],'k');
    plot([0 0],[-90 90],'k--');
    
    axis equal;
    xlim([-90 90]);ylim([-90 90]);
    xlabel('h version (deg)');ylabel('v version (deg)');
    cbar = colorbar;
    caxis([0 100]);
    
    
    subplot(1,3,2); hold on;
    scatter(mat_temp(:,11),mat_temp(:,12),5,mat_temp(:,15));
    plot([-90 90],[0 0],'k--');
%     plot([-90 90],[quantile(mat_temp(:,12),0.50) quantile(mat_temp(:,12),0.50)],'r');
%     plot([-90 90],[quantile(mat_temp(:,12),0.25) quantile(mat_temp(:,12),0.25)],'b');
%     plot([-90 90],[quantile(mat_temp(:,12),0.75) quantile(mat_temp(:,12),0.75)],'b');
%     plot([-90 90],[quantile(mat_temp(:,12),0.75)+1.5*(quantile(mat_temp(:,12),0.75) - quantile(mat_temp(:,12),0.25)) quantile(mat_temp(:,12),0.75)+1.5*(quantile(mat_temp(:,12),0.75) - quantile(mat_temp(:,12),0.25))],'k');
%     plot([-90 90],[quantile(mat_temp(:,12),0.25)-1.5*(quantile(mat_temp(:,12),0.75) - quantile(mat_temp(:,12),0.25)) quantile(mat_temp(:,12),0.25)-1.5*(quantile(mat_temp(:,12),0.75) - quantile(mat_temp(:,12),0.25))],'k');

    plot([0 0],[-90 90],'k--');
    axis equal;
    xlim([-90 90]);ylim([-90 90]);
    xlabel('vergence (deg)');ylabel('h version (deg)');
    cbar = colorbar;
    caxis([0 100]);
    

    subplot(1,3,3); hold on;
    scatter(mat_temp(:,11),mat_temp(:,13),5,mat_temp(:,15));
    plot([-90 90],[0 0],'k--');
%     plot([-90 90],[quantile(mat_temp(:,13),0.50) quantile(mat_temp(:,13),0.50)],'r');
%     plot([-90 90],[quantile(mat_temp(:,13),0.25) quantile(mat_temp(:,13),0.25)],'b');
%     plot([-90 90],[quantile(mat_temp(:,13),0.75) quantile(mat_temp(:,13),0.75)],'b');
%     plot([-90 90],[quantile(mat_temp(:,13),0.75)+1.5*(quantile(mat_temp(:,13),0.75) - quantile(mat_temp(:,13),0.25)) quantile(mat_temp(:,13),0.75)+1.5*(quantile(mat_temp(:,13),0.75) - quantile(mat_temp(:,13),0.25))],'k');
%     plot([-90 90],[quantile(mat_temp(:,13),0.25)-1.5*(quantile(mat_temp(:,13),0.75) - quantile(mat_temp(:,13),0.25)) quantile(mat_temp(:,13),0.25)-1.5*(quantile(mat_temp(:,13),0.75) - quantile(mat_temp(:,13),0.25))],'k');

    plot([0 0],[-90 90],'k--');
    axis equal;
    xlim([-90 90]);ylim([-90 90]);
    xlabel('vergence (deg)');ylabel('v version (deg)');
    cbar = colorbar;
    ylabel(cbar,'percent valid pixels');
    caxis([0 100]);
    
    print(fig1,'-depsc','-r300',[dst_dir subj1 '_' taskstr1 '_boxplot.eps']);
    print(fig2,'-depsc','-r300',[dst_dir subj1 '_' taskstr1 '_scatter.eps']);

    %vergence, version in time
    fig2 = figure(); hold on;
    subplot(3,1,1); hold on; title('vergence'); xlabel('time (s)'); ylabel('deg');
    plot((mat_temp(:,1) - mat_temp(1,1))./30, mat_temp(:,11),'k')
    subplot(3,1,2); hold on; title('h version'); xlabel('time (s)'); ylabel('deg');
    plot((mat_temp(:,1) - mat_temp(1,1))./30, mat_temp(:,12),'r')
    subplot(3,1,3); hold on; title('v version'); xlabel('time (s)'); ylabel('deg');
    plot((mat_temp(:,1) - mat_temp(1,1))./30, mat_temp(:,13),'b')
    print(fig2,'-depsc','-r300',[dst_dir subj1 '_' taskstr1 '_v_in_time.eps']);
    
    %overall eye movement velocity (deg/sec)
    x_vergence = gradient(mat_temp(:,11),.0333);
    x_hversion = gradient(mat_temp(:,12),.0333);
    x_vversion = gradient(mat_temp(:,13),.0333);
    
    %saccade and smooth pursuit eye movement velocity
    x_vergence_smooth = mat_temp(:,11);
    x_vergence_smooth(mat_temp(:,14) ~= 1) = NaN;
    x_vergence_smooth = gradient(x_vergence_smooth,.0333);
    
    x_vergence_sacc = mat_temp(:,11);
    x_vergence_sacc(mat_temp(:,14) ~= 3) = NaN;
    x_vergence_sacc = gradient(x_vergence_smooth,.0333);
    
    x_hversion_smooth = mat_temp(:,12);
    x_hversion_smooth(mat_temp(:,14) ~= 1) = NaN;
    x_hversion_smooth = gradient(x_hversion_smooth,.0333);
    
    x_hversion_sacc = mat_temp(:,12);
    x_hversion_sacc(mat_temp(:,14) ~= 3) = NaN;
    x_hversion_sacc = gradient(x_hversion_smooth,.0333);
    
    x_vversion_smooth = mat_temp(:,13);
    x_vversion_smooth(mat_temp(:,14) ~= 1) = NaN;
    x_vversion_smooth = gradient(x_vversion_smooth,.0333);
    
    x_vversion_sacc = mat_temp(:,13);
    x_vversion_sacc(mat_temp(:,14) ~= 3) = NaN;
    x_vversion_sacc = gradient(x_vversion_smooth,.0333);
    
    %make individual plots
    vvv_vel = [x_vergence ; x_hversion ; x_vversion];
    vvv_label = {[repmat({'vergence velocity'},length(x_vergence),1) ; repmat({'h version velocity'},length(x_hversion),1) ; repmat({'v version velocity'},length(x_vversion),1)]};
    
    vvv_vel_label = {[repmat({[num2str(quantile(x_vergence,[0.25 0.50]),'%.1f/') num2str(quantile(x_vergence,[0.75]),'%.1f') ...
        ' otlrs:' num2str(sum(x_vergence < quantile(x_vergence,0.25) - 1.5*(quantile(x_vergence,0.75) - quantile(x_vergence,0.25))) + sum(x_vergence > quantile(x_vergence,0.75) + 1.5*(quantile(x_vergence,0.75) - quantile(x_vergence,0.25)))) ...
        '/' num2str(sum(~isnan(x_vergence)))]},...
        length(x_vergence),1) ; ...
        repmat({[num2str(quantile(x_hversion,[0.25 0.50]),'%.1f/') num2str(quantile(x_hversion,[0.75]),'%.1f') ...
        ' otlrs:' num2str(sum(x_hversion < quantile(x_hversion,0.25) - 1.5*(quantile(x_hversion,0.75) - quantile(x_hversion,0.25))) + sum(x_hversion > quantile(x_hversion,0.75) + 1.5*(quantile(x_hversion,0.75) - quantile(x_hversion,0.25)))) ...
        '/' num2str(sum(~isnan(x_hversion)))]},...
        length(x_hversion),1) ; ...
        repmat({[num2str(quantile(x_vversion,[0.25 0.50]),'%.1f/') num2str(quantile(x_vversion,[0.75]),'%.1f') ...
        ' otlrs:' num2str(sum(x_vversion < quantile(x_vversion,0.25) - 1.5*(quantile(x_vversion,0.75) - quantile(x_vversion,0.25))) + sum(x_vversion > quantile(x_vversion,0.75) + 1.5*(quantile(x_vversion,0.75) - quantile(x_vversion,0.25)))) ...
        '/' num2str(sum(~isnan(x_vversion)))]},...
        length(x_vversion),1) ]};
    
    %vergence, version velocity box plots
    fig1 = figure(); hold on; ylabel('deg/sec'); title(regexprep([subj1 ' : ' taskstr1],'_',' '));
    h = boxplot(vvv_vel,{vvv_label vvv_vel_label});
    delete(h(end,:)); axis auto;
    plot([.5 3.5],[0 0],'k--');
    print(fig1,'-depsc','-r300',[dst_dir subj1 '_' taskstr1 '_velocity_boxplot.eps']);
    
    %calculate percentiles and outliers, make into labels
    verg_vel_label = [verg_vel_label ; repmat({[num2str(quantile(x_vergence,[0.25 0.50]),'%.1f/') num2str(quantile(x_vergence,[0.75]),'%.1f') ...
        ' otlrs:' num2str(sum(x_vergence < quantile(x_vergence,0.25) - 1.5*(quantile(x_vergence,0.75) - quantile(x_vergence,0.25))) + sum(x_vergence > quantile(x_vergence,0.75) + 1.5*(quantile(x_vergence,0.75) - quantile(x_vergence,0.25)))) ...
        '/' num2str(sum(~isnan(x_vergence)))]},...
        length(x_vergence),1)];
    hvers_vel_label = [hvers_vel_label ; repmat({[num2str(quantile(x_hversion,[0.25 0.50]),'%.1f/') num2str(quantile(x_hversion,[0.75]),'%.1f') ...
        ' otlrs:' num2str(sum(x_hversion < quantile(x_hversion,0.25) - 1.5*(quantile(x_hversion,0.75) - quantile(x_hversion,0.25))) + sum(x_hversion > quantile(x_hversion,0.75) + 1.5*(quantile(x_hversion,0.75) - quantile(x_hversion,0.25)))) ...
        '/' num2str(sum(~isnan(x_hversion)))]},...
        length(x_hversion),1)];
    vvers_vel_label = [vvers_vel_label ; repmat({[num2str(quantile(x_vversion,[0.25 0.50]),'%.1f/') num2str(quantile(x_vversion,[0.75]),'%.1f') ...
        ' otlrs:' num2str(sum(x_vversion < quantile(x_vversion,0.25) - 1.5*(quantile(x_vversion,0.75) - quantile(x_vversion,0.25))) + sum(x_vversion > quantile(x_vversion,0.75) + 1.5*(quantile(x_vversion,0.75) - quantile(x_vversion,0.25)))) ...
        '/' num2str(sum(~isnan(x_vversion)))]},...
        length(x_vversion),1)];
    
    %fixations
    fig3 = figure(); hold on;
    subplot(2,1,1); hold on; title('top view'); xlabel('z (m)'); ylabel('x (m)');
    %plot(mat_temp(:,4)./100, mat_temp(:,2)./100,'k.'); axis equal;
    scatter(mat_temp(:,4)./100,mat_temp(:,2)./100,15,mat_temp(:,15));
    cbar = colorbar; axis equal;
    ylabel(cbar,'percent valid pixels');
    caxis([0 100]);

    subplot(2,1,2); hold on; title('side view'); xlabel('z (m)'); ylabel('y (m)');
    %plot(mat_temp(:,4)./100, mat_temp(:,3)./100,'k.'); axis equal;
    scatter(mat_temp(:,4)./100,mat_temp(:,3)./100,15,mat_temp(:,15));
    cbar = colorbar; axis equal;
    ylabel(cbar,'percent valid pixels');
    caxis([0 100]);
    
    print(fig3,'-depsc','-r300',[dst_dir subj1 '_' taskstr1 '_fixations_m.eps']);
    
    %fixation distance histogram
    fig4 = figure(); hold on;
    scrsz = get(0,'ScreenSize');
    set(fig4,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
    subplot(2,2,1); hold on;
    hist(mat_temp(:,4)./100,logspace(log10(min(mat_temp(:,4)./100)),log10(max(mat_temp(:,4)./100))));
    set(gca,'XScale','log');
    xlabel({'fixation distance (m)','inf or greater: ' num2str(sum(isnan(mat_temp(:,4))))});
    
    subplot(2,2,2); hold on;
    probplot('lognormal',mat_temp(:,4)./100)
    
    subplot(2,2,3:4); hold on;
    hist(mat_temp(:,4)./100,100);
    xlabel({'fixation distance (m)','inf or greater: ' num2str(sum(isnan(mat_temp(:,4))))});
    
    print(fig4,'-depsc','-r300',[dst_dir subj1 '_' taskstr1 '_fixation_distances_hist.eps']);

    
    
%     %bin fixation counts to make 2d histogram of version/version frequency
%     bin_edges = linspace(-90,90,30);
%     num_bins = numel(bin_edges) - 1;
%     % hold counts
%     C = zeros(num_bins,num_bins);
%     % do counts
%     for i = 1:num_bins
%         for j = 1:num_bins
%             C(i,j) = length(find(mat_temp(:,12) >= bin_edges(j) & mat_temp(:,12) < bin_edges(j+1) & ...
%                 mat_temp(:,13) >= bin_edges(i) & mat_temp(:,13) < bin_edges(i+1)));
%         end
%     end
%     [x,y] = meshgrid(linspace(-90,90,num_bins),linspace(-90,90,num_bins));
%     %s1 = surfc(x,y,C);alpha(s1,'color');
%     [P,h] = contour(x,y,C,4,'color','k','linewidth',1);
%     clabel(P,h,'rotation',0);
    
    %print(fig3,'-depsc','-r300',[dst_dir subj1 '_' taskstr1 '_fixations_ecc.eps']);
    
    
    %2d histogram of fixation eccentricity
    fig4 = hist2d(mat_temp(:,12), mat_temp(:,13),50); axis equal; 

    print(fig4,'-depsc','-r300',[dst_dir subj1 '_' taskstr1 '_fixations_ecc_hist.eps']);

    close all;
    
end

if isempty(taskstr)
    %all vergence, version box plots
    fig1 = figure(); hold on;
    scrsz = get(0,'ScreenSize');
    set(fig1,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
    
    subplot(3,1,1); hold on; title('vergence'); ylabel('deg');
    h = boxplot(mat(:,11),{subj_label task_label verg_label});
    delete(h(end,:)); axis auto;
    plot([0 7],[0 0],'k--');
    subplot(3,1,2); hold on; title('h version'); ylabel('deg');
    h = boxplot(mat(:,12),{subj_label task_label hvers_label});
    delete(h(end,:)); axis auto;
    plot([0 7],[0 0],'k--');
    subplot(3,1,3); hold on; title('v version'); ylabel('deg');
    h = boxplot(mat(:,13),{subj_label task_label vvers_label});
    delete(h(end,:)); axis auto;
    plot([0 7],[0 0],'k--');
    
    set(fig1,'PaperPositionMode','auto');
    print(fig1,'-depsc','-r300',[dst_dir subj '_' taskstr '_boxplot.eps']);
    
    
    %all eye movement velocity plots
    x_vergence_all = gradient(mat(:,11),.0333);
    x_hversion_all = gradient(mat(:,12),.0333);
    x_vversion_all = gradient(mat(:,13),.0333);
    
    fig1 = figure(); hold on;
    scrsz = get(0,'ScreenSize');
    set(fig1,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
    
    subplot(3,1,1); hold on; title('vergence velocity'); ylabel('deg/sec');
    h = boxplot(x_vergence_all,{subj_label task_label verg_vel_label});
    delete(h(end,:)); axis auto;
    plot([0 7],[0 0],'k--');
    subplot(3,1,2); hold on; title('h version velocity'); ylabel('deg/sec');
    h = boxplot(x_hversion_all,{subj_label task_label hvers_vel_label});
    delete(h(end,:)); axis auto;
    plot([0 7],[0 0],'k--');
    subplot(3,1,3); hold on; title('v version velocity'); ylabel('deg/sec');
    h = boxplot(x_vversion_all,{subj_label task_label vvers_vel_label});
    delete(h(end,:)); axis auto;
    plot([0 7],[0 0],'k--');
    
    set(fig1,'PaperPositionMode','auto');
    print(fig1,'-depsc','-r300',[dst_dir subj '_' taskstr '_velocity_boxplot.eps']);
    
    %vergence/version 2d hist for all

    fig4 = hist2d(mat(:,12), mat(:,13),50); axis equal;     
    %scrsz = get(0,'ScreenSize');
    %set(fig4,'Position',[1 scrsz(4)/2 scrsz(3)/1.75 scrsz(4)]);
    %set(fig4,'PaperPositionMode','auto');
    print(fig4,'-depsc','-r300',[dst_dir subj '_' taskstr '_fixations_ecc_hist.eps']);
    
end


keyboard




