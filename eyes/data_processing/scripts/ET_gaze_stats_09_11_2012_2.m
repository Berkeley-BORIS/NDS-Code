function [] = ET_gaze_stats(subj,taskstr)
%subj = 'bwsnat5'
%task = 'task_nature_walk_1'

dst_dir = ['../data/' subj '/plots/'];
mkdir(dst_dir);

if strcmp('bws',subj) && isempty(taskstr)
    tasks = [{'bwsnat5','task_nature_walk_1'},{'bwsnat5','task_campus_walk'},...
        {'bwsure1','task_nature_walk_2'},{'bwsdrdp6','task_walk_1'},...
        {'bwscafe1','task_ordering_coffee'},{'bwssand1','task_sandwich'}];
elseif strcmp('ges',subj) && isempty(taskstr)
    tasks = [{'gesnat1','task_nature_walk_1'},...
        {'gesure1','task_nature_walk_2'},{'gesdrdp1','task_walk_1'},...
        {'gescafe2','task_ordering_coffee'},{'gessand1','task_sandwich'}];
end
if isempty(taskstr)
    mat = [];
    subj_label = [];
    task_label = [];
    verg_label = []; verg_vel_label = [];
    hvers_label = []; hvers_vel_label = [];
    vvers_label = []; vvers_vel_label = [];
    
    for t = 1:2:length(tasks)-1
        subj1 = tasks{t};
        taskstr1 = tasks{t+1};
        load(['../data/' subj1 '/' subj1 '_eye_movements.mat']);
        vvv = [];
        vvv_label = []; vvv_meas_label = [];
        
        [x] = textread(['/users/natdispstats/Documents/cameras/disparity_stats/data/' subj1 '/' subj1 '_' taskstr1 '_disparity_params.txt'], ...
            '%*s %*s %s',3);
        
        start_frame = str2num(cell2mat(x(2)));
        end_frame = str2num(cell2mat(x(3)));
        
        if strfind(subj1,'drdp')
            mat_temp = walk_1(start_frame+1:end_frame+1,:);
        else
            mat_temp = task(start_frame+1:end_frame+1,:);
        end
        mat = [mat ; mat_temp];
        
        %calculate percentiles and outliers, make into labels
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
        
        %make individual plots
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
        
        %vergence, version box plots
        fig1 = figure(); hold on; ylabel('deg'); title(regexprep([subj1 ' : ' taskstr1],'_',' '));
        h = boxplot(vvv,{vvv_label vvv_meas_label});
        delete(h(end,:)); axis auto;
        plot([.5 3.5],[0 0],'k--');
        print(fig1,'-depsc','-r300',[dst_dir subj1 '_' taskstr1 '_boxplot.eps']);
        
        %vergence, version in time
        fig2 = figure(); hold on;
        subplot(3,1,1); hold on; title('vergence'); xlabel('time (s)'); ylabel('deg');
        plot((mat_temp(:,1) - mat_temp(1,1))./30, mat_temp(:,11),'k')
        subplot(3,1,2); hold on; title('h version'); xlabel('time (s)'); ylabel('deg');
        plot((mat_temp(:,1) - mat_temp(1,1))./30, mat_temp(:,12),'r')
        subplot(3,1,3); hold on; title('v version'); xlabel('time (s)'); ylabel('deg');
        plot((mat_temp(:,1) - mat_temp(1,1))./30, mat_temp(:,13),'b')
        print(fig2,'-depsc','-r300',[dst_dir subj1 '_' taskstr1 '_v_in_time.eps']);
        
        %eye movement velocity (deg/sec)
        x_vergence = gradient(mat_temp(:,11),.0333);
        x_hversion = gradient(mat_temp(:,12),.0333);
        x_vversion = gradient(mat_temp(:,13),.0333);
        
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
        fig1 = figure(); hold on; ylabel('deg'); title(regexprep([subj1 ' : ' taskstr1],'_',' '));
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
    end
    
    close all;
    
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
    
    subplot(3,1,1); hold on; title('vergence velocity'); ylabel('deg');
    h = boxplot(x_vergence_all,{subj_label task_label verg_vel_label});
    delete(h(end,:)); axis auto;
    plot([0 7],[0 0],'k--');
    subplot(3,1,2); hold on; title('h version velocity'); ylabel('deg');
    h = boxplot(x_hversion_all,{subj_label task_label hvers_vel_label});
    delete(h(end,:)); axis auto;
    plot([0 7],[0 0],'k--');
    subplot(3,1,3); hold on; title('v version velocity'); ylabel('deg');
    h = boxplot(x_vversion_all,{subj_label task_label vvers_vel_label});
    delete(h(end,:)); axis auto;
    plot([0 7],[0 0],'k--');
    
    set(fig1,'PaperPositionMode','auto');
    print(fig1,'-depsc','-r300',[dst_dir subj '_' taskstr '_velocity_boxplot.eps']);

    
else
    
    load(['../data/' subj '/' subj '_eye_movements.mat']);
    
    [x] = textread(['/users/natdispstats/Documents/cameras/disparity_stats/data/' subj '/' subj '_' taskstr '_disparity_params.txt'], ...
        '%*s %*s %s',3);
    
    start_frame = str2num(cell2mat(x(2)));
    end_frame = str2num(cell2mat(x(3)));
    
    if strfind(subj,'drdp')
        mat = walk_1(start_frame+1:end_frame+1,:);
    else
        mat = task(start_frame+1:end_frame+1,:);
    end
    
   % mat = task(task(:,1) >= start_frame & task(:,1) <= end_frame,:);
    vvv = [mat(:,11) ; mat(:,12) ; mat(:,13)];
    vvv_label = {[repmat({'vergence'},length(mat(:,11)),1) ; repmat({'h version'},length(mat(:,12)),1) ; repmat({'v version'},length(mat(:,13)),1)]};  
    
    vvv_meas_label = {[repmat({[num2str(quantile(mat(:,11),[0.25 0.50]),'%.1f/') num2str(quantile(mat(:,11),[0.75]),'%.1f') ...
            ' otlrs:' num2str(sum(mat(:,11) < quantile(mat(:,11),0.25) - 1.5*(quantile(mat(:,11),0.75) - quantile(mat(:,11),0.25))) + sum(mat(:,11) > quantile(mat(:,11),0.75) + 1.5*(quantile(mat(:,11),0.75) - quantile(mat(:,11),0.25)))) ...
            '/' num2str(sum(~isnan(mat(:,11))))]},...
            length(mat(:,11)),1) ; ...
            repmat({[num2str(quantile(mat(:,11),[0.25 0.50]),'%.1f/') num2str(quantile(mat(:,11),[0.75]),'%.1f') ...
            ' otlrs:' num2str(sum(mat(:,11) < quantile(mat(:,11),0.25) - 1.5*(quantile(mat(:,11),0.75) - quantile(mat(:,11),0.25))) + sum(mat(:,11) > quantile(mat(:,11),0.75) + 1.5*(quantile(mat(:,11),0.75) - quantile(mat(:,11),0.25)))) ...
            '/' num2str(sum(~isnan(mat(:,11))))]},...
            length(mat(:,11)),1) ; ...
            repmat({[num2str(quantile(mat(:,11),[0.25 0.50]),'%.1f/') num2str(quantile(mat(:,11),[0.75]),'%.1f') ...
            ' otlrs:' num2str(sum(mat(:,11) < quantile(mat(:,11),0.25) - 1.5*(quantile(mat(:,11),0.75) - quantile(mat(:,11),0.25))) + sum(mat(:,11) > quantile(mat(:,11),0.75) + 1.5*(quantile(mat(:,11),0.75) - quantile(mat(:,11),0.25)))) ...
            '/' num2str(sum(~isnan(mat(:,11))))]},...
            length(mat(:,11)),1) ]};
end

if ~isempty(taskstr)
    
    %vergence, version box plots
    fig1 = figure(); hold on; ylabel('deg'); title(regexprep([subj ' : ' taskstr],'_',' '));
    h = boxplot(vvv,{vvv_label vvv_meas_label});
    delete(h(end,:)); axis auto;
    plot([.5 3.5],[0 0],'k--');
    print(fig1,'-depsc','-r300',[dst_dir subj '_' taskstr '_boxplot.eps']);

    %vergence, version in time
    fig2 = figure(); hold on;
    subplot(3,1,1); hold on; title('vergence');
    plot(mat(:,1) - start_frame, mat(:,11),'k')
    subplot(3,1,2); hold on; title('h version');
    plot(mat(:,1) - start_frame, mat(:,12),'r')
    subplot(3,1,3); hold on; title('v version');
    plot(mat(:,1) - start_frame, mat(:,13),'b')
    print(fig2,'-depsc','-r300',[dst_dir subj '_' taskstr '_v_in_time.eps']);
    
    keyboard
end

keyboard
%fixations
fig3 = figure(); hold on;
subplot(2,1,1); hold on; title('top view');
plot(mat(:,4), mat(:,2),'k.')
subplot(2,1,2); hold on; title('side view');
plot(mat(:,4), mat(:,3),'k.')
print(fig3,'-depsc','-r300',[dst_dir 'fixations.eps']);


%change in eye angle over time
x_vergence = gradient(mat(:,11),.0333);
x_hversion = gradient(mat(:,12),.0333);
x_vversion = gradient(mat(:,13),.0333);
fig4 = figure(); hold on; 
subplot(3,1,1); hold on; title('vergence deg/sec');
plot(mat(:,1) - start_frame, x_vergence,'k')
subplot(3,1,2); hold on; title('h version deg/sec');
plot(mat(:,1) - start_frame, x_hversion,'r')
subplot(3,1,3); hold on; title('v version deg/sec');
plot(mat(:,1) - start_frame, x_vversion,'b')
print(fig4,'-depsc','-r300',[dst_dir 'speed.eps']);

%vergence, version sped hist
fig1 = figure(); hold on;
subplot(3,1,1); hold on; title('vergence');
hist(x_vergence)
subplot(3,1,2); hold on; title('h version');
hist(x_hversion)
subplot(3,1,3); hold on; title('v version');
hist(x_vversion)
print(fig1,'-depsc','-r300',[dst_dir 'histograms_speed.eps']);

