function [] = ET_gaze_stats(subj,taskstr)
%subj = 'bwsnat5'
%task = 'task_nature_walk_1'

dst_dir = ['../data/' subj '/plots/'];
mkdir(dst_dir);

if strcmp('bws',subj) && isempty(taskstr)
    tasks = [{'bwsnat5','task_nature_walk_1'},{'bwsnat5','task_campus_walk'},...
            {'bwsure1','task_nature_walk_2'},{'bwsdrdp6','task_walk_1'},...
            {'bwssand1','task_sandwich'},{'bwscafe1','task_ordering_coffee'}];

    mat = [];
    vvv = []; vvv_label = [];
    for t = 1:2:length(tasks)-1
        subj1 = tasks{t};
        taskstr1 = tasks{t+1};
        load(['../data/' subj1 '/' subj1 '_eye_movements.mat']);

        [x] = textread(['/users/natdispstats/Documents/cameras/disparity_stats/data/' subj1 '/' subj1 '_' taskstr1 '_disparity_params.txt'], ...
            '%*s %*s %s',3);

        start_frame = str2num(cell2mat(x(2)));
        end_frame = str2num(cell2mat(x(3)));

        if strfind(subj1,'drdp')
           % mat_temp = walk_1(walk_1(:,1) >= start_frame & walk_1(:,1) <= end_frame,:);
            mat_temp = walk_1(start_frame+1:end_frame+1,:);
        else
            %mat_temp = task(task(:,1) >= start_frame & task(:,1) <= end_frame,:);
            mat_temp = task(start_frame+1:end_frame+1,:);
        end
        mat = [mat ; mat_temp];
        
        vvv = [vvv ; mat_temp(:,11) ; mat_temp(:,12) ; mat_temp(:,13)];
        vvv_label = [vvv_label ; {[repmat({'vergence'},length(mat_temp(:,11)),1) ; ...
            repmat({'h version'},length(mat_temp(:,12)),1) ;...
            repmat({'v version'},length(mat_temp(:,13)),1)] ...
            repmat({subj1},3*length(mat_temp(:,11)),1)}];
    end
    c1 = [vvv_label{:,1}]; c1 = c1(:);
    c2 = [vvv_label{:,2}]; c2 = c2(:);
    b = cell(1,2);
    b{1,1} = c1;
    b{1,2} = c2;
    vvv_label = b;
else
    
    load(['../data/' subj '/' subj '_eye_movements.mat']);
    
    [x] = textread(['/users/natdispstats/Documents/cameras/disparity_stats/data/' subj '/' subj '_' taskstr '_disparity_params.txt'], ...
        '%*s %*s %s',3);
    
    start_frame = str2num(cell2mat(x(2)));
    end_frame = str2num(cell2mat(x(3)));
    
    mat = task(task(:,1) >= start_frame & task(:,1) <= end_frame,:);
    vvv = [mat(:,11) ; mat(:,12) ; mat(:,13)];
    vvv_label = {[repmat({'vergence'},length(mat(:,11)),1) ; repmat({'h version'},length(mat(:,12)),1) ; repmat({'v version'},length(mat(:,13)),1)] ...
        repmat({subj},3*length(mat(:,11)),1)};
    %[j-1 fixation href_le href_re vergence version_h version_v data_flag];
    
end

%vergence, version box plots
fig1 = figure(); hold on;
boxplot(vvv,vvv_label);
print(fig1,'-depsc','-r300',[dst_dir subj '_' taskstr '_boxplot.eps']);


if ~isempty(taskstr)
    %vergence, version in time
    fig2 = figure(); hold on;
    subplot(3,1,1); hold on; title('vergence');
    plot(mat(:,1) - start_frame, mat(:,11),'k')
    subplot(3,1,2); hold on; title('h version');
    plot(mat(:,1) - start_frame, mat(:,12),'r')
    subplot(3,1,3); hold on; title('v version');
    plot(mat(:,1) - start_frame, mat(:,13),'b')
    print(fig2,'-depsc','-r300',[dst_dir subj '_' taskstr '_v_in_time.eps']);
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

