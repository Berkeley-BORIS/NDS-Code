function [] = calc_and_save_stats_by_gaze_ecc(data_mat_path)

data = load(data_mat_path);
subj = strtok(fliplr(strtok(fliplr(data_mat_path),'/')),'_');
load(['../../../eyes/data_processing/data/' subj '/' subj '_eye_movements.mat']);

[x] = textread(['..' strtok(data_mat_path,'.') '_params.txt'], ...
    '%*s %*s %s',3);

start_frame = str2num(cell2mat(x(2)));
end_frame = str2num(cell2mat(x(3)));

if strfind(subj,'drdp')
    if strfind(data_mat_path,'walk_1')
        fixations_mat = walk_1(start_frame+1:end_frame+1,:);
    elseif strfind(data_mat_path,'walk_2')
        fixations_mat = walk_2(start_frame+1:end_frame+1,:);
    elseif strfind(data_mat_path,'walk_3')
        fixations_mat = walk_3(start_frame+1:end_frame+1,:);
    end
else
    fixations_mat = task(start_frame+1:end_frame+1,:);
end

fixations_mat = fixations_mat(~isnan(fixations_mat(:,11)),:);

datamat.vergence.bins = [-inf -11.25:2.5:11.25 inf];
versionbins = [-inf -17.5:7:17.5 inf];
%datamat.version.hbins = repmat(versionbins(1:end-1),1,11);
%datamat.version.vbins = [repmat(versionbins(1),1,11) repmat(versionbins(2),1,11) repmat(versionbins(3),1,11) ...
%    repmat(versionbins(4),1,11) repmat(versionbins(5),1,11) repmat(versionbins(6),1,11) ...
%    repmat(versionbins(7),1,11) repmat(versionbins(8),1,11) repmat(versionbins(9),1,11) repmat(versionbins(10),1,11) ...
%    repmat(versionbins(11),1,11)];
datamat.vergence.data = {};
datamat.version.data = {};

for vg = 1:length(datamat.vergence.bins)-1
    datamat.vergence.data{vg} = [];
end
for vs = 1:(length(versionbins)-1)^2
    datamat.version.data{vs} = [];
end

%datacnt = 0;
    
display(['processing ' data_mat_path]);
tic
for n = 1:size(fixations_mat,1)
    
    for vg = 1:length(datamat.vergence.bins)-1
        if fixations_mat(n,11) > datamat.vergence.bins(vg) && fixations_mat(n,11) <= datamat.vergence.bins(vg+1)
            datamat.vergence.data{vg}(:,:,end+1) = data.disparity(:,:,n);
        end
    end
    vscnt = 1;
    
    for vvs = 1:length(versionbins)-1
        for hvs = 1:length(versionbins)-1
            if fixations_mat(n,12) > versionbins(hvs) && ...
                    fixations_mat(n,12) <= versionbins(hvs+1) && ...
                    fixations_mat(n,13) > versionbins(vvs) && ...
                    fixations_mat(n,13) <= versionbins(vvs+1)
                datamat.version.data{vscnt}(:,:,end+1) = data.disparity(:,:,n);
                %datacnt = datacnt +1;
                %display(datacnt)
            end
            datamat.version.vbins(vscnt) = versionbins(vvs);
            datamat.version.hbins(vscnt) = versionbins(hvs);
            vscnt = vscnt+1;
        end
    end
   % if datacnt == 0
   %     keyboard
   % end
end

%clear out empty first bin
for vg = 1:length(datamat.vergence.bins)-1
    if sum(sum(datamat.vergence.data{vg}(:,:,1))) == 0
        datamat.vergence.data{vg} = datamat.vergence.data{vg}(:,:,2:end);
    end
end
for vs = 1:vscnt-1
    if sum(sum(datamat.version.data{vs}(:,:,1))) == 0
        datamat.version.data{vs} = datamat.version.data{vs}(:,:,2:end);
    end
end

for vg = 1:length(datamat.vergence.bins)-1
    if ~isempty(datamat.vergence.data{vg})
        for j = 1:207
            for k = 1:207
                datamat.vergence.quartile1_mat{vg}(j,k) = quantile(datamat.vergence.data{vg}(j,k,~isnan(datamat.vergence.data{vg}(j,k,:))),0.25,3);
                datamat.vergence.quartile3_mat{vg}(j,k) = quantile(datamat.vergence.data{vg}(j,k,~isnan(datamat.vergence.data{vg}(j,k,:))),0.75,3);
            end
        end
        datamat.vergence.median_mat{vg} = nanmedian(datamat.vergence.data{vg},3);
        datamat.vergence.var_mat{vg} = datamat.vergence.quartile3_mat{vg} - datamat.vergence.quartile1_mat{vg};
        datamat.vergence.count_mat{vg} = sum(~isnan(datamat.vergence.data{vg}),3);
        datamat.vergence.percent_mat{vg} = sum(~isnan(datamat.vergence.data{vg}),3)./size(datamat.vergence.data{vg},3);
        
        datamat.vergence.mean_mat{vg} = nanmean(datamat.vergence.data{vg},3);
        datamat.vergence.std_mat{vg} = nanstd(datamat.vergence.data{vg},0,3);
        datamat.vergence.sterr_mat{vg} = datamat.vergence.std_mat{vg}./sqrt(datamat.vergence.count_mat{vg});
    else
        datamat.vergence.quartile1_mat{vg} = [];
        datamat.vergence.quartile3_mat{vg} = [];
        datamat.vergence.median_mat{vg} = []; 
        datamat.vergence.var_mat{vg} = [];
        datamat.vergence.count_mat{vg} = [];
        datamat.vergence.percent_mat{vg} = [];
        datamat.vergence.mean_mat{vg} = [];
        datamat.vergence.std_mat{vg} = [];
        datamat.vergence.sterr_mat{vg} = [];
    end
end

for vs = 1:length(datamat.version.hbins)
    if ~isempty(datamat.version.data{vs})
        for j = 1:207
            for k = 1:207
                datamat.version.quartile1_mat{vs}(j,k) = quantile(datamat.version.data{vs}(j,k,~isnan(datamat.version.data{vs}(j,k,:))),0.25,3);
                datamat.version.quartile3_mat{vs}(j,k) = quantile(datamat.version.data{vs}(j,k,~isnan(datamat.version.data{vs}(j,k,:))),0.75,3);
            end
        end
        datamat.version.median_mat{vs} = nanmedian(datamat.version.data{vs},3);
        datamat.version.var_mat{vs} = datamat.version.quartile3_mat{vs} - datamat.version.quartile1_mat{vs};
        datamat.version.count_mat{vs} = sum(~isnan(datamat.version.data{vs}),3);
        datamat.version.percent_mat{vs} = sum(~isnan(datamat.version.data{vs}),3)./size(datamat.version.data{vs},3);
        
        datamat.version.mean_mat{vs} = nanmean(datamat.version.data{vs},3);
        datamat.version.std_mat{vs} = nanstd(datamat.version.data{vs},0,3);
        datamat.version.sterr_mat{vs} = datamat.version.std_mat{vs}./sqrt(datamat.version.count_mat{vs});
    else
        datamat.version.quartile1_mat{vs} = [];
        datamat.version.quartile3_mat{vs} = [];
        datamat.version.median_mat{vs} = []; 
        datamat.version.var_mat{vs} = [];
        datamat.version.count_mat{vs} = [];
        datamat.version.percent_mat{vs} = [];
        datamat.version.mean_mat{vs} = [];
        datamat.version.std_mat{vs} = [];
        datamat.version.sterr_mat{vs} = [];
    end
    
end

clear data.disparity;
datamat.vergence.data = [];
datamat.version.data = [];
save(strrep(data_mat_path,'.mat','_stats_by_gaze_ecc.mat'),'datamat');
toc

