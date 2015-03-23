function [] = ET_frame_sync_drdp_eyemovements(filename)

%load in data
subj = filename;
load(['/Users/natdispstats/Documents/eyes/data_processing/data/' subj '/' subj '_proc.mat']);
subj = filename;



all_frames = [];

%intialize matrices for syncing
radials_100_1 = []; radials_50_1 = []; radials_450_1 = [];
walk_1 = [];
radials_450_2 = []; radials_100_2 = []; radials_50_2 = [];
walk_2 = [];
radials_50_3 = []; radials_100_3 = []; radials_200_3 = []; radials_450_3 = [];
walk_3 = [];
radials_450_4 = []; radials_100_4 = []; radials_50_4 = [];


%initialize matrices for target frames
radials_100_1_frames = []; radials_50_1_frames = []; radials_450_1_frames = [];

radials_450_2_frames = []; radials_100_2_frames = []; radials_50_2_frames = [];

radials_50_3_frames = []; radials_100_3_frames = []; radials_200_3_frames = []; radials_450_3_frames = [];

radials_450_4_frames = []; radials_100_4_frames = []; radials_50_4_frames = [];


%column vec of time stamps for each frame capture
frames = B(B(:,3) == 1);

%for each frame capture
for j = 1:length(frames)
    
    timestamp = frames(j);
    
    %grab indices of time points within 30ms window of fram capture
    ind = find(DM(:,1) >= timestamp - 15 & DM(:,1) <= timestamp + 15);
    
    %if there is at least 1 valid fixation within window
    if ~isempty(ind)
        
        %grab fixation location at the time point closest to shutter open
        ind = ind(abs(timestamp - DM(ind, 1 )) == min(abs(timestamp - DM(ind, 1 ))));
        
        fixation = DM(ind(1), 14:16 );
        href_le = DM(ind(1), 28:30 );
        href_re = DM(ind(1), 31:33);
        vergence = DM(ind(1), 21);
        version_h = DM(ind(1), 22);
        version_v = DM(ind(1), 23);
        
        %add data type flags
        data_flag = 1; %fixation
        
        SRdiff = SR-timestamp;
        SLdiff = SL-timestamp;
        if any(SRdiff(:,1) < 0 & SRdiff(:,2) > 0) || any(SLdiff(:,1) < 0 & SLdiff(:,2) > 0)
            data_flag = 3; %saccade
        end
        
    else
        fixation = [ NaN NaN NaN ];
        href_le = [ NaN NaN NaN ];
        href_re = [ NaN NaN NaN ];
        
        vergence = [NaN];
        version_h = [NaN];
        version_v = [NaN];
        
        %if there's no fixation location, classify as blink or missing data
        BRdiff = BR-timestamp;
        BLdiff = BL-timestamp;
        if any(BRdiff(:,1) < 0 & BRdiff(:,2) > 0) || any(BLdiff(:,1) < 0 & BLdiff(:,2) > 0)
            data_flag = 2; %blink
        else
            data_flag = 4; %missing data
        end
    end
    

    %save to data matrix
    all_frames(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag];
    
    if TG1(1,3) <= timestamp && timestamp < TG2(1,3)
       radials_100_1(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag];
    elseif TG1(2,3) <= timestamp && timestamp < TG2(2,3)
       radials_50_1(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag];
    elseif TG1(3,3) <= timestamp && timestamp < TG2(3,3)
       radials_450_1(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag]; 
    elseif TG1(4,3) <= timestamp && timestamp < TG2(4,3)
       radials_450_2(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag];
    elseif TG1(5,3) <= timestamp && timestamp < TG2(5,3)
       radials_100_2(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag];
    elseif TG1(6,3) <= timestamp && timestamp < TG2(6,3)
       radials_50_2(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag]; 
    elseif TG1(7,3) <= timestamp && timestamp < TG2(7,3)
       radials_50_3(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag]; 
    elseif TG1(8,3) <= timestamp && timestamp < TG2(8,3)
       radials_100_3(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag]; 
    elseif TG1(9,3) <= timestamp && timestamp < TG2(9,3)
       radials_200_3(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag];
    elseif TG1(10,3) <= timestamp && timestamp < TG2(10,3)
       radials_450_3(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag];   
    elseif TG1(11,3) <= timestamp && timestamp < TG2(11,3)
       radials_450_4(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag];   
    elseif TG1(12,3) <= timestamp && timestamp < TG2(12,3)
       radials_100_4(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag];    
    elseif TG1(13,3) <= timestamp && timestamp < TG2(13,3)
       radials_50_4(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag];   
    elseif TK1(1,1) <= timestamp && timestamp < TK2(1,1)
       walk_1(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag];  
    elseif TK1(2,1) <= timestamp && timestamp < TK2(2,1)
       walk_2(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag]; 
    elseif TK1(3,1) <= timestamp && timestamp < TK2(3,1)
       walk_3(end+1,:) = [j-1 fixation href_le href_re vergence version_h version_v data_flag];    
    end
end


save(['/Users/natdispstats/Documents/eyes/data_processing/data/' subj '/' subj '_eye_movements.mat'],'all_frames','radials_50_1','radials_50_2','radials_50_3','radials_50_4','radials_100_1','radials_100_2','radials_100_3','radials_100_4','radials_200_3','radials_450_1','radials_450_2','radials_450_3','radials_450_4','walk_1','walk_2','walk_3')

