%%% Emily Cooper, Banks Lab, Nov 2011
%%% function to process eye tracking data (Emily, Nov 2011)
% as of 7/30/2012, this function uses D_orig, to maintain timepoints with
% saccades (used to only use fixation time points)
function [] = ET_process(filename)

% filename  = name of parsed mat file

%what does this function do?
% (1) load the parsed session data (removed blinks, saccades, grabbed flags, etc)
% (2) store target location info for immediate and drift targets
% (3) checks for, grabs and separates out any repeated radial targets used for drift correction
% (4) calculate radial target locations
% (5) converts arbitrary HREF coords to HREFCM coords in known coordinate system
% (6a) fit epipolar plane to gaze at each time point, so the eyes rays intersect
% (6b) build up a nice big data matrix with great information about gaze at each time point, like vergence, azimuth, 3D fixation pt, and more! (calculates torsion too)
% (6c) calculate gaze and error for radial targets
% (6d) if there are drift targets, calculate change in accuracy over time
% (7) perform immediate drift correction using overall error in radial targets
% (8) if there are ending drift targets, perform linear-in-time drift correction
% (9) save workspace
% (10) display warnings

%%NOTE: DM, DM2, DM3 are the primary output of this function. DM is the
%%data without any drift correction applied, DM2 and DM3 have some level of
%%drift corrections to them

%load in subject ipd in cm
subj_name = filename(1:3);
ipd = importdata(['/Users/natdispstats/Documents/eyes/ipds/' subj_name '.txt']);

href_dist = 100; %distance to calibration screen (cm) hardcoded now because it doesn't seem to change
L = [-ipd/2 0 0]; R = [ipd/2 0 0]; %translation vectors from head reference coordinate system to each eye
file_ext = 'proc'; %file extention to use when saving processed data
warnings = {}; %initialize structure to contain warnings
warning off all; %suppress warnings (the only warning you should get would be "matrix poorly conditioned" or "matrix almost singular." these are fine)
screen_cm = [121 68]; %screen dimensions in cm
screen_pix = [1920 1080]; %screen dimensions in pixels


%(1) LOAD IN PARSED DATA
tic
load(['/Users/natdispstats/Documents/eyes/data_processing/data/' filename '/' filename '.mat'])
display('file loaded')
toc
%D         = HREF fixation data (the important stuff, with blinks and saccades removed) (*timestamp* *left eye HREF x* *y* *pupil area* *right eye HREF x* *y* *pupil area*)
%D_orig    = HREF all data (nothing removed) (same as D) CURRENTLY WE USE D_orig FOR DATA PROCESSING TO RETAIN SACCADE TIME POINTS
%SL/SR     = saccade event data (left/right eyes) (*start timestamp* *endtimestamp* not sure about the rest, need to check the python code and fill it in)
%FL/FR     = fixation event data (same as SL/SR)
%BL/BR     = blink event data (same as SL/SR)
%T         = radial target data (MSG *timestamp* RADIAL_TARGET *view_dist* *slant(deg)* *tilt(deg)* *1 or 0 for start or end*)
%G         = luminance ramping data used to apply pupil size dependent gaze correction (MSG *timestamp* PUPIL_CORRECT *x pixels* *y pixels*)


%(2) STORE LOCATION INFO FOR ALL TARGETS, INCLUDING DRIFT
tic
[TargetsAll] = target_locations(T,ipd,L,R);
display('found all radial targets')
toc
%TargetsAll = thorough target location information 
% Columns = 
% (1) x ; (2) y ; (3) z(cm) (up, right, forward = positive) ; 
% (4) horiz vergence ; (5) horiz version ; (6) vert version (deg)
% (7) ; (8) NaNs
% (9) ; (10) NaNs
% (11) ; (12) tilt and slant of target from cyclopean eye (deg)
% (13) ; (14) NaNs
% (15) ; (16) starting and ending timestamps

%(3) CHECK FOR AND GRAB REPEATED TARGETS (indicates a drift correction)
tic
[T,Td] = drift_check(T,TG1,TG2);
display('checked for drift')
toc
%T = radial target location data for initial radial targets
%Td = radial target location data for repeated drift correction radial targets


%(4) FIND RADIAL TARGET LOCATIONS NON DRIFT
tic
[Targets] = target_locations(T,ipd,L,R);
display('found radial targets')
toc
%Targets = thorough target location information 
% Columns = 
% (1) x ; (2) y ; (3) z(cm) (up, right, forward = positive) ; 
% (4) horiz vergence ; (5) horiz version ; (6) vert version (deg)
% (7) ; (8) NaNs
% (9) ; (10) NaNs
% (11) ; (12) tilt and slant of target from cyclopean eye (deg)
% (13) ; (14) NaNs
% (15) ; (16) starting and ending timestamps


%(5) CONVERT HREF to HREFCM
% (a) find scaling between HREF and cm
% (b) find x,y translation between center of HREF plane and 0,0,href_dist in cm
tic
[Dcm] = convert_href(href_dist,Targets,T,D_orig);
display('data converted to cm')
toc
%Dcm = new D matrix with x,y,z of gaze in cm, in head reference (cyclopean eye) coordinates
% Columns:
% timestamp LEx(HREFCM) LEy LEz REx REy REz LEpupil REpupil


%(6) PROCESS DATA AND BUILD UP DATA MATRIX
% (a) MAKE EYE VECTORS INTERSECT
% (b) BUILD UP DATA MATRIX
% (c) GRAB EYE TRACKING DATA FOR EACH RADIAL TARGET AND CALCULATE ACCURACY
% (d) IF THERE ARE DRIFT TARGETS, CALCULATE ACCURACY OVER TIME
tic
[Dcmint,DM,Dtarget,Derr,DriftDat,Targetsd,Dtargetd,Derrd,warnings] = analysis_steps(ipd,href_dist,Dcm,Targets,Td,filename,warnings,L,R);
display('processed data')
toc
%Dcmint = Dcm with plane adjustment
%DM = full data matrix from Dcmint (columns)
% (1) timestamp 
% (2) (3) (4) HREFCM LE x y z (normed)
% (5) (6) (7) HREFCM RE x y z (normed)
% (8) (9) (10) EYEREFCM LE x y z (normed)
% (11) (12) (13) EYEREFCM RE x y z (normed)
% (14) (15) (16) FIXATION PT in HREFCM
% (17) (18) NANs
% (19) (20) NANs
% (21) horizontal vergence
% (22) horizontal version
% (23) vertical version
% (24) left eye pupil area
% (25) right eye pupil area
% (26) NANs
% (27) NANs
% (28) (29) (30) HREFCM LE x y z (unnormed, ie. projection of eye vector into plane at 100cm)
% (31) (32) (33) HREFCM RE x y z (unnormed)
%Dtarget = cell array of DM data associated with each target presentation
%Derr (columns)
% (1) target index
% (2) (3) (4) left, right, and mean of both eyes angular error between fixation and target
% (5) (6) (7) vergence, horizontal version, vertical version error
% (8) cm distance between fixation and target 
% (9) (10) NaNs
% (11) (12) NaNs
%DriftDat first row for radial targets, each additional row for each drift test
% (1) time(min)  (2) left eye error  (3) right eye error   (4) average both eyes error (deg)
%Targetsd = same as Targets, but for last drift targets only
%Dtargetsd = same are Dtargets, but for last drift targets only
%Derrd = Derr, except with one cell for each drift measurement


%(7) USE RADIAL TARGETS TO REDUCE OVERALL ERROR (IMMEDIATE DRIFT CORRECTION)
tic
[correction_le, correction_re, Dcm2] = immediate_drift(Targets,Dtarget,Dcm,href_dist,L,R);
%[correction_le, correction_re, Dcm2] = immediate_drift(Targets,Dtarget,Dcm,href_dist,L,R);
%Dcm2 = same as Dcm, but with first drift correction
[Dcmint2,DM2,Dtarget2,Derr2,DriftDat2,Targetsd2,Dtargetd2,Derrd2,warnings] = analysis_steps(ipd,href_dist,Dcm2,Targets,Td,filename,warnings,L,R);
display('processed first corrected data')
toc

%(8) IF THERE IS A ENDING RADIAL TARGETS SET, USE THIS TO APPLY LINEAR DRIFT
%CORRECTION OVER TIME AND REPEAT PROCESS DATA
tic
if ~isempty(Td)
[TargetsEND, DtargetEND, DerrEND, Dcm3, DriftDat3, Dcmint3, DM3, Dtarget3, Derr3, driftmat, warnings] ...
    = ending_drift(Td,ipd,L,R,DM,warnings,href_dist,correction_le,correction_re,Dcm,Dtarget,Targets,Targetsd,filename);
end
display('processed drift corrected data')
toc

%(9) SAVE WORKSPACE
save(['/Users/natdispstats/Documents/eyes/data_processing/data/' filename '/' filename '_' file_ext '.mat']);

%(10) DISPLAY WARNINGS
for j = 1:length(warnings)
    display(warnings{j})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Targets] = target_locations(T,ipd,L,R)
% Targets = thorough target location information 
% Columns = 
% (1) x ; (2) y ; (3) z(cm) (up, right, forward = positive) ; 
% (4) horiz vergence ; (5) horiz version ; (6) vert version (deg)
% (7) ; (8) NaNs
% (9) ; (10) NaNs
% (11) ; (12) tilt and slant of target from cyclopean eye (deg)
% (13) ; (14) NaNs
% (15) ; (16) starting and ending timestamps

Targets = [];

if ~isempty(T)
    %separate matrix for target starts and ends
    Tstart  = T(T(:,5) == 1,:);
    Tend    = T(T(:,5) == 0,:);

    %%for each target, determine coordinates, version/vergence, and az/el (helmholtz)
    for t = 1:length(Tstart(:,1))
        
        %xyz coordinates
        l = Tstart(t,2).*tand(Tstart(t,3)); %straight line distance from screen center to target
        x = l*cosd(Tstart(t,4));    Targets(t,1) = x; %convert to x and y using tilt
        y = l*sind(Tstart(t,4));    Targets(t,2) = y;
        z = Tstart(t,2);            Targets(t,3) = z; %z is just distance to screen
        F = [x y z]; %fixation pt vector

        % vergence
        Targets(t,4) = acosd( min(1,((L-F)/norm(L-F) * (R'-F')/norm(R'-F'))) ); %disconjugate angle between le and re vectors. (due to a floating point error in matlab, we don't let this value be less than one)
        
        % horizontal version
        N               = cross( L/norm(L), F/norm(F) ); % vector orthogonal to plane P
        N               = N / norm(N); %normalize
        A               = -cross(L/norm(L),N/norm(N)) + (N/norm(N) * L'/norm(L'))*N; % N rotated into plane containing fixation (epipolar plane)
        A               = A / norm(A); %normalize
        Targets(t,5)    = acosd( min(1,(A/norm(A) * F'/norm(F'))) ); %angle between A and fixation
        if( sign(F(1)*F(3)) < 0 )
            Targets(t,5) = -Targets(t,5); % leftward is negative
        end
        
        % vertical version
        Targets(t,6) = acosd( min(1, (A/norm(A) * [0 0 1]')) ); %angle between epipolar plane and z plane
        if( sign(F(2)*F(3)) < 0 )
            Targets(t,6) = -Targets(t,6); % down is negative
        end

        % azimuth and elevation in helmholtz - NAN        
        Targets(t,7) = NaN;
        Targets(t,8) = NaN;
        Targets(t,9) = NaN;
        Targets(t,10) = NaN;

        
        %tilt and slant
        Targets(t,11) = Tstart(t,4);
        Targets(t,12) = Tstart(t,3);
        
        %cyclopean eye elevation and azimuth - NAN
        Targets(t,13) = NaN;
        Targets(t,14) = NaN;
        
        
        %starting and ending time stamps
        Targets(t,15) = Tstart(t,1);
        Targets(t,16) = Tend(t,1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dcm] = convert_href(href_dist,Targets,T,D)

%Dcm 
% timestamp LEx LEy LEz REx REy REz LEpupil REpupil

%convert href to cms
href2cm  = (href_dist/15000); 
%convert HREF units to cm scaling
Dcm = [D(:,1) D(:,2).*href2cm D(:,3).*href2cm repmat(href_dist,length(D),1) D(:,5).*href2cm D(:,6).*href2cm repmat(href_dist,length(D),1) D(:,4) D(:,7)];

%To finish converting, I have to find the coordinates of the center of
%the HREF plane. I take the average coordinates for fixationat the center of the calibration screen

%find index of fixation at center of HREF plane (ie. target is at HREF distance and x and y equal zero
ind = find(Targets(:,3) == href_dist & Targets(:,1) == 0 & Targets(:,2) == 0);

%grab sample data between button down and button up
start_pos = T(ind,1);
start_delta = D(:,1) - start_pos;
start_delta( find(start_delta<0) ) = 1e9; % set negative values to large positive value
[minval,start_ind] = min( start_delta );

end_pos = T(ind+1,1);
end_delta = end_pos - D(:,1);
end_delta( find(end_delta<0) ) = 1e9; % set negative values to large positive value
[minval,end_ind] = min( end_delta );

%build matrix
target_center = D(start_ind:end_ind,:);

%estimate the origin of the HREF plane using the HREF coords of center fixation
%average for both eyes
center_x_HREF = median([([target_center(:,2)]) ; ([target_center(:,5)])]);
center_y_HREF = median([([target_center(:,3)]) ; ([target_center(:,6)])]);
center_x = center_x_HREF*href2cm;
center_y = center_y_HREF*href2cm;

%currenly unused, but it might be preferable to have a different center for each eye
le_center_x_HREF = median([target_center(:,2)]);
le_center_y_HREF = median([target_center(:,3)]);
le_center_x = le_center_x_HREF*href2cm;
le_center_y = le_center_y_HREF*href2cm;

re_center_x_HREF = median([target_center(:,5)]);
re_center_y_HREF = median([target_center(:,6)]);
re_center_x = re_center_x_HREF*href2cm;
re_center_y = re_center_y_HREF*href2cm;

%translate HREF coords to 0,0
Dcm(:,2) = (Dcm(:,2) - center_x); 
Dcm(:,3) = -(Dcm(:,3) - center_y);% flip U/D HREF coords to get into proper coordinate system (pos/neg: up/down, right/left)
Dcm(:,5) = (Dcm(:,5) - center_x);
Dcm(:,6) = -(Dcm(:,6) - center_y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DM] = build_data_mat(ipd,href_dist,Dcmint,filename,L,R)
%DM output (columns)
% (1) timestamp 
% (2) (3) (4) HREFCM LE x y z (normed)
% (5) (6) (7) HREFCM RE x y z (normed)
% (8) (9) (10) EYEREFCM LE x y z (normed)
% (11) (12) (13) EYEREFCM RE x y z (normed)
% (14) (15) (16) FIXATION PT in HREFCM
% (17) (18) EL and AZ of LE in HH
% (19) (20) EL and AZ of RE in HH
% (21) horizontal vergence
% (22) horizontal version
% (23) vertical version
% (24) left eye pupil area
% (25) right eye pupil area
% (26) TR of LE in HH
% (27) TR of RE in HH
% (28) (29) (30) HREFCM LE x y z (unnormed, ie. projection of eye vector into plane at 100cm)
% (31) (32) (33) HREFCM RE x y z (unnormed)

%ADD PUPIL AREA for LE and RE
Dcmint(:,24) = Dcmint(:,8);
Dcmint(:,25) = Dcmint(:,9);

%ADD EYEREF TO HREF COORDS
Dcmint(:,8:10) = bsxfun(@minus,Dcmint(:,2:4),L);
Dcmint(:,11:13) = bsxfun(@minus,Dcmint(:,5:7),R);

%normalize all vectors
n = sqrt(sum(Dcmint(:,2:4).^2,2)); n = n(:,ones(1,3));
Dcmint(:,2:4) = Dcmint(:,2:4)./n;
n = sqrt(sum(Dcmint(:,5:7).^2,2)); n = n(:,ones(1,3));
Dcmint(:,5:7) = Dcmint(:,5:7)./n;
n = sqrt(sum(Dcmint(:,8:10).^2,2)); n = n(:,ones(1,3));
Dcmint(:,8:10) = Dcmint(:,8:10)./n;
n = sqrt(sum(Dcmint(:,11:13).^2,2)); n = n(:,ones(1,3));
Dcmint(:,11:13) = Dcmint(:,11:13)./n;


%FIND INTERSECTION OF EYE VECTORS (FIXATION PT)
for k = 1:length(Dcmint)
    t = [Dcmint(k,8) -Dcmint(k,11);...
        Dcmint(k,9) -Dcmint(k,12);...
        Dcmint(k,10) -Dcmint(k,13)]...
        \[R(1) - L(1) ;  R(2) - L(2) ; R(3) - L(3)]; %matrix solution given eye direction vectors and eye translation vectors
    Dcmint(k,14) = L(1) + (Dcmint(k,8)*t(1)); %apply solution to left eye vector to get intersectio pt
    Dcmint(k,15) = L(2) + (Dcmint(k,9)*t(1));
    Dcmint(k,16) = L(3) + (Dcmint(k,10)*t(1));
end

%CALCULATE AZIMUTH AND ELEVATION IN HH - NANs
Dcmint(:,17) = NaN*ones(length(Dcmint(:,1)),1);
Dcmint(:,18) = NaN*ones(length(Dcmint(:,1)),1);
Dcmint(:,19) = NaN*ones(length(Dcmint(:,1)),1);
Dcmint(:,20) = NaN*ones(length(Dcmint(:,1)),1);

%CALCULATE VERGENCE AND VERSION BETWEEN EYES
LM = repmat(L./norm(L),length(Dcmint),1); % repeat left eye translation vector by length of data matric
n = sqrt(sum(Dcmint(:,14:16).^2,2)); n = n(:,ones(1,3)); %normalize fixation pt vector
F = Dcmint(:,14:16)./n; %normalized fixation pt vector

% VERGENCE
Dcmint(:,21) = acosd(min(1,dot(Dcmint(:,11:13),Dcmint(:,8:10),2))); %vergence angle between two vectors (bug in matlab makes some normed vectors slightly greater than one. cap it at one to avoid this)
Dcmint((Dcmint(:,16) < 0),21) = -Dcmint((Dcmint(:,16) < 0),21); %vergence is negative if rays diverge
% HORIZONTAL VERSION
N = cross( LM, F,2 ); % vector orthogonal to epipolar plane
n = sqrt(sum(N.^2,2)); n = n(:,ones(1,3)); %normalize vector
N = N./n;
A = -cross(LM,N,2) + bsxfun(@times, dot(N,LM,2), N); % N rotated into plane
n = sqrt(sum(A.^2,2)); n = n(:,ones(1,3)); %normalize vector
A = A./n;
Dcmint(:,22) = acosd( min(1,dot(A,F,2))); %version
Dcmint(((F(:,1).*F(:,3)) < 0),22) = -Dcmint(((F(:,1).*F(:,3)) < 0),22); %give version direction
% VERTICAL VERSION
Dcmint(:,23) = acosd( min(1,dot(A,repmat([0 0 1],length(Dcmint),1),2)));
Dcmint(((F(:,3)) < 0),23) = 180 - Dcmint(((F(:,3)) < 0),23); %correct version angle when eyes diverge
Dcmint(((F(:,2).*F(:,3)) < 0),23) = -Dcmint(((F(:,2).*F(:,3)) < 0),23); %give version up/down direction


%IF RAYS DIVERGE, WE DON'T WANT A FIXATION PT. REMOVE IT IN THIS CASE
Dcmint(Dcmint(:,16) < 0,14) = NaN;
Dcmint(Dcmint(:,16) < 0,15) = NaN;
Dcmint(Dcmint(:,16) < 0,16) = NaN;

%rename
DM = Dcmint;

%%ADD TORSION IN HELMHOLTZ COORDINATES
DM(:,26) = NaN*ones(length(DM(:,1)),1);
DM(:,27) = NaN*ones(length(DM(:,1)),1);

%add on unnormed HREFCM for LE and RE (ie. projection of each eye's vector into a plante at 100cm)
href_le = bsxfun(@times,(href_dist./DM(:,10)),DM(:,8:10));
DM(:,28:30) = bsxfun(@plus,href_le,[-ipd/2 0 0]);
href_re = bsxfun(@times,(href_dist./DM(:,13)),DM(:,11:13));
DM(:,31:33) = bsxfun(@plus,href_re,[ipd/2 0 0]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dtarget,Derr,warnings] = depth_test(ipd,Targets,DM,warnings,L,R)

Dtarget = [];
Derr = [];

if size(Targets,1) > 1
    Derr = [];
    for k = 1 : length(Targets(:,1))
        %grab sample data between button down and button up
        start_pos = Targets(k,15);
        start_delta = DM(:,1) - start_pos;
        start_delta( find(start_delta<0) ) = 1e9; % set negative values to large positive value
        [minval,start_ind] = min( start_delta );

        end_pos = Targets(k,16);
        %end_pos = start_pos+500;
        end_delta = end_pos - DM(:,1);
        end_delta( find(end_delta<0) ) = 1e9; % set negative values to large positive value
        [minval,end_ind] = min( end_delta );

        %build matrix
        Dtarget{k} = DM(start_ind:end_ind,:);
        
        %BUILD DATA MATRIX WITH ERROR METRICS
        % pure anglular error between target and fixation for le and re
        le_ang   = acosd(min(1,dot(Dtarget{k}(:,8:10),repmat([Targets(k,1:3) - L]./norm([Targets(k,1:3) - L]),size(Dtarget{k},1),1),2)));
        re_ang   = acosd(min(1,dot(Dtarget{k}(:,11:13),repmat([Targets(k,1:3) - R]./norm([Targets(k,1:3) - R]),size(Dtarget{k},1),1),2)));
        mean_ang = mean([le_ang re_ang],2);
        % vergence and version
        hverg    = bsxfun(@minus,Dtarget{k}(:,21),Targets(k,4));
        hvers    = bsxfun(@minus,Dtarget{k}(:,22),Targets(k,5));
        vvers    = bsxfun(@minus,Dtarget{k}(:,23),Targets(k,6));
        % cm distance
        dist     = bsxfun(@minus,Dtarget{k}(:,14:16),Targets(k,1:3));
        dist     = sqrt(dist(:,1).^2 + dist(:,2).^2 + dist(:,3).^2);
        
        %error in helmholtz elevation and azimuth - NaNs
        nan1 = NaN*ones(length(dist),1);
        nan2 = NaN*ones(length(dist),1);
        nan3 = NaN*ones(length(dist),1);
        nan4 = NaN*ones(length(dist),1);

        Derr = [Derr ; repmat(k,size(Dtarget{k},1),1) le_ang re_ang mean_ang hverg hvers vvers dist nan1 nan2 nan3 nan4];
    end
end
%Derr
% (1) target index
% (2) (3) (4) left, right, and mean of both eyes angular error between fixation and target
% (5) (6) (7) vergence, horizontal version, vertical version error
% (8) cm distance between fixation and target 
% (9) (10) NaNs
% (11) (12) NaNs

%look for missing data during target fixations and issue a warning if found
skip = 0;
for j = 1:length(Dtarget)
    if isempty(Dtarget{j})
        w = ['missing radial target: ' num2str(Targets(j,1:3))];
        %don't repeat warnings
        if isempty(warnings)
            warnings{end+1} = w;
        else
            for k = 1:length(warnings)
                if strcmp(warnings{k}, w)
                    skip = 1;
                end
            end
            if skip == 0
                warnings{end+1} = w;
            elseif skip == 1
                skip = 0;
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dcmint,DM,Dtarget,Derr,DriftDat,Targetsd,Dtargetd,Derrd,warnings] = analysis_steps(ipd,href_dist,Dcm,Targets,Td,filename,warnings,L,R)

% (7) MAKE VECTORS INTERSECT
tic
[Dcmint] = do_epipolar(Dcm, ipd,L,R);
display('epipolar planed')
toc

% (8) BUILD UP FIXATION DATA MATRIX
tic
[DM] = build_data_mat(ipd,href_dist,Dcmint,filename,L,R);
display('built data matrix')
toc

% (9) GRAB EYE TRACKING DATA FOR EACH RADIAL TARGET AND CALCULATE ACCURACY
tic
[Dtarget,Derr,warnings] = depth_test(ipd,Targets,DM,warnings,L,R);
display('calculated target info')
toc

% (10) IF THERE ARE DRIFT TARGETS CALCULATE DRIFT IN ACCURACY OVER TIME
tic
if ~isempty(Td)
    [DriftDat, Targetsd, Dtargetd, Derrd, warnings] = analyze_drift(Td,ipd,DM,Targets,Dtarget,href_dist, warnings, L, R, Derr);
else
    DriftDat = [];
    Targetsd = [];
    Dtargetd = [];
    Derrd = [];
end
display('drift analyzed')
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,Td]  = drift_check(T,TG1,TG2)
% T = radial target data for initial radial targets
% Td = radial target data for repeated drift correction radial targets

tmp = T;
Td = [];
cnt = 0;
for t = 1:length(unique(TG1(:,2)))
    Target_starts = TG1(TG1(:,2) == t,3);
    Target_ends = TG2(TG2(:,2) == t,3);
    if t == 1
        T = tmp(tmp(:,1) >= Target_starts(1) & tmp(:,1) <= Target_ends(end),:);
    else
        cnt = cnt + 1;
        Td{cnt} = tmp(tmp(:,1) >= Target_starts(1) & tmp(:,1) <= Target_ends(end),:);
    end
end

%Td = T(T(:,1) > TG2(3,2),:);
%T = T(T(:,1) >= TG1(1,2) & T(:,1) <= TG2(3,2),:);

%Td = [];%intialize a new Target matrix to fill up with drift targets


%nt = length(TG1)/length(unique(TG1(:,1))); %are there repeated targets at all?
%if nt > 1 %if so...
%    tarnum = size(unique(T(:,2:end),'rows'),1); %number of unique targets
%    temp = T; %temporary copy of T
%    T = temp(1:tarnum,:); %original target matrix contains just first instances of unique targets
%    temp = temp(tarnum+1:end,:); %drift targets = everything else
%    ntd = length(temp)/length(unique(temp(:,2:end),'rows')); %how many repeats are there in the drift targets?
%    tarnumd = size(unique(temp(:,2:end),'rows'),1); %how many targets are being repeats?
%    for m = 1:ntd %fill up drift target matrix with repeated targets
%        Td(:,:,m) = temp(1:tarnumd,:);
%        temp = temp(tarnumd + 1:end,:);
%    end
%end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p2_new,q2_new] = fitplane(ipd,p2,q2,L,R)

%%% DATA
p1 = L; % define translation
p2 = p2 - p1; %define direction
p2 = p2 / norm(p2); % make unit vector from direction
p2 = p2 + p1; % translate back

q1 = R;
q2 = q2 - q1; %define direction
q2 = q2 / norm(q2); % make unit vector from direction
q2 = q2 + q1; %translate back


%%% SOLVE FOR PLANE THROUGH p1 AND p2, WHILE MINIMIZING DISTANCE TO
%%% q1 AND q2: CONSTRAINED LEAST-SQUARES (solve: y = bz)
X      = [p1(1) p2(1) q1(1) q2(1)];
Y      = [p1(2) p2(2) q1(2) q2(2)];
Z      = [p1(3) p2(3) q1(3) q2(3)];
A(:,1) = Z';
u      = inv(A'*A)*A'*Y';
u      = [0 u 0]; % equation of plane y = ax + bz + c


%%% PROJECT VECTORS ONTO PLANE
p2_new  = [p2(1) u(2)*p2(3) p2(3)]; % project
q2_new  = [q2(1) u(2)*q2(3) q2(3)]; % project

q2_new  = q2_new - q1; % translate to origin
q2_new  = q2_new / norm(q2_new); % make unit vector
q2_new  = q2_new + q1; % translate 

p2_new  = p2_new - p1; % translate to origin
p2_new  = p2_new / norm(p2_new); % make unit vector
p2_new  = p2_new + p1; % translate 

P       = [u(1) -1 u(2)]; % vector normal to plane
P       = P / norm(P);
thp     = asind( P*p2' ); % angle between original 
thq     = asind( P*q2' ); %  vector and plane
thp_new = asind( P*p2_new' ); % angle between original 
thq_new = asind( P*q2_new' ); %  vector and plane


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [DriftDat, Targetsd, Dtargetd, Derrd, warnings] = analyze_drift(Td,ipd,DM,Targets,Dtarget,href_dist, warnings, L, R, Derr)

%start drift structure with error for initial radial targets
DriftDat(1,:) = [0 nanmedian(Derr(:,2)) nanmedian(Derr(:,3)) nanmedian(Derr(:,4))];

%if there are drift measurements
if ~isempty(Td)
    %for each drift measurements
    for d = 1:length(Td)
        %get target info
        [Targetsd] = target_locations(Td{d},ipd,L,R);
        %get gaze and error info for targets
        [Dtargetd,Derrd_temp,warnings] = depth_test(ipd,Targetsd,DM,warnings,L,R);
        %save error into structure
        Derrd{d} = Derrd_temp;
        %fill in drift error info
        DriftDat(end+1,:) = [(mean(Dtargetd{1}(:,1)) - mean(Dtarget{1}(:,1)))/60000 nanmedian(Derrd{d}(:,2)) nanmedian(Derrd{d}(:,3)) nanmedian(Derrd{d}(:,4))];
        %DriftDat, first row for radial targets, each additional row for each drift test
        % Columns:
        % (1) time(min)
        % (2) median left eye error (deg)
        % (3) median right eye error (deg)
        % (4) median both eyes error (deg)
        
        %look for missing data during target fixations and issue a warning if found
        for j = 1:length(Dtargetd)
            if isempty(Dtargetd{j})
                warnings{end+1} = ['missing drift radial target: ' num2str(Targetsd(j,1:3)) ' for drift #' num2str(d)];
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dcmint] = do_epipolar(Dcm,ipd,L,R)
%assume that vertical vergence is zero, and any vertical difference between the two eye vectors 
%is due to error in the system. To remove it, project the two vectors vertically onto the best-fit plane 
%(the plane is constrained to contain the interocular axis, ie. x axis)

Dcmint(:,1) = Dcm(:,1);
for j = 1:length(Dcm)
    [Dcmint(j,2:4), Dcmint(j,5:7)] = fitplane(ipd,Dcm(j,2:4),Dcm(j,5:7),L,R);
end

Dcmint(:,8) = Dcm(:,8);
Dcmint(:,9) = Dcm(:,9);

% Dcmint = Dcm with plane adjustment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [correction_le, correction_re, DcmCor] = immediate_drift(Targets,DtargetOrig,DcmOrig,href_dist,L,R)

%intialize
gaze_href_le = []; gaze_href_re = []; tar_loc_href = [];

%for each target
for p = 1:length(Targets(:,1))
    %only use targets at href distance for correction
    if Targets(p,3) == href_dist
        %take mean x and y locations of each eyes gaze within href plane (cm)
        gaze_href_le(end+1,:) = [mean(DtargetOrig{p}(:,28)) mean(DtargetOrig{p}(:,29))];
        gaze_href_re(end+1,:) = [mean(DtargetOrig{p}(:,31)) mean(DtargetOrig{p}(:,32))];
        
        %store target location in href plane
        tar_loc_href(end+1,:) = Targets(p,1:2);
    end
end

%calculate average x,y error between gaze and target for each eye
correction_le = nanmean(gaze_href_le - tar_loc_href);
correction_re = nanmean(gaze_href_re - tar_loc_href);

%use targets at all distances to apply correction
if(1)
    %initialize data vectors
    gaze_href_le = []; %HREFCM for left and right eyes (projection of each eye vector into plan 100cm away)
    gaze_href_re = [];
    tar_loc_href_le = []; %HREFCM for left and right targets
    tar_loc_href_re = [];
    %Targets_temp = [];
    
    %for each target
    for p = 1:length(Targets(:,1))
        %take mean x and y locations of each eyes gaze with href plane
        gaze_href_le(end+1,:) = [mean(DtargetOrig{p}(:,28)) mean(DtargetOrig{p}(:,29))];
        gaze_href_re(end+1,:) = [mean(DtargetOrig{p}(:,31)) mean(DtargetOrig{p}(:,32))];
        
        %calculate intersection of each predicted target vector with href plane
        tar_loc_href_le_temp = bsxfun(@minus,Targets(p,1:3),L);
        tar_loc_href_le_temp = bsxfun(@times,(href_dist./tar_loc_href_le_temp(3)),tar_loc_href_le_temp);
        tar_loc_href_le(end+1,:) = bsxfun(@plus,tar_loc_href_le_temp,L);
        
        tar_loc_href_re_temp = bsxfun(@minus,Targets(p,1:3),R);
        tar_loc_href_re_temp = bsxfun(@times,(href_dist./tar_loc_href_re_temp(3)),tar_loc_href_re_temp);
        tar_loc_href_re(end+1,:) = bsxfun(@plus,tar_loc_href_re_temp,R);
        
        %Targets_temp(end+1,:) = Targets(p,:);
    end
    %calculate average x,y error between vector and target for each eye
    correction_le = nanmean(gaze_href_le - tar_loc_href_le(:,1:2));
    correction_re = nanmean(gaze_href_re - tar_loc_href_re(:,1:2));
end


%apply correction to HREFCM data matrix
DcmCor = bsxfun(@minus,DcmOrig,[0 correction_le 0 correction_re 0 0 0]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Targetsd, DtargetOrigEND, DerrOrigEND, DcmCor, DriftDat3, Dcmint3, DM3, Dtarget3, Derr3, driftmat, warnings] ...
    = ending_drift(Td,ipd,L,R,DM,warnings,href_dist,correction_le,correction_re,DcmOrig,DtargetOrig,Targets,Targetsd,filename)

if ~isempty(Td)
    
    %initialize
    gaze_href_leEND = []; gaze_href_reEND = []; tar_loc_hrefEND = [];
    
    %grab gaze data only during last targets
    [DtargetOrigEND,DerrOrigEND,warnings] = depth_test(ipd,Targetsd,DM,warnings,L,R);
    
    %for each target
    for p = 1:length(Targetsd(:,1))
        %only use targets at href distance for correction
        if Targetsd(p,3) == href_dist
            %take mean x and y locations of each eyes gaze within href plane (cm)
            gaze_href_leEND(end+1,:) = [mean(DtargetOrigEND{p}(:,28)) mean(DtargetOrigEND{p}(:,29))];
            gaze_href_reEND(end+1,:) = [mean(DtargetOrigEND{p}(:,31)) mean(DtargetOrigEND{p}(:,32))];
            
            %store target location in href plane
            tar_loc_hrefEND(end+1,:) = Targetsd(p,1:2);
        end
    end
    
    %average error between gaze and target in x and y
    correction_leEND = nanmean(gaze_href_leEND - tar_loc_hrefEND);
    correction_reEND = nanmean(gaze_href_reEND - tar_loc_hrefEND);
    
    %regression
    %amt of time between first and last targets
    driftmat.run = (mean(DtargetOrigEND{1}(:,1)) - mean(DtargetOrig{1}(:,1))); 
    
    driftmat.le_driftx = (correction_leEND(1) - correction_le(1))./driftmat.run;
    driftmat.le_offsetx = correction_le(1);
    driftmat.le_drifty = (correction_leEND(2) - correction_le(2))./driftmat.run;
    driftmat.le_offsety = correction_le(2);
    
    driftmat.re_driftx = (correction_reEND(1) - correction_re(1))./driftmat.run;
    driftmat.re_offsetx = correction_re(1);
    driftmat.re_drifty = (correction_reEND(2) - correction_re(2))./driftmat.run;
    driftmat.re_offsety = correction_re(2);
    
    driftmat.le_end = correction_leEND;
    driftmat.re_end = correction_reEND;
    
    le_Dx_adj = ((DcmOrig(:,1) - mean(DtargetOrig{1}(:,1))).*driftmat.le_driftx) + driftmat.le_offsetx;
    le_Dy_adj = ((DcmOrig(:,1) - mean(DtargetOrig{1}(:,1))).*driftmat.le_drifty) + driftmat.le_offsety;
    re_Dx_adj = ((DcmOrig(:,1) - mean(DtargetOrig{1}(:,1))).*driftmat.re_driftx) + driftmat.re_offsetx;
    re_Dy_adj = ((DcmOrig(:,1) - mean(DtargetOrig{1}(:,1))).*driftmat.re_drifty) + driftmat.re_offsety;
    
    
    DcmCor = DcmOrig;
    DcmCor(:,2) = DcmOrig(:,2) - le_Dx_adj;
    DcmCor(:,3) = DcmOrig(:,3) - le_Dy_adj;
    DcmCor(:,5) = DcmOrig(:,5) - re_Dx_adj;
    DcmCor(:,6) = DcmOrig(:,6) - re_Dy_adj;
    
    
    [Dcmint3,DM3,Dtarget3,Derr3,DriftDat3,Targetsd3,Dtargetd3,Derrd3,warnings] = analysis_steps(ipd,href_dist,DcmCor,Targets,Td,filename,warnings,L,R);
end