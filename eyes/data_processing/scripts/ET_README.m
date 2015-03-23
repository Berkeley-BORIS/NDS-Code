%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ET_parse

%%%%WHAT DOES THIS FUNCTION DO?

%%%%OUTPUT

% D         = HREF fixation data (the important stuff, with blinks and saccades removed) (*timestamp* *left eye HREF x* *left eye HREF y* *left eye pupil area* *right eye HREF x* *right eye HREF y* *right eye pupil area*)
% D_orig    = HREF all data (nothing removed) (same as D)
% SL        = saccade event data (left eye) (*start timestamp* *end timestamp* not sure about the rest, need to check the python code and fill it in)
% SR        = saccade event data (right eye) (same as SL)
% FL        = fixation event data (left eye) (same as SL)
% FR        = fixation event data (right eye) (same as SL)
% BL        = blink event data (left eye) (same as SL)
% BR        = blink event data (right eye) (same as SL)
% T         = radial target data (MSG *timestamp* RADIAL_TARGET *view_dist* *slant(deg)* *tilt(deg)* *1 or 0 for start or end*)
% C         = cyclovergence data (MSG *timestamp* TORSION_RESULTS *view_dist* *slant(deg)* *tilt(deg* *torsion(half of full cyclovergence angle in deg))
% P         = fixation disparity data (MSG *timestamp* FD_RESULTS *view_dist* *horz_offset(deg)* *vertical_offset(deg)*)
% G         = luminance ramping data used to apply pupil size dependent gaze correction (MSG *timestamp* PUPIL_CORRECT *x pixels* *y pixels*)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ET_process

%%%%WHAT DOES THIS FUNCTION DO?

% (1) load the parsed session data (removed blinks, saccades, grabbed flags, etc)
% (2) checks for, grabs and separates out any repeated radial targets used for drift correction
% (3) calculate radial target locations (Targets)
% (4) calculate target locations for Pauling (PTargets)
% (5) converts arbitrary HREF coords to HREFCM coords in known coordinate system (Dcm)
% (6) do Paul pupil correction on Dcm and output various info to use to diagnose if this went well (DcmPaul, datPaul_LE, datPaul_RE, tarsPaulLE, tarPaulRE, paul_loo_err_LE, paul_loo_err_RE)
% (7) fit epipolar plane to gaze at each time point, so the eyes rays intersect (Dcmint)
% (8) build up a nice big data matrix (DM) with great information about gaze at each time point, like vergence, azimuth, 3D fixation pt, and more! (calculates torsion too)
% (9) calculate gaze and error for radial targets (DTargets, Derr)
% (10) if there are drift targets, calculate change in accuracy over time (Targetsd, Dtargetd, Derrd, DriftDat, DriftDatAll)
% (11) repeat steps 7-10 for Pauled data (DcmintPaul, DMPaul, DTargetsPaul, DerrPaul, DriftDatPaul, DriftDatAllPaul)
% (12) perform immediate drift correction using overall error in radial targets (DcmPaul2)
% (13) repeat steps 7-10 for Pauled and immediate drift corrected data (DcmintPaul2, DMPaul2, DTargetsPaul2, DerrPaul2, DriftDatPaul2, DriftDatAllPaul2)

%start work here, also, go over torsion
% (14) if there are drift targets, do the same for them
% (15) save workspace
% (16) display warnings

%%%%OUTPUT

%RADIAL TARGETS:

% T = target location data from EDF for initial radial targets (MSG *timestamp* RADIAL_TARGET *view_dist* *slant(deg)* *tilt(deg)* *1 or 0 for start or end*)
% Td = target location data from EDF for repeated drift correction radial targets

% Targets = thorough target location information 

% Columns = 
% (1) x ; (2) y ; (3) z(cm) (up, right, forward = positive) ; 
% (4) horiz vergence ; (5) horiz version ; (6) vert version (deg)
% (7) ; (8) left eye elevation and azimuth (helmholtz: down,left,clockwise = positive)
% (9) ; (10) right eye same
% (11) ; (12) tilt and slant of target from cyclopean eye (deg)
% (13) ; (14) cyclopean eye elevation and azimuth (deg)
% (15) ; (16) starting and ending timestamps


%PAUL TARGETS

% PTargets = target locations for pauling targets 

% Columns = 
% (1) x ; (2) y ; (3) z(cm) cyclopean ref (up, right, forward = positive)
% (4) x ; (5) y ; (6) z(cm) left eye ref (normed)
% (7) x ; (8) y ; (9) z(cm) right eye ref (normed)
% (10) ; (11) starting and ending timestamps

% DPTarget.LE/RE = gaze location and pupil size info for each paul target
% Columns =
% LEx(HREFCM) LEy LEpupil PTARGETx PTARGETy PTARGETindex(for each location/luminance combination)

% tarsPaul.LE/RE = x and y coords of paul targets in HREFCM

% paulfit. 
% p_LE, p_RE = structures of pupil sizes in column, each cell is a target location; 
% dx_LE, dy_LE, dx_RE, dy_RE = structures of of gaze error in column, each cell is a target location; 
% fx_LE, fy_LE, fx_RE, fy_RE = matrices of 3 polynomial coefficients for gaze error as a function of pupil size, each target is a row

% paul_loo_err_LE & RE (columns) for each timepoint of paul data (datPaul_LE and datPaul_RE)
% (1) original error between target and fixation(deg) 
% (2) same error after correction


%DATA MATRICES

% Dcm = D matrix with HREF converted to HREFCM, cyclopean ref 
% DcmPaul = Dcm, but pauled
% Dcmint = Dcm with plane adjustment

% Columns = 
% timestamp LEx LEy LEz REx REy REz

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
% (28) (29) HREFCM LE x y (unnormed, ie. projection of eye vector into plane at 100cm)
% (30) (31) HREFCM RE x y (unnormed)






