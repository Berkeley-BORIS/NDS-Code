function [] = average_horopters

%HORIZONTAL HOROPTER DISPARITIES

%initialize H values from literature
H = [];
H = [.22 .24 .43];          %from Hillis & Banks 2001 (where accomodation effect was illiminated)
H = [H .076 .045 .05 .035]; %from Ogle farest viewings distance (76 cm) p34
H = [H .005];               %Ogle p45 (76cm)
H = [H .07 .1 .04];         %Helmholtz (71cm), Lau (150cm), and Lieberman (97cm) theifed from Ogle
H = [H .25 .36 -.11];       %from Schreiber et al.
H = [H .13 .07];            %Amigo 1967 (666 mm vd);

Have = mean(H);

%left eye azimuths
aLh      = -10:1:10; %azimuths to plot (based on vertical horopter data points)

%for each azimuth in the left eye, calculate azimuth in right eye for the average H value
aRh = acotd(cotd(aLh)+Have);

hh_disparity = aLh-aRh;

%invert sign because of image inversion in disparity data
hh_disparity = -hh_disparity;
hh_ecc = aLh;

save('average_horizontal_horopter.mat','hh_disparity','hh_ecc')

%%VERTICAL HOROPTER
%empty matrices
aLv = [];
aRv = [];

% %load in azimuth values for each subject, correcting for rotation matrix
% %conventions (all rotations are counterclockwise)
% counter = 1;
% for subj = {'AEL','ERB','GXK','HNS','JDB','KKD','LAT','MHA','MST','NSG','SCC','SPT','JMW','JAI','MDW','RYL','MEP','RXB','CCS','SAS','BFB','XMP','NBP','BTL','SBM','JJC','SEC','JMP'}
% load(['/Users/emily/Documents/CPointsv2/HeightData/correctedfits/Corrected-' num2str(cell2mat(subj)) '-MainExp-VD114-GA0-HA0-001.mat']);
%     for t = 1:length([D.shearValueDegN])
%             aLv(counter,t) = -D.shearValueDegN(t);
%             aRv(counter,t) = D.shearValueDegN(t);
%     end
%     counter = counter + 1;
% end
% 
% vh_disparity = mean(aLv - aRv,1);
% 
% %invert sign because of image inversion in disparity data
% 
% vh_disparity = -vh_disparity;

vh_ecc = [-8 -7 -6 -5 -4 -3 -2 2 3 4 5 6 7 8];
vh_disparity = [0.1264 0.1218 0.1550 0.1038 0.0894 0.0478 0.0464 -0.0475 -0.0765 -0.1303 -0.2106 -0.2677 -0.2734 -0.3235];

save('average_vertical_horopter.mat','vh_disparity','vh_ecc')



