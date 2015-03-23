function [] = downsample_disparity_mat(filename)

subjname = strtok(filename,'_');
loadpath = ['../data/' subjname '/' filename];
savepath = ['../data/' subjname '/' strtok(filename,'.') '_downsampled.mat'];
disparity = load(loadpath);

disparity = disparity.disparity;

%allocate space for downsampled disparity maps
disparity_mat = NaN*ones(69,69,size(disparity,3));

%create grid of indices for every third pixel
[x,y] = meshgrid(2:3:207,2:3:207);

%create gaussian filter to sample around every third pixel
filter = fspecial('gaussian',5,(2/3));

%for m = 1:length(disparity)
for m = 1:size(disparity,3)
    
    
    im = disparity(:,:,m);
    im2 = NaN*ones(69,69);
    
    %downsample by a factor of 3, taking the weighted average of 5x5 blocks of
    %pixels, ignoring nans, so long as at least 3 valid pixels are found in
    %a 3x3 central sub-block
    for j = 1 : length(y)
        for k = 1:length(x)
            %if this is an edge pixel without a 5x5 block around it, use
            %only 3x3 block for averaging
            if j == 1 || k == 1 || j == 69 || k == 69
                blk5 = im(y(j,k)-1:y(j,k)+1,x(j,k)-1:x(j,k)+1);
                blk5 = cat(2,[NaN ; NaN ; NaN], blk5, [NaN ; NaN ; NaN]);
                blk5 = cat(1,[NaN NaN NaN NaN NaN], blk5, [NaN NaN NaN NaN NaN]);
                ind5 = find( ~isnan(blk5) );
            else
                %5x5 block for filtering
                blk5 = im(y(j,k)-2:y(j,k)+2,x(j,k)-2:x(j,k)+2);
                ind5 = find( ~isnan(blk5) );
            end
            
            %3x3 block for nan thresholding
            blk3 = im(y(j,k)-1:y(j,k)+1,x(j,k)-1:x(j,k)+1);
            ind3 = find( ~isnan(blk3) );
             
            if( length(ind3) >= 3 )
               im2(j,k) = nansum( filter(ind5)/sum(filter(ind5)).*blk5(ind5));
            end
        end
    end

%     if m == 15
%         
%         fig1 = figure(); hold on;
%         max_abs_disparity = max([abs(max(max(im))) abs(min(min(im)))]);
%         sc(flipud([im]),'diff',[-max_abs_disparity max_abs_disparity],[0 0 0]);
%         cbar = colorbar;
%         ylabel(cbar,'disparity (deg)');
%         
%         print(fig1,'-depsc','-r300',['~/Desktop/9415notsmoothed.eps']);
% 
%         
%         keyboard
%         
%     end
    disparity_mat(:,:,m) = im2;

end

keyboard

save(savepath,'disparity_mat');

