%interpolate missing pixels and add 10 deg circle to eye images
function [] = interp_pixels(data_dir)

color_circle = 0;

dst_dir = [data_dir 'interped_ana/'];
mkdir(dst_dir);
list_files = [data_dir '*.png'];
imfiles = dir(list_files);

R = 103^2;


    
%for k = 1:1
for f = 1:length(imfiles)
    filename = [data_dir imfiles(f).name];
    %im1 = flipdim(imread(filename),1);
    im1 = imread(filename);
    if size(im1,3) == 3
        im = im1(:,:,3);
    else
        im = im1(:,:,1);
    end

    % replace NaNs
    im2 = im;
    [Y,X] = find( im == 0 );
    for k = 1 : length(X)
        if Y(k) > 2 && Y(k) < 549 && X(k) > 2 && X(k) < 549
            blk = im(Y(k)-2:Y(k)+2,X(k)-2:X(k)+2);
            ind = find( blk ~= 0 );
            if( length(ind) >= 3 )
                im2(Y(k),X(k)) = sum( blk(ind) )/length(ind);
            end
        end
    end
    
    im3 = imadjust(im2,[.01 .99],[]);
    
    if color_circle
        im2 = cat(3,im2,im2,im2);
        im3 = cat(3,im3,im3,im3);
        
        im3(268:284,275:277,1) = 255;
        im3(268:284,275:277,2) = 0;
        im3(268:284,275:277,3) = 0;
        
        im3(275:277,268:284,1) = 255;
        im3(275:277,268:284,2) = 0;
        im3(275:277,268:284,3) = 0;
        
        
        for j = 1:550
            for k = 1:550
                if round((j-275).^2 + (k-275).^2) <= R+200 && round((j-275).^2 + (k-275).^2) >= R-200
                    im3(k,j,1) = 255;
                    im3(k,j,2) = 0;
                    im3(k,j,3) = 0;
                end
            end
        end
    else
        im3(268:284,275:277) = 0;  
        im3(275:277,268:284) = 0;
        
        %crop to 15 deg square
        im3 = im3(25:526,25:526);
        %if red
        %    im3 = cat(3,im3,zeros(314,314),zeros(314,314));
        %elseif green
        %    im3 = cat(3,zeros(314,314),im3,zeros(314,314));
        %end
            

    end
    
    imwrite(im3,[dst_dir imfiles(f).name]);
    % display original and interpolated
    %imagesc( [im3] ); axis image off;
end

keyboard