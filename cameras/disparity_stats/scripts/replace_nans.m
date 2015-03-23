clear;

% test image
if(0)
    N = 206;
    R = 102^2;
    im = ones(N,N);
    ind = round( N^2*rand(1,1500) );
    im(ind) = NaN;
else
    load disparity1.mat;
    im = disparity;
    N = size(im,1);
    R = (N/2-2)^2;
end

% replace NaNs
im2 = im;
[Y,X] = find( isnan(im) );
for k = 1 : length(X)
    if( ((X(k)-N/2).^2 + (Y(k)-N/2).^2) <= R )
        blk = im(Y(k)-1:Y(k)+1,X(k)-1:X(k)+1);
        ind = find( ~isnan(blk) );
        if( length(ind) >= 4 )
            im2(Y(k),X(k)) = sum( blk(ind) )/length(ind);
        end
    end
end

% display original and interpolated
imagesc( [im im2] ); axis image off;