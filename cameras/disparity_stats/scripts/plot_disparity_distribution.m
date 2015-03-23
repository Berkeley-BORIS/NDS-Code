function [] = plot_disparity_distribution(mat,color)

%bws_nature_walk_1 = load('../data/bwsnat5/bwsnat5_task_nature_walk_1_disparity.mat');

%cut down to 8deg radius circle

degperpixel = 2*atand(.5/583);
rad = ceil(8./degperpixel);

mat = mat(104-82:104+82,104-82:104+82,:);

y = 1:size(mat,1);
x = 1:size(mat,2);
[x y] = meshgrid(x,y);

[row col] = find((y-83).^2 + (x-83).^2 > rad.^2);

for t = 1:length(row)
   mat(row(t),col(t),:) = NaN;
end

%normalize each frame
for f = 1:size(mat,3)
    mat(:,:,f) = mat(:,:,f)-mat(83,83,f);
end

mat = mat(:);
mat = mat(~isnan(mat));

[n x] = hist(mat,100);
%bar(x,n./trapz(x,n),color);
xi = linspace(min(mat),max(mat),201);
[f,xi] = ksdensity(mat);
plot(xi,f,color,'linewidth',1)



