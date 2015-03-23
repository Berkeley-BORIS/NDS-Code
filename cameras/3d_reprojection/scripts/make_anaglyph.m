function [] = make_anaglyph(left_path,right_path)


dst_dir = [left_path 'anaglyph/'];
mkdir(dst_dir);
list_files_left = [left_path '*.png'];
imfiles_left = dir(list_files_left);

list_files_right = [right_path '*.png'];
imfiles_right = dir(list_files_right);

for f = 1:length(imfiles_left)
    filename_left = [left_path imfiles_left(f).name];
    im_left = imread(filename_left);
    
    filename_right = [right_path imfiles_right(f).name];
    im_right = imread(filename_right);
    
    %im_ana = cat(3,im_left,im_right,im_right);
    im_ana = cat(3,im_right,im_left,im_left);
    
    imwrite(im_ana,[dst_dir imfiles_left(f).name]);
end