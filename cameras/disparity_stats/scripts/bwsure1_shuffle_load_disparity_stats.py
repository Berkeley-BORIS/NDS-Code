from create_disparity_stats_matrix import *

subj = 'bwsure1'
flag = '_shuffle'
activity = 'task_nature_walk_2'
src_dir = '../../3d_reprojection/data/' + subj + flag + '/' + subj + '_' + activity + '_eye_images_transform1/10deg/disparityeye/'
dst_dir = '../data/' + subj + flag + '/'
start_frame = 14974
duration_sec = 120

#average_disparities(src_dir=src_dir,dst_dir=dst_dir,frames=frames)
create_disparity_stats_matrix(subj=subj,activity=activity,src_dir=src_dir,dst_dir=dst_dir,start_frame=start_frame,duration_sec=duration_sec)