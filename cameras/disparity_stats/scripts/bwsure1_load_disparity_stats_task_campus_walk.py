from create_disparity_stats_matrix import *

subj = 'bwsure1'
activity = 'task_campus_walk'
src_dir = '../../3d_reprojection/data/' + subj + '/' + subj + '_' + activity + '_eye_images_transform1/10deg/disparityeye/'
dst_dir = '../data/' + subj + '/'
start_frame = 4851
duration_sec = 120

#average_disparities(src_dir=src_dir,dst_dir=dst_dir,frames=frames)
create_disparity_stats_matrix(subj=subj,activity=activity,src_dir=src_dir,dst_dir=dst_dir,start_frame=start_frame,duration_sec=duration_sec)