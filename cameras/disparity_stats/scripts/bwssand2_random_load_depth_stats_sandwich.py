import os
from create_depth_stats_matrix import *

def find_start_frame(src_dir):
    """
    Finds the first frame number
    """
    imlist = os.listdir(src_dir)
    imlist.sort()
    imname_components = imlist[0].strip('.npy').split('_')

    return int(imname_components[-1])

subj = 'bwssand2'
fixation_flag = '_random'
if fixation_flag and not fixation_flag.startswith('_'):
    fixation_flag = '_' + fixation_flag

activity = 'task_sandwich'

src_tail = '3d_reprojection/data/' + subj + fixation_flag + '/' + subj + '_' + activity + '_eye_images_transform1/deptheye/'
src_hd1 = '../../'
src_hd2 = '/Volumes/Macintosh HD 2/cameras2'
if os.path.exists(os.path.join(src_hd1, src_tail)):
    src_dir = os.path.join(src_hd1, src_tail)
else:
    src_dir = os.path.join(src_hd2, src_tail)
dst_dir = '../data/' + subj + fixation_flag + '/'
start_frame = find_start_frame(src_dir)
duration_sec = 120

#average_disparities(src_dir=src_dir,dst_dir=dst_dir,frames=frames)
create_depth_stats_matrix(subj=subj,activity=activity,src_dir=src_dir,dst_dir=dst_dir,start_frame=start_frame,duration_sec=duration_sec)