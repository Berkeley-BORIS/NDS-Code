import numpy as np
import os
import math
import bottleneck as bn
import  matplotlib.pyplot as plt
import matplotlib.cm as cm
import cv
import scipy.io as sio
import sys

def create_vertical_disparity_stats_matrix(subj,activity,src_dir,dst_dir,start_frame,duration_sec):
    '''take in disparity values from frames in a directory and calculate disparity distribution
    '''

    #create destination directory if it doesn't exist
    if not os.path.exists(dst_dir):
        os.mkdir(dst_dir)

    #intrinsics
    focal_length_pix = 583

    #mat_size = 1842

    #convert duration in milliseconds to a end frame
    num_frames = int((duration_sec*1000)/33.333)
    end_frame = start_frame + num_frames

    #iterate through source directory of disparity maps and grab desired frames
    imlist = os.listdir(src_dir)
    frame_num_list = np.ones((len(imlist),1))*np.nan
    cnt = 0
    for im in imlist:
        frame_num_list[cnt] = int(im.split('_')[2].strip('.npy'))
        cnt += 1

    if start_frame > 0:
        #imlist = np.asarray(imlist)[np.where(frame_num_list >= start_frame)[0]].tolist()
        imlist = np.asarray(imlist)[np.where(np.logical_and(frame_num_list >= start_frame, frame_num_list <=  end_frame))[0]].tolist()

    mat_size = len(imlist)
    print 'number of frames = ' + str(mat_size)

    if mat_size > 5000:
        print 'mat is toooo big!'
        sys.exit(1)

    #allocate data matrix for disparity frames
    disparity_mat = np.ones((207,207,mat_size))*np.nan

    #start counter
    cnt = 0
    #load in each disparity map for averaging
    for im in imlist:
        #load disparity map
        disparity_frame = np.load(src_dir + im)

        #if not np.isnan(np.nanmax(disparity_frame)):
        current_im_string = "%d of %d (%s)" % (cnt, mat_size, im)
        update_progress_bar(current_im_string)
        #convert disparities from pixels to degrees
        disparity_frame = (2*np.arctan2((disparity_frame/2),focal_length_pix))*(180/math.pi)

        disparity_mat[:,:,cnt] = disparity_frame
        cnt += 1


    #save full disparity matrix w/ 32 bit precision
    disparity_mat = disparity_mat.astype('float32', casting='same_kind')
    np.savez_compressed(dst_dir + subj + '_' + activity + '_' +'vertical_disparity.npz',disparity=disparity_mat)
    sio.savemat(dst_dir + subj + '_' + activity + '_' + 'vertical_disparity.mat', {'disparity':disparity_mat}, do_compression=True)

def update_progress_bar(update_string):
    sys.stdout.write('\r' + update_string)
    sys.stdout.flush()