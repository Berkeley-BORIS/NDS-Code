import numpy as np
import os
import math
import bottleneck as bn
import  matplotlib.pyplot as plt
import matplotlib.cm as cm
import cv
import scipy.io as sio
import sys
from glob import glob

def create_disparity_stats_matrix(subj,activity,src_dir,dst_dir,start_frame,duration_sec):
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
        #disparity_frame = (2*np.arctan2((disparity_frame/2),focal_length_pix))*(180/math.pi)

        disparity_mat[:,:,cnt] = disparity_frame
        cnt += 1


    #save full disparity matrix w/ 32 bit precision
    disparity_mat = disparity_mat.astype('float32', casting='same_kind')
    np.savez_compressed(dst_dir + subj + '_' + activity + '_' +'disparity.npz', disparity=disparity_mat)
    sio.savemat(dst_dir + subj + '_' + activity + '_' + 'disparity.mat', {'disparity':disparity_mat}, do_compression=True)

    #calculate gaze statistics
    #load frame status matrix
    frame_status_mat = load_frame_status_mat(src_dir, subj, activity)
    #frame_status_mat = np.load('../../3d_reprojection/data/' + subj + '/' + subj + '_' + activity + '_eye_images_transform1/' + subj + '_' + activity + '_frame_status.npy')

    frame_status_mat = frame_status_mat[np.logical_and(frame_status_mat[:,0] >= start_frame,frame_status_mat[:,0] <= end_frame)]

    print 'total frames: ' + str(len(frame_status_mat))
    print 'frames used: ' + str(mat_size)
    print 'good fixations:' + '%.2f' % (np.true_divide(len(frame_status_mat[frame_status_mat == 1]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[frame_status_mat == 1]))) + ')'
    print 'not full frame:' + '%.2f' % (np.true_divide(len(frame_status_mat[frame_status_mat == 5]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[frame_status_mat == 5]))) + ')'
    print 'saccades:' + '%.2f' % (np.true_divide(len(frame_status_mat[frame_status_mat == 3]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[frame_status_mat == 3]))) + ')'
    print 'blinks:' + '%.2f' % (np.true_divide(len(frame_status_mat[frame_status_mat == 2]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[frame_status_mat == 2]))) + ')'
    print 'missing data:' + '%.2f' % (np.true_divide(len(frame_status_mat[frame_status_mat == 4]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[frame_status_mat == 4]))) + ')'
    print 'Nans?:' + '%.2f' % (np.true_divide(len(frame_status_mat[np.isnan(frame_status_mat)]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[np.isnan(frame_status_mat)]))) + ')'

    #save parameters into text file in destination directory
    params = open(dst_dir + subj + '_' + activity + '_' + "disparity_params.txt", "w")
    params.write("\n disparity directory: " + src_dir)
    params.write("\n start frame: " + str(start_frame))
    params.write("\n end frame: " + str(end_frame))
    params.write("\n duration (sec): " + str(duration_sec))
    params.write("\n total frames: " + str(len(frame_status_mat)))
    params.write("\n frames used: " + str(mat_size))
    params.write("\n good fixations: " + '%.2f' % (np.true_divide(len(frame_status_mat[frame_status_mat == 1]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[frame_status_mat == 1]))) + ')')
    params.write("\n not full frame: " + '%.2f' % (np.true_divide(len(frame_status_mat[frame_status_mat == 5]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[frame_status_mat == 5]))) + ')')
    params.write("\n saccadess: " + '%.2f' % (np.true_divide(len(frame_status_mat[frame_status_mat == 3]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[frame_status_mat == 3]))) + ')')
    params.write("\n blinks: " + '%.2f' % (np.true_divide(len(frame_status_mat[frame_status_mat == 2]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[frame_status_mat == 2]))) + ')')
    params.write("\n missing data: " + '%.2f' % (np.true_divide(len(frame_status_mat[frame_status_mat == 4]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[frame_status_mat == 4]))) + ')')
    params.write("\n NaNs?: " + '%.2f' % (np.true_divide(len(frame_status_mat[np.isnan(frame_status_mat)]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[np.isnan(frame_status_mat)]))) + ')')
    params.close()

def update_progress_bar(update_string):
    sys.stdout.write('\r' + update_string)
    sys.stdout.flush()

def load_frame_status_mat(src_dir, subj, activity):
    frame_status_dpath = os.path.join(src_dir, '..', '..')
    frame_status_fpath = glob(os.path.join(frame_status_dpath, '*frame_status.npy'))[0]
    return np.load(frame_status_fpath)