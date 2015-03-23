import numpy as np
import os
import math
import bottleneck as bn
import  matplotlib.pyplot as plt
import matplotlib.cm as cm
import cv
import scipy.io as sio
import sys

def apply_circle_mask(img,circle_center):
    '''make 8 deg circle mask for eye images and disparity maps
    '''

    radius = 81. #10deg pixel radius with 583 fl

    start_pixel = circle_center-radius
    end_pixel = circle_center+radius+1

    #crop image to radius
    if len(img.shape) == 3:
        img = img[start_pixel:end_pixel,start_pixel:end_pixel,:]
    else:
        img = img[start_pixel:end_pixel,start_pixel:end_pixel]

    x,y = np.meshgrid(range(0,(2*int(radius))+1),range(0,(2*int(radius))+1))
    mask = (  (((x-radius)/(radius))**2 + ((y-radius)/(radius))**2) < 1 )
    img[mask == False] = np.nan

    return img


def create_depth_stats_matrix(subj,activity,src_dir,dst_dir,start_frame,duration_sec):
    '''take in disparity values from frames in a directory and calculate disparity distribution
    '''

    #create destination directory if it doesn't exist
    if not os.path.exists(dst_dir):
        os.mkdir(dst_dir)

    #convert duration in milliseconds to a end frame
    num_frames = int((duration_sec*1000)/33.333)
    end_frame = start_frame + num_frames

    #iterate through source directory of depth maps and grab desired frames
    imlist = os.listdir(src_dir)
    frame_num_list = np.ones((len(imlist),1))*np.nan
    cnt = 0
    for im in imlist:
        frame_num_list[cnt] = int(im.split('_')[2].strip('.npy'))
        cnt += 1

    if start_frame > 0:
        imlist = np.asarray(imlist)[np.where(np.logical_and(frame_num_list >= start_frame, frame_num_list <=  end_frame))[0]].tolist()

    mat_size = len(imlist)
    print 'number of frames = ' + str(mat_size)

    if mat_size > 5000:
        print 'mat is toooo big!'
        sys.exit(1)

    #allocate data matrix for disparity frames
    depth_mat = np.ones((163,163,mat_size))*np.nan

    #start counter
    cnt = 0
    frames_used = []
    #load in each disparity map for averaging
    for im in imlist:
        # add this frame number to the frames_used list
        frame_number = int(im.strip('.npy').split('_')[-1])
        frames_used.append(frame_number)
        #load disparity map
        depth_frame = np.load(src_dir + im)

        depth_frame = apply_circle_mask(depth_frame,275)

        #if not np.isnan(np.nanmax(disparity_frame)):
        current_im_string = "%d of %d (%s)" % (cnt, mat_size, im)
        update_progress_bar(current_im_string)
        #convert disparities from pixels to degrees
        #disparity_frame = (2*np.arctan2((disparity_frame/2),focal_length_pix))*(180/math.pi)

        depth_mat[:,:,cnt] = depth_frame
        cnt += 1


    #save full depth matrix

    #depth_mat = np.round(depth_mat)
    #depth_mat = np.int8(depth_mat)

    # save depth mat w/ 32 bit precision
    depth_mat = depth_mat.astype('float32', casting='same_kind')
    np.savez_compressed(dst_dir + subj + '_' + activity + '_' +'depth.npz',depth=depth_mat)
    sio.savemat(dst_dir + subj + '_' + activity + '_' + 'depth.mat', {'depth':depth_mat}, do_compression=True)
    np.save(os.path.join(dst_dir, subj+'_'+activity+'_frames_used.npy'), np.array(frames_used))

def update_progress_bar(update_string):
    sys.stdout.write('\r' + update_string)
    sys.stdout.flush()