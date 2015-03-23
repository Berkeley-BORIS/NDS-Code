import cv
import os
import cv2
import numpy
import math
import bottleneck

def make_destination_directory_tree(dst_dir,SAD_window_size,save_images,denoise,time_smooth,three_d):
	'''make a tree of directories for storing disparity map results
	'''
	
	#main directoryy
	os.mkdir(dst_dir)	

	#create destination sub directories if they don't exist	
	#for win in SAD_window_size:
	#	if not os.path.exists(dst_dir + "SAD" + str(win) + "/"):
	#		os.mkdir(dst_dir + "SAD" + str(win) + "/")	
	#	if save_images == True:
	#		if not os.path.exists(dst_dir + "image_SAD" + str(win) + "/"):
	#			os.mkdir(dst_dir + "image_SAD" + str(win) + "/")
				
	if denoise == True:
		if not os.path.exists(dst_dir + "denoised/"):
			os.mkdir(dst_dir + "denoised/")
		if save_images == True:
			if not os.path.exists(dst_dir + "image_denoised/"):
					os.mkdir(dst_dir + "image_denoised/")	
				
	if time_smooth == True:
		if not os.path.exists(dst_dir + "tsmooth_denoised/"):
			os.mkdir(dst_dir + "tsmooth_denoised/")
		if save_images == True:
			if not os.path.exists(dst_dir + "image_tsmooth_denoised/"):
				os.mkdir(dst_dir + "image_tsmooth_denoised/")
			if not os.path.exists(dst_dir + "image_tsmooth_denoised_registered_frames/"):
				os.mkdir(dst_dir + "image_tsmooth_denoised_registered_frames/")
			if not os.path.exists(dst_dir + "image_tsmooth_denoised_registered_dmaps/"):
				os.mkdir(dst_dir + "image_tsmooth_denoised_registered_dmaps/")
					
	if three_d == True:
		#if not os.path.exists(dst_dir + "three_d_tsmooth_denoised/"):
		#	os.mkdir(dst_dir + "three_d_tsmooth_denoised/")
		if not os.path.exists(dst_dir + "three_d_denoised/"):
			os.mkdir(dst_dir + "three_d_denoised/")

def store_disparity_mat_parameters(dst_dir,src_dir,min_distance_mm,max_pixel_disparity,SAD_window_size,disp_12_maxdiff,uniqueness_ratio,speckle_window_size,speckle_range,smoothing,full_DP,denoise,time_smooth):
	'''store all of the parameters used for disparity algorithm
	'''
	
	params = open(dst_dir+'disparity_algorithm_params.txt', 'w')
	params.write("\n frames source directory: " + src_dir)
	params.write("\n minimum distance in mm (estimated): " + str(min_distance_mm))
	params.write("\n minimum distance in mm (actual): " + str(32.5/math.tan(6*math.atan((max_pixel_disparity/2)/(1000*4.5)))))
	params.write("\n max pixel disparity: " + str(max_pixel_disparity))
	params.write("\n SAD window size in pixels: " + str(SAD_window_size))
	params.write("\n disparity12maxdiff: " + str(disp_12_maxdiff))
	params.write("\n uniqueness ratio: " + str(uniqueness_ratio))
	params.write("\n speckle window size: " + str(speckle_window_size))
	params.write("\n speckle range: " + str(speckle_range))
	params.write("\n smoothing: " + str(smoothing))
	params.write("\n Full DP: " + str(full_DP))
	params.write("\n denoised: " + str(denoise))
	params.write("\n temporal smoothing: " + str(time_smooth))
	params.close()	

def save_disparity_image(disparity_mat,black_offset=20):
	'''convert disparity matrix to grey scale disparity map
	'''
	
	disparity_image = disparity_mat.copy()
	#image will saturate at 100 pixel disparity
	disparity_image = (disparity_image)/100
	disparity_image[disparity_image > 1] = 1
	
	#set min disparity value to be noticably different from black
	disparity_image[disparity_image > 0] = ((255-black_offset)*disparity_image[disparity_image > 0]) + black_offset

	disparity_image = numpy.dstack((disparity_image, disparity_image, disparity_image))
	disparity_image[disparity_image[:,:,2] < 0,0] = 0 #set invalid pixels to yellow
	disparity_image[disparity_image[:,:,2] < 0,1] = 255
	disparity_image[disparity_image[:,:,2] < 0,2] = 255
	disparity_image = numpy.ceil(disparity_image).astype("uint8")
	
	return disparity_image

def temporal_median_smoothing(denoise_mat,tsmooth_disp_mat_temp,tsmooth_affine_mat_temp,tsmooth_cnt,three_d,save_images,im,dst_dir):
	'''takes three consecutive frames, applies an affine transform to closely register them, and applies a median filter to denoise
	'''
	
	#fill up data matrices with depth maps and image data
	tsmooth_disp_mat_temp[:,:,tsmooth_cnt] = denoise_mat[:,:,0].copy()
	tsmooth_affine_mat_temp[:,:,tsmooth_cnt] = cv2.imread(src_dir + im,flags=0)
	tsmooth_cnt += 1
	
	#if there are 3 frames
	if tsmooth_cnt == 3:
		
		#copy data matrices over for manipulation
		tsmooth_disp_mat = tsmooth_disp_mat_temp.copy()
		tsmooth_affine_mat = tsmooth_affine_mat_temp.copy()
	
		#estimate the affine transformation from preceeding/following frames to center frame
		afftrans_preceeding_frame = cv2.estimateRigidTransform(tsmooth_affine_mat[:,:,0].astype('uint8'), tsmooth_affine_mat[:,:,1].astype('uint8'), fullAffine=True)
		afftrans_following_frame = cv2.estimateRigidTransform(tsmooth_affine_mat[:,:,2].astype('uint8'), tsmooth_affine_mat[:,:,1].astype('uint8'), fullAffine=True)
	
		#convert disparity maps to integer matrices
		disp_preceeding_toregister = tsmooth_disp_mat[:,:,0].copy()
		disp_following_toregister = tsmooth_disp_mat[:,:,2].copy()
	
		#set invalid pixels to nans so they don't contaminate interpolated pixels
		disp_preceeding_toregister[disp_preceeding_toregister == -1] = numpy.nan
		disp_following_toregister[disp_following_toregister == -1] = numpy.nan	
	
		#initialize results
		disp_preceeding_registered = numpy.zeros((480,640))
		disp_following_registered = numpy.zeros((480,640))
	
		#apply affine transformation
		disp_preceeding_registered = cv2.warpAffine(src=disp_preceeding_toregister, M=afftrans_preceeding_frame, dsize=(640,480), dst=disp_preceeding_registered)
		disp_following_registered = cv2.warpAffine(src=disp_following_toregister, M=afftrans_following_frame, dsize=(640,480), dst=disp_following_registered)
	
		#reset invalid pixels to -1
		disp_preceeding_registered[numpy.isnan(disp_preceeding_registered)] = -1
		disp_following_registered[numpy.isnan(disp_following_registered)] = -1
		tsmooth_disp_mat[:,:,0] = disp_preceeding_registered
		tsmooth_disp_mat[:,:,2] = disp_following_registered
	
		#if there are invalid (-1) pixels in the target frame, first try to replace with with values from the preceeding frame, then from the following frame
		tsmooth_disp_mat[tsmooth_disp_mat[:,:,1] == -1,1] = tsmooth_disp_mat[tsmooth_disp_mat[:,:,1] == -1,0]
		tsmooth_disp_mat[tsmooth_disp_mat[:,:,1] == -1,1] = tsmooth_disp_mat[tsmooth_disp_mat[:,:,1] == -1,2]

		#take median of new matrix
		disp_median = bottleneck.median(tsmooth_disp_mat,axis=2)
	
		#save array
		frame_num_orig = "0"*(6-len(str(int(im.split("frame_")[1].split(".bmp")[0])))) + str(int(im.split("frame_")[1].split(".bmp")[0]))
		frame_num_new = "0"*(6-len(str(int(im.split("frame_")[1].split(".bmp")[0])-1))) + str(int(im.split("frame_")[1].split(".bmp")[0])-1)
		numpy.save(dst_dir + "tsmooth_denoised/disp_tsmooth_denoised" + im.replace(frame_num_orig,frame_num_new).replace("bmp","npy"), disp_median)
	
		if three_d == True:
			depth_mat = numpy.zeros((disp_median.shape[0],disp_median.shape[1],3))
			depth_mat = cv2.reprojectImageTo3D(disparity=disp_median.astype('float32'), Q=Q, _3dImage=depth_mat, handleMissingValues=False, ddepth=-1)
			#set invalid pixels to (0,0,0)
			depth_mat[disp_median == -1] = 0
			numpy.save(dst_dir + "three_d_tsmooth_denoised/depth_tsmooth_denoised" + im.replace(frame_num_orig,frame_num_new).replace("bmp","npy"), depth_mat)
		
		if save_images == True:
			#convert to image format for easy visualization
			disp_median_image = save_disparity_image(diso_median,black_offset=20)
			cv.SaveImage(dst_dir + "image_tsmooth_denoised/disp_smooth_denoised" + im.replace(frame_num_orig,frame_num_new).replace("bmp","png"), cv.fromarray(disp_median_image.copy()))
					
		#set tsmooth arrays up for next frame
		tsmooth_disp_mat_temp[:,:,0] = tsmooth_disp_mat_temp[:,:,1]
		tsmooth_disp_mat_temp[:,:,1] = tsmooth_disp_mat_temp[:,:,2]
		tsmooth_disp_mat_temp[:,:,2] = numpy.zeros((480,640))
		tsmooth_affine_mat_temp[:,:,0] = tsmooth_affine_mat_temp[:,:,1]
		tsmooth_affine_mat_temp[:,:,1] = tsmooth_affine_mat_temp[:,:,2]
		tsmooth_affine_mat_temp[:,:,2] = numpy.zeros((480,640))
		tsmooth_cnt = 2
	
	return tsmooth_cnt
	
def mask_frames(src_dir,dst_dir,mask_pixel_width=40):
	''' mask messy edges from disparity map frames
	example call: disparity_map_utils.mask_frames('../data/disparity_02_15_12_sgbm_smooth/','../data/disparity_02_15_12_sgbm_smooth_masked/')
	'''
	
	#create destination directory if it doesn't exist
	if not os.path.exists(dst_dir):
		os.mkdir(dst_dir)

	#create masking bars
	bar_vert = numpy.zeros( (480-(2*mask_pixel_width),mask_pixel_width) )
	bar_horz = numpy.zeros( (mask_pixel_width, 640) )
	#iterate through disparity maps in source directory
	imlist = os.listdir(src_dir)
	for im in imlist:
		#if it's an image and not a params file
		if im.find(".png") > 0:
			#import pdb; pdb.set_trace()
			im_orig = cv2.imread(src_dir + im,flags=0)

			#apply mask
			im_masked = im_orig[mask_pixel_width:im_orig.shape[0]-mask_pixel_width,mask_pixel_width:im_orig.shape[1]-mask_pixel_width]
			im_masked = numpy.append(bar_vert,im_masked,axis=1); im_masked = numpy.append(im_masked,bar_vert,axis=1)
			im_masked = numpy.append(bar_horz,im_masked,axis=0); im_masked = numpy.append(im_masked,bar_horz,axis=0)
			cv.SaveImage(dst_dir + "masked_" + im, cv.fromarray(im_masked.copy()))
		
	#save parameters into text file in destination directory
	params = open(dst_dir+'frame_mask_params.txt', 'w')
	params.write("\n disparity maps source directory: " + src_dir)
	params.write("\n mask pixel width: " + str(mask_pixel_width))
	params.close()