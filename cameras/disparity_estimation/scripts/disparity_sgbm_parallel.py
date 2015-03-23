import cv
import os
import fnmatch
import cv2
import numpy
import math
import bottleneck
import string
import pp
import disparity_sgbm_parallel_utils


def disparity_sgbm(min_distance_mm,SAD_window_size,disp_12_maxdiff,uniqueness_ratio,speckle_filter,smoothing,speckle_window_size,speckle_range,prefilter_cap,full_DP,save_images,denoise,time_smooth,three_d,subj,parallel_input):
	''' perform SGBM disparity estimation
	
	optionally: denoise = use larger SAD window sizes to filter out false matches from smallest window size and fill in missing pixels
				time_smooth = register and combine a window of 3 frames in temporal sequence to denoise depth map 
				
	http://code.google.com/p/tjpstereovision/source/browse/code/Python/stereorectify.py?spec=svnba9b04d5082c96c46860edca01e9133c93c53ead&r=ba9b04d5082c96c46860edca01e9133c93c53ead
	example call: disparity_map_utils.disparity_sgbm('../data/frames_rect/','../data/disparity_sgbm/')
	'''

	src_dir = "../../image_rectification/data/" + subj + "/" + subj + "_" + parallel_input[0] + "_frames_rect/"
	dst_dir = "../data/" + subj + "/" + subj + "_" + parallel_input[0] + "_disparity_maps/"
	frames = parallel_input[1]

	#grab calbration frames directory
	cam_dir = "../../stereo_calibration/data/"  + subj + "/"
	for dir in os.listdir(cam_dir):
		if fnmatch.fnmatch(dir,'calibration_frames*'):
			cam_dir = cam_dir + dir + '/'
			break
	
	#create destination directory tree if it doesn't exist
	if not os.path.exists(string.join(string.split(dst_dir,'/')[0:-2],'/')):
		os.mkdir(string.join(string.split(dst_dir,'/')[0:-2],'/'))
	if not os.path.exists(dst_dir): 
		disparity_sgbm_parallel_utils.make_destination_directory_tree(dst_dir,SAD_window_size,save_images,denoise,time_smooth,three_d)
		
	#load perspective transformation matrix from stereo rectification, for 3d reconstruction
	Q = cv.Load(cam_dir + "Disp2depth_matrix.xml")
	Q = numpy.asarray(Q)
												
	#turn speckle filtering off if given flag
	if speckle_filter == False:
		speckle_window_size = 0
		speckle_range = 0
		
	#calculate max pixel disparity to look for based on minimum distance we want to measure
	half_disparity_angle = math.atan(32.5/min_distance_mm)
	half_disparity_pixels = (1000*4.5*math.tan(half_disparity_angle))/6
	max_pixel_disparity   = int(2*half_disparity_pixels)
	#make this value a multiple of 16
	max_pixel_disparity = 16*((max_pixel_disparity/16)+1)
	
	#save parameters into text file in destination directory
	disparity_sgbm_parallel_utils.store_disparity_mat_parameters(dst_dir,src_dir,min_distance_mm,max_pixel_disparity,SAD_window_size,disp_12_maxdiff,uniqueness_ratio,speckle_window_size,speckle_range,smoothing,full_DP,denoise,time_smooth)
	
	#iterate through images in source directory and crop imlist if number of frames is given
	imlist = numpy.array(fnmatch.filter(os.listdir(src_dir), '*.bmp'))
	#if len(frames) == 2 and not frames[0] == 0 and not frames[1] == 0:
	#	frames = range(frames[0], frames[1]+1)  # list of frames to calculate disparity maps on INCLUSIVE of last frame number
	#imlist = imlist[frames]

	if len(frames) == 2:
	 	if not frames[1] == 0:
			frames = range(frames[0], frames[1]+1)  # list of frames to calculate disparity maps on INCLUSIVE of last frame number
		else:
			imlist = imlist
	else:
		imlist = imlist[frames]
	
	#initialize time smoothing matrix
	tsmooth_disp_mat_temp = numpy.zeros((480,640,3))
	tsmooth_affine_mat_temp = numpy.zeros((480,640,3))
	tsmooth_cnt = 0
	
	#make disparity maps
	for im in imlist:
		
		#only perform matching once per frame, so only initiate if it is a left frame
		if im.find('cam1') > 0:
			#print im

			#grab left and right frames
			im_left = cv2.imread(src_dir + im,flags=0)
			im_right = cv2.imread(src_dir + im.replace("cam1","cam2"),flags=0)
			
			
			#initialize denoising matrix
			denoise_mat = numpy.zeros((im_left.shape[0],im_left.shape[1],len(SAD_window_size)))
			win_cnt = 0
			
			#hack: add black bar to the left side of the stereo images so the disparity map doesn't get cropped
			bar = numpy.zeros( (480,max_pixel_disparity) , dtype = 'uint8')
			im_left = numpy.append(bar,im_left,axis=1)
			im_right = numpy.append(bar,im_right,axis=1)
			
			#estimate disparities using each SAD window size
			for win in SAD_window_size:
					
				#set smoothing parameters
				if smoothing == False:
					P1 = 0
					P2 = 0
				elif smoothing == True:
					P1 = 8*win*win
					P2 = 32*win*win
				
				stereo = cv2.StereoSGBM(minDisparity=0, 
					numDisparities=max_pixel_disparity, 
					SADWindowSize=win, 
					P1=P1, 
					P2=P2, 
					disp12MaxDiff=disp_12_maxdiff, 
					preFilterCap=prefilter_cap,
					uniquenessRatio=uniqueness_ratio, 
					speckleWindowSize=speckle_window_size, 
					speckleRange=speckle_range,
					fullDP=full_DP)
			
				#disparity map
				disparity_left = stereo.compute(im_left, im_right).astype(numpy.float32)
			
				#remove hacky bar
				disparity_left = disparity_left[:,max_pixel_disparity:]
			
				#divide disparity values by 16 to get pixel disparity, pixel value of -1 mean invalid
				disparity_left = disparity_left/16
				
				#save depth map as numpy array
				
				#no longering saving individual SAD disparity maps
				#numpy.save(dst_dir + "SAD" + str(win) + "/disp_SAD" + str(win) + im.replace("bmp","npy"), disparity_left)
				
				#add to denoising matrix
				if denoise == True:
					denoise_mat[:,:,win_cnt] = disparity_left
					win_cnt += 1
					
				#if save_images == True:
					#convert to image format for easy visualization
					
					#no longer saving individual SAD images
					#disparity_left_image = disparity_sgbm_parallel_utils.save_disparity_image(disparity_left,black_offset=20)
					#cv.SaveImage(dst_dir + "image_SAD" + str(win) + "/disp_SAD" + str(win) + im.replace("bmp","png"), cv.fromarray(disparity_left_image.copy()))
					
							
			#combine frames to get denoised estimate
			if denoise == True:
				
				#false matches in the fine scale become zeros in the coarse scale, if fine scale pixel has value, but coarsest scale value is zero, set pixel to zero
				denoise_mat[numpy.logical_and(denoise_mat[:,:,0] > -1,denoise_mat[:,:,2] == 0),0] = denoise_mat[numpy.logical_and(denoise_mat[:,:,0] > -1,denoise_mat[:,:,2] == 0),2]

				#if fine scale pixel is invalid, set pixel to next coarse scale value, repeat
				denoise_mat[denoise_mat[:,:,0] == -1,0] = denoise_mat[denoise_mat[:,:,0] == -1,1]
				denoise_mat[denoise_mat[:,:,0] == -1,0] = denoise_mat[denoise_mat[:,:,0] == -1,2]
				
				#save depth map as numpy array
				numpy.save(dst_dir + "denoised/disp_denoised" + im.replace("bmp","npy"), denoise_mat[:,:,0])
				
				if three_d == True:
					depth_mat_denoise = numpy.zeros((denoise_mat[:,:,0].shape[0],denoise_mat[:,:,0].shape[1],3))
					depth_mat_denoise = cv2.reprojectImageTo3D(disparity=denoise_mat[:,:,0].astype('float32'), Q=Q, _3dImage=depth_mat_denoise, handleMissingValues=False, ddepth=-1)
					#set invalid pixels to (0,0,0)
					depth_mat_denoise[denoise_mat[:,:,0] == -1] = 0
					numpy.save(dst_dir + "three_d_denoised/depth_denoised" + im.replace("bmp","npy"), depth_mat_denoise)
				
				if save_images == True:
					#convert to image format for easy visualization
					denoise_mat = denoise_mat[:,:,0].copy()
					denoise_mat_image = disparity_sgbm_parallel_utils.save_disparity_image(denoise_mat,black_offset=20)
					cv.SaveImage(dst_dir + "image_denoised/disp_denoised" + im.replace("bmp","png"), cv.fromarray(denoise_mat_image.copy()))
			
			#use preceeding and following frames to fill in missing pixels and smooth each frame
			if time_smooth == True:
				tsmooth_cnt = disparity_sgbm_parallel_utils.temporal_median_smoothing(denoise_mat,tsmooth_disp_mat_temp,tsmooth_affine_mat_temp,tsmooth_cnt,three_d,save_images,im,dst_dir)	
					

	