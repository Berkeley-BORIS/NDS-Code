import cv
import os.path
import cv2
import numpy
import scipy.io
import math
import os
import string
import camera_to_eye_utils_parallel
import pp

def camera_image_to_eye_image(subj,transform_src_file,fixation_src_file,ipd_cm,frames,targets,drift,parallel_input, fixed_fixations=''):
	'''transforms left camera image into left and right eye images
	'''
	
	fixation_src_mat = parallel_input

        if fixed_fixations and not fixed_fixations.startswith('_'):
            fixed_fixations = '_'+fixed_fixations
	
	if drift == 0:
		dst_dir = '../data/' + subj + fixed_fixations + '/' + subj + '_' + fixation_src_mat + '_eye_images_transform' + transform_src_file.split('_')[3].strip('.npz') + '/'
	elif drift == 1:
		dst_dir = '../data/' + subj + fixed_fixations + 'drift/' + subj + 'drift_' + fixation_src_mat + '_eye_images_transform' + transform_src_file.split('_')[3].strip('.npz') + '/'

	
	# Normal source directory operation
	depth_src_dir = '../../disparity_estimation/data/' + subj + '/' + subj + '_' + fixation_src_mat + '_disparity_maps/three_d_denoised/'
	image_src_dir = '../../image_rectification/data/' + subj + '/' + subj + '_' + fixation_src_mat + '_frames_rect/'
	
	#create destination directory if it doesn't exist
	if not os.path.exists(string.join(string.split(dst_dir,'/')[0:-2],'/')):
		os.mkdir(string.join(string.split(dst_dir,'/')[0:-2],'/'))
	if not os.path.exists(dst_dir):
		camera_to_eye_utils_parallel.make_destination_directory_tree(dst_dir)
	if not os.path.exists(dst_dir+'10deg/'):
		camera_to_eye_utils_parallel.make_destination_directory_tree(dst_dir+'10deg/')
	#if not os.path.exists(dst_dir+'10degfull/'):
	#	camera_to_eye_utils_parallel.make_destination_directory_tree(dst_dir+'10degfull/')		

	#save parameters into text file in destination directory
	camera_to_eye_utils_parallel.save_reprojection_parameters(dst_dir,transform_src_file,depth_src_dir,image_src_dir,fixation_src_file)
		
	#load in fixation location at each frame
	fixation_mat = scipy.io.loadmat(fixation_src_file)
	
	if fixation_src_mat.find('_orig') > -1 or fixation_src_mat.find('_cont') > -1:
		fixation_src_mat = fixation_src_mat[0:-5]
		
	if targets == 0:
		
		if fixation_src_mat == 'task_walk_1':
			fixation_info = fixation_mat['walk_1']
		elif fixation_src_mat == 'task_walk_2':
			fixation_info = fixation_mat['walk_2']
		elif fixation_src_mat == 'task_walk_3':
			fixation_info = fixation_mat['walk_3']	
		else:
			fixation_info = fixation_mat['task']

		#fixation_info = fixation_mat['walk_1']
	else:
		fixation_info = fixation_mat[fixation_src_mat]
		
	fixation_coords_le = fixation_info[:,(4,5,6)]
	fixation_coords_re = fixation_info[:,(7,8,9)]
	#and the ground truth target locations in cyclopean coordinates
	target_coords = fixation_info[:,(10,11,12)]
	target_tilt_slant = fixation_info[:,(13,14)]
	#and flags
	eye_flag      = fixation_info[:,(15)]
	#eye flags: 1 = found fixation, 2 = blink, 3 = saccade, 4 = missing data, to add 5 = not full frame

	#load in pixel locations of targets in camera frames if present
	if os.path.isfile(image_src_dir + 'E_coords.pydict'):
		target_camera_pixels = numpy.load(image_src_dir + 'E_coords.pydict')
	
	#load in transformation from cyclopean-eye-to-camera
	#and calculate inverse rotation and translation for camera-to-cyclopean-eye
	R_eye_to_cam,R_cam_to_eye,tvec_eye_to_cam,tvec_cam_to_eye,intrinsic_mat = camera_to_eye_utils_parallel.load_and_invert_transforms(transform_src_file)
	
	#left camera intrinsics
	focal_length_y = intrinsic_mat[1][1]
	focal_length_x = intrinsic_mat[0][0]
	center_y = intrinsic_mat[1][2]
	center_x = intrinsic_mat[0][2]

	#eye intrinsics (550x550 pixel images cover 25x25deg field of view), 3.5mm fl 583 pixel units
	intrinsic_eye = numpy.zeros( (3,3) )
	intrinsic_eye[0,0] = 583
	intrinsic_eye[1,1] = 583
	intrinsic_eye[0,2] = 275
	intrinsic_eye[1,2] = 275
	intrinsic_eye[2,2] = 1.0
	dist_coeffs_eye = numpy.zeros( (5,1) )

	#iterate through source directory of 3d coordinates and grab desired frames
	imlist = os.listdir(depth_src_dir)
	if frames[0] > 0 or frames[1] > 0:
		imlist = imlist[frames[0]:frames[1]]
	
	#initialize matrix with flags to identify status of each frame	
	frame_status_mat = numpy.zeros((len(imlist),2)) #like eye flags, but for non-target trials
	frame_cnt = 0
					
	#allocate error matrix if using targets
	if targets == 1:
		error_mat = numpy.ones((len(imlist),21))*numpy.nan
		
	#for each frame matrix of 3d camera coordinates, calculate eye images and disparity maps	
	for im in imlist:
		
		#make sure necessary files exist
		if os.path.isfile(depth_src_dir + im) and os.path.isfile(image_src_dir + 'rect_cam1_frame_' + im.split('_frame_')[1].strip('.npy') + '.bmp'):
			
			#frame number
			frame = int(im.split('_frame_')[1].strip('.npy'))
			#print str(frame)
			
			#frame status
			frame_status_mat[frame_cnt,0] = frame
			frame_status_mat[frame_cnt,1] = eye_flag[frame]
			
			if targets == 1:
				#add target pts and store error info if present
				error_mat[frame_cnt,0] = frame #frame number
				error_mat[frame_cnt,(1,2,3)] = target_coords[frame,:] #target cyclopean x,y,z
				error_mat[frame_cnt,4] = target_tilt_slant[frame,0] #tilt
				error_mat[frame_cnt,5] = target_tilt_slant[frame,1] #slant	
				error_mat[frame_cnt,6] = int(fixation_src_mat.split('_')[2]) #which distance repeat		
			
			#check the eye flag for this frame. If this frame is not a fixation or saccade ( i.e. a blink or missing data) skip reprojection
			if frame_status_mat[frame_cnt,1] == 2 or frame_status_mat[frame_cnt,1] == 4:
				if frame_status_mat[frame_cnt,1] == 2:
					#blink: images are all black
					lefteye_image,righteye_image,disparity_image,vdisparity_image,disparity_mat,vdisparity_mat,depth_mat = camera_to_eye_utils_parallel.initialize_images_and_mats()
				elif frame_status_mat[frame_cnt,1] == 4:
					#missing data: images are all grey
					lefteye_image,righteye_image,disparity_image,vdisparity_image,disparity_mat,vdisparity_mat,depth_mat = camera_to_eye_utils_parallel.initialize_images_and_mats()
					lefteye_image = lefteye_image+127
					righteye_image = righteye_image+127	
					disparity_image = disparity_image+127	
					vdisparity_image = vdisparity_image+127	
				
				#save files
				#print "No fixation, Saving files...\n"
				cv.SaveImage(dst_dir + 'lefteye/lefteye_frame_' + im.split('_frame_')[1].strip('.npy') + '.png', cv.fromarray(lefteye_image.copy()))
				cv.SaveImage(dst_dir + 'righteye/righteye_frame_' + im.split('_frame_')[1].strip('.npy') + '.png', cv.fromarray(righteye_image.copy()))
				cv.SaveImage(dst_dir + 'disparityeye_image/disparityeye_image_frame_' + im.split('_frame_')[1].strip('.npy') + '.png', cv.fromarray(disparity_image.copy()))
				numpy.save(dst_dir + 'disparityeye/disparityeye_frame_' + im.split('_frame_')[1], disparity_mat)
				cv.SaveImage(dst_dir + 'vdisparityeye_image/vdisparityeye_image_frame_' + im.split('_frame_')[1].strip('.npy') + '.png', cv.fromarray(vdisparity_image.copy()))
				numpy.save(dst_dir + 'vdisparityeye/vdisparityeye_frame_' + im.split('_frame_')[1], vdisparity_mat)	
				#numpy.save(dst_dir + 'deptheye/deptheye_frame_' + im.split('_frame_')[1], depth_mat)			
		
				#increment frame cnt
				frame_cnt += 1
				
				continue
					
			#load camera 3d coordinates
			camera_coords = numpy.load(depth_src_dir + im)
			camera_coords = camera_coords.astype('float64')
			camera_coords_orig = camera_coords.copy()
		
			#remove noisy edges
			camera_coords = camera_to_eye_utils_parallel.remove_image_edges(camera_coords,25)
		
			if targets == 1:
				#grab ground truth target location in cyclopean coordinates 
				target_point = camera_to_eye_utils_parallel.ground_truth_target_cyclopean_coords(target_coords,frame)
				#grab camera pixel target location in camera coordinates
				target_camera_coords = camera_to_eye_utils_parallel.camera_target_camera_coords(target_point,target_camera_pixels,camera_coords,center_x,center_y,focal_length_x,focal_length_y,frame)

			#for points with inf distance, convert to just very very large distance
			inf_coords = numpy.zeros((numpy.where(camera_coords[:,:,2] == numpy.inf)[0].size,3))
			inf_coords[:,0] = -1e100 * (-(center_x-numpy.where(camera_coords[:,:,2] == numpy.inf)[1])/focal_length_x)
			inf_coords[:,1] = -1e100 * (-(center_y-numpy.where(camera_coords[:,:,2] == numpy.inf)[0])/focal_length_y)
			inf_coords[:,2] = -1e100*numpy.ones((1,numpy.where(camera_coords[:,:,2] == numpy.inf)[0].size))	
			camera_coords[numpy.where(camera_coords[:,:,2] == numpy.inf)] = inf_coords
				
			#for points at 0,0,0 placeholder replace with nans
			camera_coords[numpy.where(camera_coords[:,:,2] == 0)] = numpy.nan

			#reshape camera coords into 2d array with each column an x,y,z point
			camera_mat = numpy.concatenate([numpy.reshape(camera_coords[:,:,0],(1,307200)),numpy.reshape(camera_coords[:,:,1],(1,307200)),numpy.reshape(camera_coords[:,:,2],(1,307200))])	
					
			#grab fixation pt in cyclopean coords for each eye
			fixation_point_le = fixation_coords_le[frame,:]
			fixation_point_re = fixation_coords_re[frame,:]
			
			#print frame
			#print fixation_point_le
			#print fixation_point_re
		
			#print 'Reprojecting points...'
			#calculate R and t for each eye to fixation pt
			R_left,t_left,R_right,t_right = camera_to_eye_utils_parallel.find_eyes_Rt(numpy.matrix(numpy.append(numpy.reshape(fixation_point_le,(3,1)),numpy.reshape(fixation_point_re,(3,1)),axis=1)),ipd_cm)
			helmholtz_coords = camera_to_eye_utils_parallel.fixation_to_Helmholtz(numpy.append(numpy.reshape(fixation_point_le,(3,1)),numpy.reshape(fixation_point_re,(3,1)),axis=1), ipd_cm, G=0.8)
			
			R_left = camera_to_eye_utils_parallel.R_from_Helmholtz(phi=helmholtz_coords[0,0], theta=helmholtz_coords[1,0], psi=helmholtz_coords[2,0])
			R_right = camera_to_eye_utils_parallel.R_from_Helmholtz(phi=helmholtz_coords[0,1], theta=helmholtz_coords[1,1], psi=helmholtz_coords[2,1])
					
			#cast as arrays
			R_left = numpy.asarray(R_left)
			R_right = numpy.asarray(R_right)
			t_left = numpy.asarray(t_left)
			t_right = numpy.asarray(t_right)
	
			#transform camera coordinates into cyclopean coordinates, then into left and right eye coordinates
			lefteye_3d_coords_cyclo = R_cam_to_eye.dot(camera_mat + numpy.reshape(tvec_cam_to_eye,(3,1)))
			righteye_3d_coords_cyclo = R_cam_to_eye.dot(camera_mat + numpy.reshape(tvec_cam_to_eye,(3,1)))
			lefteye_3d_coords = R_left.T.dot(lefteye_3d_coords_cyclo + t_left)
			righteye_3d_coords = R_right.T.dot(righteye_3d_coords_cyclo + t_right)
	
			#store cyclopean radial distance
			cyclopean_radial_distance = camera_to_eye_utils_parallel.project_to_pixels(intrinsic_eye,lefteye_3d_coords_cyclo,calc_radial_distance=1)
			
			#project 3d points to pixels
			righteye_pixels = camera_to_eye_utils_parallel.project_to_pixels(intrinsic_eye,righteye_3d_coords)
			lefteye_pixels = camera_to_eye_utils_parallel.project_to_pixels(intrinsic_eye,lefteye_3d_coords)
		
			#get coords of fixation pt
			lefteye_fixation = R_left.T.dot(numpy.reshape(fixation_point_le,(3,1)) + t_left)
			righteye_fixation = R_right.T.dot(numpy.reshape(fixation_point_re,(3,1)) + t_right)
			#project fixation pt to pixels
			lefteye_fixation_pixels = camera_to_eye_utils_parallel.project_to_pixels(intrinsic_eye,lefteye_fixation)
			righteye_fixation_pixels = camera_to_eye_utils_parallel.project_to_pixels(intrinsic_eye,righteye_fixation)
			
			if targets == 1:
				#get coords of GT target
				lefteye_target = R_left.T.dot(numpy.reshape(target_point,(3,1)) + t_left)
				righteye_target = R_right.T.dot(numpy.reshape(target_point,(3,1)) + t_right)
				#project GT target to pixels
				lefteye_target_pixels = camera_to_eye_utils_parallel.project_to_pixels(intrinsic_eye,lefteye_target)
				righteye_target_pixels = camera_to_eye_utils_parallel.project_to_pixels(intrinsic_eye,righteye_target)
		
				#get coords of camera target pt
				#first transform from camera coordinates into cyclopean coordinates, then into left and right eye coordinates
				lefteye_target_camera_cyclo = R_cam_to_eye.dot(numpy.reshape(target_camera_coords,(3,1)) + numpy.reshape(tvec_cam_to_eye,(3,1)))
				righteye_target_camera_cyclo = R_cam_to_eye.dot(numpy.reshape(target_camera_coords,(3,1)) + numpy.reshape(tvec_cam_to_eye,(3,1)))
				lefteye_target_camera = R_left.T.dot(lefteye_target_camera_cyclo + t_left)
				righteye_target_camera = R_right.T.dot(righteye_target_camera_cyclo + t_right)
				#project camera target to pixels
				
			#	print 'left'
			#	print tvec_cam_to_eye
				#print 'right'
				#print righteye_target_camera_cyclo
				
				lefteye_camera_target_pixels = camera_to_eye_utils_parallel.project_to_pixels(intrinsic_eye,lefteye_target_camera)
				righteye_camera_target_pixels = camera_to_eye_utils_parallel.project_to_pixels(intrinsic_eye,righteye_target_camera)
		
			#create new eye images
			lefteye_image,righteye_image,disparity_image,vdisparity_image,disparity_mat,vdisparity_mat,depth_mat = camera_to_eye_utils_parallel.initialize_images_and_mats()

			#load in camera image
			cam_image = cv2.imread(image_src_dir + 'rect_cam1_frame_' + im.split('_frame_')[1].strip('.npy') + '.bmp', flags=0)
	
			#print "Filling eye images..."
			
			#for each pixel in the camera image, relocate it to the eye images
			lefteye_image = camera_to_eye_utils_parallel.create_eye_image(lefteye_pixels,cam_image,lefteye_image)
			righteye_image = camera_to_eye_utils_parallel.create_eye_image(righteye_pixels,cam_image,righteye_image)
			
			disparity_mat,vdisparity_mat = camera_to_eye_utils_parallel.create_disparity_mats(lefteye_pixels,righteye_pixels,disparity_mat,vdisparity_mat)
				
			#for each radial distance in the camera image, relocate it to the left eye
			depth_mat = camera_to_eye_utils_parallel.create_eye_image(lefteye_pixels,cyclopean_radial_distance,depth_mat,distance = 1)
			
			
			
			#make disparity image
			disparity_image = camera_to_eye_utils_parallel.make_colormap_disparity(disparity_image,disparity_mat)
			vdisparity_image = camera_to_eye_utils_parallel.make_colormap_disparity(vdisparity_image,vdisparity_mat)
		
			#add yellow fixation pt
			righteye_image = camera_to_eye_utils_parallel.draw_fixation_point(numpy.dstack((righteye_image,righteye_image,righteye_image)),righteye_fixation_pixels)
			lefteye_image = camera_to_eye_utils_parallel.draw_fixation_point(numpy.dstack((lefteye_image,lefteye_image,lefteye_image)),lefteye_fixation_pixels)
			disparity_image = camera_to_eye_utils_parallel.draw_fixation_point(disparity_image,lefteye_fixation_pixels)
			vdisparity_image = camera_to_eye_utils_parallel.draw_fixation_point(vdisparity_image,lefteye_fixation_pixels)
		
			if targets == 1:
				#add target pts
				
				#if numpy.isnan(target_point[2]) == 0:
				#	#cyan is projection of ground truth target to eyes
				#	righteye_image = draw_fixation_point(righteye_image,righteye_target_pixels,rgb=(0,255,255))
				#	lefteye_image = draw_fixation_point(lefteye_image,lefteye_target_pixels,rgb=(0,255,255))
			
				if numpy.isnan(target_camera_coords[2]) == 0:
					#magenta is projection of target pixels in camera to eyes
					righteye_image = camera_to_eye_utils_parallel.draw_fixation_point(righteye_image,righteye_camera_target_pixels,rgb=(255,0,255))
					lefteye_image = camera_to_eye_utils_parallel.draw_fixation_point(lefteye_image,lefteye_camera_target_pixels,rgb=(255,0,255))
		


		#save files
		#print "Saving files...\n"
		cv.SaveImage(dst_dir + 'lefteye/lefteye_frame_' + im.split('_frame_')[1].strip('.npy') + '.png', cv.fromarray(lefteye_image.copy()))
		cv.SaveImage(dst_dir + 'righteye/righteye_frame_' + im.split('_frame_')[1].strip('.npy') + '.png', cv.fromarray(righteye_image.copy()))
		cv.SaveImage(dst_dir + 'disparityeye_image/disparityeye_image_frame_' + im.split('_frame_')[1].strip('.npy') + '.png', cv.fromarray(disparity_image.copy()))
		numpy.save(dst_dir + 'disparityeye/disparityeye_frame_' + im.split('_frame_')[1], disparity_mat)
		cv.SaveImage(dst_dir + 'vdisparityeye_image/vdisparityeye_image_frame_' + im.split('_frame_')[1].strip('.npy') + '.png', cv.fromarray(vdisparity_image.copy()))
		numpy.save(dst_dir + 'vdisparityeye/vdisparityeye_frame_' + im.split('_frame_')[1], vdisparity_mat)
		numpy.save(dst_dir + 'deptheye/deptheye_frame_' + im.split('_frame_')[1], depth_mat)
	
		#save truncated files (10deg circular window for analysis)
		lefteye_image_10deg = camera_to_eye_utils_parallel.apply_circle_mask(lefteye_image,intrinsic_eye[0,2])
		righteye_image_10deg = camera_to_eye_utils_parallel.apply_circle_mask(righteye_image,intrinsic_eye[0,2])
		#disparity_image_10deg = camera_to_eye_utils_parallel.apply_circle_mask(disparity_image,intrinsic_eye[0,2])
		#vdisparity_image_10deg = camera_to_eye_utils_parallel.apply_circle_mask(vdisparity_image,intrinsic_eye[0,2])
		disparity_mat_10deg = camera_to_eye_utils_parallel.apply_circle_mask(disparity_mat,intrinsic_eye[0,2])
		vdisparity_mat_10deg = camera_to_eye_utils_parallel.apply_circle_mask(vdisparity_mat,intrinsic_eye[0,2])
		depth_mat_10deg = camera_to_eye_utils_parallel.apply_circle_mask(depth_mat,intrinsic_eye[0,2])
		
		cv.SaveImage(dst_dir + '10deg/' + 'lefteye/lefteye_frame_' + im.split('_frame_')[1].strip('.npy') + '.png', cv.fromarray(lefteye_image_10deg.copy()))
		cv.SaveImage(dst_dir + '10deg/' + 'righteye/righteye_frame_' + im.split('_frame_')[1].strip('.npy') + '.png', cv.fromarray(righteye_image_10deg.copy()))
		#cv.SaveImage(dst_dir + '10deg/' + 'disparityeye_image/disparityeye_image_frame_' + im.split('_frame_')[1].strip('.npy') + '.png', cv.fromarray(disparity_image_10deg.copy()))
		numpy.save(dst_dir + '10deg/' + 'disparityeye/disparityeye_frame_' + im.split('_frame_')[1], disparity_mat_10deg)
		#cv.SaveImage(dst_dir + '10deg/' + 'vdisparityeye_image/vdisparityeye_image_frame_' + im.split('_frame_')[1].strip('.npy') + '.png', cv.fromarray(vdisparity_image_10deg.copy()))
		numpy.save(dst_dir + '10deg/' + 'vdisparityeye/vdisparityeye_frame_' + im.split('_frame_')[1], vdisparity_mat_10deg)
		numpy.save(dst_dir + '10deg/' + 'deptheye/deptheye_frame_' + im.split('_frame_')[1], depth_mat_10deg)
		
		#check for missing edges
		#if any edges contain only nans, the data is missing, add flag
		if numpy.any(numpy.isnan(disparity_mat_10deg[0:207,1])==0) == False or  numpy.any(numpy.isnan(disparity_mat_10deg[0:207,-2])==0) == False or numpy.any(numpy.isnan(disparity_mat_10deg[1,0:207])==0) == False or numpy.any(numpy.isnan(disparity_mat_10deg[-2,0:207])==0) == False:
			frame_status_mat[frame_cnt,1] = 5
			
		#separately save only frames with no missing edges	
		#else:
		#	cv.SaveImage(dst_dir + '10degfull/' + 'lefteye/lefteye_frame_' + im.split('_frame_')[1].strip('.npy') + '.png', cv.fromarray(lefteye_image_10deg.copy()))
		#	cv.SaveImage(dst_dir + '10degfull/' + 'righteye/righteye_frame_' + im.split('_frame_')[1].strip('.npy') + '.png', cv.fromarray(righteye_image_10deg.copy()))
		#	cv.SaveImage(dst_dir + '10degfull/' + 'disparityeye_image/disparityeye_image_frame_' + im.split('_frame_')[1].strip('.npy') + '.png', cv.fromarray(disparity_image_10deg.copy()))
		#	numpy.save(dst_dir + '10degfull/' + 'disparityeye/disparityeye_frame_' + im.split('_frame_')[1], disparity_mat_10deg)
		#	cv.SaveImage(dst_dir + '10degfull/' + 'vdisparityeye_image/vdisparityeye_image_frame_' + im.split('_frame_')[1].strip('.npy') + '.png', cv.fromarray(vdisparity_image_10deg.copy()))
		#	numpy.save(dst_dir + '10degfull/' + 'vdisparityeye/vdisparityeye_frame_' + im.split('_frame_')[1], vdisparity_mat_10deg)
				
		
		if targets == 1:
			#store error info
			error_mat[frame_cnt,(7,8)] = lefteye_fixation_pixels #should just be center of the sensor (275,275)
			error_mat[frame_cnt,(9,10)] = lefteye_target_pixels #projection of ground truth target from camera to eye
			error_mat[frame_cnt,(11,12)] = lefteye_camera_target_pixels #projection of camera target to eye
			
			error_mat[frame_cnt,(13,14)] = righteye_fixation_pixels
			error_mat[frame_cnt,(15,16)] = righteye_target_pixels
			error_mat[frame_cnt,(17,18)] = righteye_camera_target_pixels
			
			#disparity of target pixel
			if not numpy.isnan(lefteye_camera_target_pixels[0]):
				camera_target_disparity = disparity_mat[-lefteye_camera_target_pixels[1],lefteye_camera_target_pixels[0]]
				error_mat[frame_cnt,(19)] = camera_target_disparity

			#eye flag
			error_mat[frame_cnt,(20)] = frame_status_mat[frame_cnt,1]
			
			numpy.save(dst_dir + dst_dir.split('/')[2] + '_' + fixation_src_mat + '_target_errors.npy', error_mat)
		
		numpy.save(dst_dir + dst_dir.split('/')[2] + '_' + fixation_src_mat + '_frame_status.npy' , frame_status_mat)
		
		#increment frame cnt
		frame_cnt += 1
		
	print dst_dir.split('/')[2] + '-' + fixation_src_mat
	print 'total frames: ' + str(len(frame_status_mat))
	print 'good fixations:' + '%.2f' % (numpy.true_divide(len(frame_status_mat[frame_status_mat == 1]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[frame_status_mat == 1]))) + ')'
	print 'blinks:' + '%.2f' % (numpy.true_divide(len(frame_status_mat[frame_status_mat == 2]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[frame_status_mat == 2]))) + ')'
	print 'saccades:' + '%.2f' % (numpy.true_divide(len(frame_status_mat[frame_status_mat == 3]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[frame_status_mat == 3]))) + ')'
	print 'missing data:' + '%.2f' % (numpy.true_divide(len(frame_status_mat[frame_status_mat == 4]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[frame_status_mat == 4]))) + ')'
	print 'not full frame:' + '%.2f' % (numpy.true_divide(len(frame_status_mat[frame_status_mat == 5]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[frame_status_mat == 5]))) + ')'
	print 'Nans?:' + '%.2f' % (numpy.true_divide(len(frame_status_mat[numpy.isnan(frame_status_mat)]),len(frame_status_mat))*100) + '% (' + str((len(frame_status_mat[numpy.isnan(frame_status_mat)]))) + ')'
