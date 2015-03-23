import cv
import os
import string
import pp
import fnmatch

def rectify_frames(subj,parallel_input):
	''' load stereo calibration and rectify images in source dir
	http://code.google.com/p/tjpstereovision/source/browse/code/Python/stereorectify.py?spec=svnba9b04d5082c96c46860edca01e9133c93c53ead&r=ba9b04d5082c96c46860edca01e9133c93c53ead
	example call:  disparity_map_utils.rectify_frames("../data/frames/","../data/frames_rect/")
	'''
	
	#source directory for unrectified camera frames
	src_dir = "../../raw_data/" + subj + "_images/" + parallel_input + "/"
	dst_dir = "../../image_rectification/data/" + subj + "/" + subj + "_" + parallel_input + "_frames_rect/"
	
	#grab calbration frames directory
	cam_dir = "../../stereo_calibration/data/"  + subj + "/"
	for dir in os.listdir(cam_dir):
		if fnmatch.fnmatch(dir,'calibration_frames*'):
			cam_dir = cam_dir + dir + '/'
			break
			
	#load stereo camera calibration
	map1x = cv.Load(cam_dir + "Rectification_map_cam1x.xml")
	map1y = cv.Load(cam_dir + "Rectification_map_cam1y.xml")
	map2x = cv.Load(cam_dir + "Rectification_map_cam2x.xml")
	map2y = cv.Load(cam_dir + "Rectification_map_cam2y.xml")
		
	#create destination directory if it doesn't exist
	if not os.path.exists(string.join(string.split(dst_dir,'/')[0:-2],'/')):
		os.mkdir(string.join(string.split(dst_dir,'/')[0:-2],'/'))
	if not os.path.exists(dst_dir):
		os.mkdir(dst_dir)
		
	#iterate through images in source directory
	imlist = os.listdir(src_dir)
	
	for im in imlist:
		
		if im.find(".bmp") > 0:
			
			im_tmp = cv.LoadImage(src_dir + im, cv.CV_LOAD_IMAGE_GRAYSCALE)
			r = cv.CloneImage(im_tmp)
		
			#apply left or right camera rectification
			if im.find('cam1') == 0:
				cv.Remap(im_tmp, r, map1x, map1y)
			elif im.find('cam2') == 0:
				cv.Remap(im_tmp, r, map2x, map2y)
		
			#fix numbering of frames in source directory to add leading zeros. this will be helpful down the line
			im_new = "rect_" + im.split((im.split("frame")[1].strip(".bmp")))[0] + "_" + "0"*(6-len(im.split("_frame_")[1].strip(".bmp"))) + im.split("_frame_")[1]
		
			cv.SaveImage(dst_dir + im_new,r)
			
	#save parameters into text file in destination directory
	params = open(dst_dir+"rectification_params.txt", "w")
	params.write("\n calibration directory: " + cam_dir)
	params.write("\n frames source directory: " + src_dir)
	params.close()