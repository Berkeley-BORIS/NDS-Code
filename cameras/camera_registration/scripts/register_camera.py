import os.path
import os
import sys
import math
import cv
import cv2
import numpy as np
import CornerPicker
import fnmatch

def load_XML(path_to_params):

    dist_coeffs = cv.Load(os.path.join(path_to_params, "Distortion_cam1.xml"))
    cam_mat = cv.Load(os.path.join(path_to_params, "Intrinsics_cam1.xml"))

    return np.array(cam_mat), np.array(dist_coeffs).squeeze()

def get_board_directories(path, trial=1):
    """Returns a list of paths to the directories containing the chessboards in the order 50cm
    100cm, 450cm."""

    distances = ['50', '100', '450']
    board_directories = []
    for distance in distances:
        directory_path = os.path.join(path, "_".join(["circles", distance, str(trial)]))
        board_directories.append(directory_path)

    return board_directories

def find_circle_centers(board_directories):
    """Returns an array of chessboard corners in an image from each distance of 50cm, 100cm
    and 450cm in that order."""

    print "Finding board corners..."
    pattern_size = (4,11)
    image_coords = []
    corner_images = []
    for directory in board_directories:
        img_name = 'cam1_frame_1.bmp'
        print "Finding circles in", os.path.join(directory, img_name)
        board_image = cv2.imread(os.path.join(directory, img_name),1)  # CHANGED: Loading as RGB
        [pattern_was_found,corners] = cv2.findCirclesGridDefault(board_image,pattern_size,flags=cv2.CALIB_CB_ASYMMETRIC_GRID)
        corners = corners.squeeze()  # CHANGED: moved squeeze to before corner checking instead of after
        if not pattern_was_found and corners == None:
            try:
                corners = np.loadtxt(os.path.join(directory, 'cam1_frame_1_corners.txt'),delimiter=',')
                pattern_was_found = True
                corners = corners.astype('float32')
            except IndexError:
                print 'No corners found! Please find them yourself in Photoshop'
                sys.exit(-1)
        if not pattern_was_found and not corners==None:
            print "Not all corners found! Find them yourself!"
            corners = CornerPicker.main(board_image.copy(), corners, pattern_size)
        corner_image = board_image.copy()
        cv2.drawChessboardCorners(corner_image, pattern_size, corners, pattern_was_found)
        cv2.imwrite('./tmp.png',corner_image)
        image_coords.append(corners)
        corner_images.append(corner_image)

    return np.concatenate(image_coords), corner_images

def find_board_corners(board_directories):
    """Returns an array of chessboard corners in an image from each distance of 50cm, 100cm
    and 450cm in that order."""

    print "Finding board corners..."
    flags = cv2.CALIB_CB_ADAPTIVE_THRESH
    pattern_size = (7,5)
    image_coords = []
    corner_images = []
    for directory in board_directories:
        img_name = 'cam1_frame_1.bmp'
        print "Finding corners in", os.path.join(directory, img_name)
        #board_image = cv2.imread(os.path.join(directory, img_name), cv2.CV_LOAD_IMAGE_COLOR)
        board_image = cv2.imread(os.path.join(directory, img_name),0)
        (pattern_was_found, corners) = cv2.findChessboardCorners(board_image, pattern_size, flags=flags)
        if pattern_was_found:
            cv2.cornerSubPix(board_image, corners, (4, 4),  (-1, -1), (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 100, 1e-7))
        if not pattern_was_found and corners == None:
            try:
                corners = np.loadtxt(os.path.join(directory, 'cam1_frame_1_corners.txt'),delimiter=',')
                pattern_was_found = True
                corners = corners.astype('float32')
            except IndexError:
                print 'No corners found! Please find them yourself in Photoshop'
                sys.exit(-1)
        if not pattern_was_found and not corners==None:
            print "Not all corners found! Find them yourself!"
            corners = CornerPicker.main(board_image.copy(), corners)
        corner_image = board_image.copy()
        corners = corners.squeeze()
        cv2.drawChessboardCorners(corner_image, pattern_size, corners, pattern_was_found)
        image_coords.append(corners)
        corner_images.append(corner_image)

    return np.concatenate(image_coords), corner_images

def get_circle_world_coords_EMILY(dims,num_images,square_size):
    '''determine 3d object points for each image
    '''

    dims = (4,11)
    width = dims[0]
    height = dims[1]
    num_pts = width*height
    inter_circle_distance = 52.25*2.0  # Pixel distance between circles in the same row
    x_offset = 75.25  # start 75.25 pixels in from the left side of the board png
    y_offset = 69.50  # start 383.50 pixels down from the top side of the board png

    #circle pixel coords with 0,0 at top left circle
    opts_pix = np.zeros((num_pts,2))
    for i in range(height):
        for j in range(width):
                if i%2==0:
                    opts_pix[i*width+j,0] = (i*(inter_circle_distance/2.00))
                    opts_pix[i*width+j,1] = j*inter_circle_distance
                    opts_pix[i*width+j,2] = 0
                else:
                    opts_pix[i*width+j,0] = (i*(inter_circle_distance/2.00))
                    opts_pix[i*width+j,1] = (j*inter_circle_distance) + inter_circle_distance/2.00
                    opts_pix[i*width+j,2] = 0

    opts_pix[0] = opts_pix[0] + x_offset
    opts_pix[1] = opts_pix[1] + y_offset
    #opts = np.array(opts, dtype = np.float32)

    return opts
        
def get_circle_world_coords():
    """
    This returns an array of the world coordinates of the centers of the CIRCLES
    in the CIRCLEGRID. Pixel offsets were measured in photoshop and are therefor hard coded.
    This function is specific to the acircles_pattern_grey_DIST.png files used to display the
    circle grid in the experiment. If those images were not used on the screen, DO NOT USE THIS
    FUNCTION.
    """

    x_offset = 597.5  # start 597.5 pixels in from the left side of the board png
    y_offset = 69.50  # start 69.50 pixels down from the top side of the board png
    inter_circle_distance = 52.25*2.0  # Pixel distance between circles in the same row
    center_points_pix_png = []  # List to fill with png pixel positions
    shift_direction = 1  # When we move to the right to calculate the next set of centers, we need to shift the y_offset in the negative direction

    # ATTENTION! ORDERING COULD BE MESSED UP!! Looking at the calibration frames, it looks like a
    # (4,11) sized board starts with the LOWER LEFT circle and works its way up so that is what I am assuming!!
    size = (11,4)  # THIS SIZE IS NOT REFLECTIVE OF BOARD ORIENTATION! YES, I MEAN TO DO 11 X-POSITIONS AND 4 Y-POSITIONS TO CALCULATE THE WORLD COORDINATES!!
    for x_pos in xrange(size[0]):
        for y_pos in xrange(size[1]):
            next_point = (x_offset - x_pos*inter_circle_distance/2.0,  # Move over to the left half the inter-circle distance each time we want to
                          y_offset + y_pos*inter_circle_distance)  # Move down by the intercircle distance
            center_points_pix_png.append(next_point)

        y_offset += shift_direction * inter_circle_distance/2.0  # We need to change the y offset when we move to the next x position
        shift_direction *= -1  # Next time we'll want to shift the y_offset in the opposite direction

    center_points_pix_png = np.array(center_points_pix_png)  # Turn this into an array


    center_points_cm_world = np.zeros((size[0]*size[1]*3, 3), dtype='float')
    start_idx = 0
    end_idx = size[0] * size[1]
    for dist in [50, 100, 450]:
        if dist == 50:
            scale = .5
            top_left = ((1920-337)/2.0, (1024-250)/2.0)  # top left corner of the 50 png in pixels with (0,0) at center, minus y is UP
        elif dist == 100:
            scale = 1.0
            top_left = ((1920-674)/2.0, (1024-500)/2.0)
        elif dist == 450:
            scale = 2.0
            top_left = ((1920-1349)/2.0, (1024-1000)/2.0)

        center_points_pix_screen = top_left + (center_points_pix_png * scale)  # shift the grid to be in the center of the screen and scale them
        #center_points_pix_screen = (center_points_pix_png * scale) - top_left
    
        center_points_pix_screen -= (1920/2.0, 1024/2.0)
        print "circle centers in screen pixels:\n", center_points_pix_screen

        cm_per_pix_x = 121.0/1920  # LG 3D TV width in cm over pixels wide
        cm_per_pix_y = 68.0/1024
        center_points_cm_world[start_idx:end_idx, 0] = cm_per_pix_x * center_points_pix_screen[:,0]  # Convert x position to cm
        center_points_cm_world[start_idx:end_idx, 1] = cm_per_pix_y * -center_points_pix_screen[:,1]  # Convert y position to cm, remember POSITIVE Y is now UP!
        center_points_cm_world[start_idx:end_idx, 2] = float(dist)

        start_idx += size[0]*size[1]
        end_idx += size[0]*size[1]

    print "circle centers in cm:\n", center_points_cm_world

    return center_points_cm_world

def get_world_coords():
    """Returns an array of the world coordinates of chessboard intersections from each distance
    in the order 50cm, 100cm, 450cm"""

    print "Generating world coordinates..."
    world_coords = []
    for dist in [50, 100, 450]:
        (model_x, model_y) = np.meshgrid(np.arange(-3,4),
                                         np.arange(2,-3,-1))
        model_coords = np.column_stack((model_x.flatten(),
                                        model_y.flatten(),
                                        np.zeros((35,1))))
        if dist == 50:
            square_size = 5.24/2
        else:
            square_size = 5.24
        model_coords = model_coords * square_size
        model_coords[:,2] = dist
        world_coords.append(model_coords)


    return np.concatenate(world_coords)

def project_to_camera_pixels(intrinsic,threedcoords):
    '''takes in camera intrinsics and 3d points in world and projects points to pixels
    '''

    #convert intrinsic focal length from pixels to cm
    focallength_x_cm = (intrinsic[0,0]*(6.0/10000.0))
    focallength_y_cm = (intrinsic[1,1]*(6.0/10000.0))

    #perspective projection
    threedcoords_x_cm = (focallength_x_cm*(threedcoords[0,:]))/threedcoords[2,:]
    threedcoords_y_cm = (focallength_y_cm*threedcoords[1,:])/threedcoords[2,:]

    #convert to pixels
    threedcoords_x_pix = threedcoords_x_cm*(10000.0/6.0) + intrinsic[0,2]
    threedcoords_y_pix = threedcoords_y_cm*(10000.0/6.0) + intrinsic[1,2]

    pixels = np.append(threedcoords_x_pix,threedcoords_y_pix,axis=1)

    return pixels

def plot_circle_centers(board_directories, save_path, world_coords, image_coords, rvec, tvec, cam_mat, dist_coeffs):
    
    start = 0
    end = 44
    for directory in board_directories:
        img = cv2.imread(os.path.join(directory, 'cam1_frame_1.bmp'), 1)
        dist_world_coords = world_coords[start:end, :]
        
        if start == 0:
            dist = 50.0
        elif start == 44:
            dist = 100.0
        elif start == 88:
            dist = 450.0

        screen_edges = np.array([[-121/2.0, 68/2.0, dist],
                                 [121/2.0, 68/2.0, dist],
                                 [-121/2.0, -68/2.0, dist],
                                 [121/2.0, -68/2.0, dist]])
        
        dist_world_coords = np.concatenate((dist_world_coords, screen_edges))
        camera_coords, jac = cv2.projectPoints(objectPoints=dist_world_coords,rvec=rvec,tvec=tvec,cameraMatrix=cam_mat,distCoeffs=dist_coeffs)
        camera_coords = camera_coords.squeeze()
        for coords in camera_coords:
            cv2.circle(img, tuple(coords), 3, (0,0,255), -1)
        for coords in image_coords[start:end]:
            cv2.circle(img, tuple(coords), 5, (180,255,0), 2)
        cv2.imwrite(os.path.join(save_path, 'reprojection_'+str(int(dist))+'.png'), img)
        
        start += 44
        end += 44
        
        

if __name__=="__main__":
    try:
        path_to_exp = os.path.abspath(sys.argv[1])
    except IndexError:
        print "Input the paths to the experiment_images and calibration_parameters directories of interest and, optionally, the trial number."
        sys.exit(-1)
    try:
        path_to_params = os.path.abspath(sys.argv[2])
        for dir in os.listdir(path_to_params):
            if fnmatch.fnmatch(dir,'calibration_frames*'):
                path_to_params = path_to_params + '/' + dir + '/'
                break
    except IndexError:
        print "Input the paths to the experiment_images and calibration_parameters directories of interest and, optionally, the trial number."
        sys.exit(-1)
    try:
        trial = sys.argv[3]
    except IndexError:
        trial = 1
    board_directories = get_board_directories(path_to_exp, trial)
    image_coords, corner_images = find_circle_centers(board_directories)
    world_coords = get_circle_world_coords()
    cam_mat, dist_coeffs = load_XML(path_to_params=path_to_params)

    retval, rvec, tvec = cv2.solvePnP(objectPoints=world_coords,
                                      imagePoints=image_coords,
                                      cameraMatrix=cam_mat,
                                      distCoeffs=dist_coeffs,
                                      useExtrinsicGuess=0)
    
	#use guess of rvec and tvec
    #rvec = np.asarray([ 0*(math.pi/180), -180*(math.pi/180) , 0*(math.pi/180)])
    #tvec = np.asarray([-3.25 ,0,  0])

    #q = np.array(world_coords, dtype=np.float32, copy=True)
    #rvec, tvec, inliers = cv2.solvePnPRansac(objectPoints=q, imagePoints=image_coords, cameraMatrix=cam_mat, distCoeffs=dist_coeffs, useExtrinsicGuess=0)
    #print "inliers:", inliers

    #rvec = np.array([5.71195205*(math.pi/180),-176.26498361*(math.pi/180),-12.66819407*(math.pi/180)])
    #tvec = np.array([-4.34928759,-7.41820741, 6.38744438])
    #rvec, tvec = cv2.solvePnP(objectPoints=world_coords, imagePoints=image_coords, cameraMatrix=cam_mat, distCoeffs=dist_coeffs,rvec = rvec, tvec=tvec, useExtrinsicGuess=1)

    camera_coords, jac = cv2.projectPoints(objectPoints=world_coords,rvec=rvec,tvec=tvec,cameraMatrix=cam_mat,distCoeffs=dist_coeffs)
    camera_coords = np.reshape(camera_coords,(132,2))
    error = np.sqrt(np.power(camera_coords[:,0]-image_coords[:,0],2) + np.power(camera_coords[:,1]-image_coords[:,1],2))
    rms = np.sqrt(np.mean(np.power(error,2)))

    print "solvepnp rms error: " + str(rms)
    print "tvec: " + str(tvec)
    print "rvec: " + str(rvec*(180/math.pi))

    print "Transformations found! Saving file..."

    # Save data
    (root, tail) = os.path.split(path_to_exp)
    experiment_name = tail.split('_')[0]
    path_to_save_data = "/Users/natdispstats/Documents/cameras/camera_registration/data/" + experiment_name + "/"
    if not os.path.exists(path_to_save_data):
        os.mkdir(path_to_save_data)
    save_file = os.path.join(path_to_save_data, experiment_name + "_transform_" + str(trial) + ".npz")
    np.savez(save_file, rvec=rvec, tvec=tvec, cam_mat=cam_mat, dist_coeffs=dist_coeffs)

    # Save chessboard images with corners drawn
    dist_names = ['50', '100', '450']
    for img, i in zip(corner_images, xrange(len(corner_images))):
        fname = "checkerboard_" + dist_names[i] + "_" + str(trial) + ".png"
        cv2.imwrite(os.path.join(path_to_save_data, fname), img)
    print "File saved to", save_file

    # save reprojection images
    plot_circle_centers(board_directories, path_to_save_data, world_coords, image_coords, rvec, tvec, cam_mat, dist_coeffs)

    # TODO Save the CLI arguments so we have a record of what
    # produced this data

    print "\nSaving rms reprojection error..."
    rms_file = open(path_to_save_data+ 'registration_rms_error.txt', 'w')
    rms_file.write("\n rms pixel error: " + str(rms))
    rms_file.close()
