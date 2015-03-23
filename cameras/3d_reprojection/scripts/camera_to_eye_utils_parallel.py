import cv
import os.path
import cv2
import numpy
import scipy.io
import math
import os

def find_pts_near_fovea(coords):

    coords_x = coords[0,:].reshape((480,640))
    coords_y = coords[1,:].reshape((480,640))
    coords_z = coords[2,:].reshape((480,640))

    xy_mag = numpy.sqrt(coords_x**2 + coords_y**2)
    foveal_ind = numpy.where(xy_mag == numpy.nanmin(xy_mag))

    pix_buffer = 1
    fovea_pts_x = coords_x[foveal_ind[0][0]-pix_buffer:foveal_ind[0][0]+pix_buffer+1, foveal_ind[1][0]-pix_buffer:foveal_ind[1][0]+pix_buffer+1]
    fovea_pts_y = coords_y[foveal_ind[0][0]-pix_buffer:foveal_ind[0][0]+pix_buffer+1, foveal_ind[1][0]-pix_buffer:foveal_ind[1][0]+pix_buffer+1]
    fovea_pts_z = coords_z[foveal_ind[0][0]-pix_buffer:foveal_ind[0][0]+pix_buffer+1, foveal_ind[1][0]-pix_buffer:foveal_ind[1][0]+pix_buffer+1]

    # coords = coords[:,~numpy.isnan(coords).all(axis=0)]

    # foveal_inds = numpy.where((numpy.abs(coords[0,:]) < 5) & (numpy.abs(coords[1,:]) < 5))

    return numpy.vstack((fovea_pts_x.ravel(), fovea_pts_y.ravel(), fovea_pts_z.ravel()))

def make_destination_directory_tree(dst_dir):
    '''makes a directory tree for results of reprojection
    '''
    os.mkdir(dst_dir)
    os.mkdir(dst_dir + 'lefteye/')
    os.mkdir(dst_dir + 'righteye/')
    os.mkdir(dst_dir + 'cycloeye/')
    # os.mkdir(dst_dir + 'disparityeye/')
    # os.mkdir(dst_dir + 'disparityeye_image/')
    # os.mkdir(dst_dir + 'vdisparityeye/')
    # os.mkdir(dst_dir + 'vdisparityeye_image/')
    os.mkdir(dst_dir + 'disparity_cycloeye/')
    os.mkdir(dst_dir + 'disparity_cycloeye_image/')
    os.mkdir(dst_dir + 'vdisparity_cycloeye/')
    os.mkdir(dst_dir + 'vdisparity_cycloeye_image/')
    os.mkdir(dst_dir + 'deptheye/')
    os.mkdir(dst_dir + 'cyclopts_info/')

def save_reprojection_parameters(dst_dir,transform_src_file,depth_src_dir,image_src_dir,fixation_src_file):
    '''stores source files used for reprojection
    '''

    params = open(dst_dir+'reprojection_params.txt', 'w')
    params.write("\n transform source file: " + transform_src_file)
    params.write("\n 3d camera points source directory: " + depth_src_dir)
    params.write("\n camera images source directory: " + image_src_dir)
    params.write("\n fixation points source file: " + fixation_src_file)
    params.close()

def load_and_invert_transforms(transform_src_file):
    '''loads in eye_to_cam transforms from registration and inverts to get cam_to_eye
    '''

    transforms = numpy.load(transform_src_file)
    rvec_eye_to_cam = transforms['rvec']
    R_eye_to_cam, jac = cv2.Rodrigues(rvec_eye_to_cam)
    R_cam_to_eye = numpy.asarray(numpy.asmatrix(R_eye_to_cam).I)
    rvec_cam_to_eye, jac  = cv2.Rodrigues(R_cam_to_eye)
    tvec_eye_to_cam = transforms['tvec']
    tvec_cam_to_eye = tvec_eye_to_cam*-1
    intrinsic_mat = transforms['cam_mat']

    return R_eye_to_cam,R_cam_to_eye,tvec_eye_to_cam,tvec_cam_to_eye,intrinsic_mat

def remove_image_edges(image,edge):
    '''replaces values at edge of image (noisy ones) with nans
    '''

    image[0:edge-1,:,:] = numpy.nan
    image[480-edge:480,:,:] = numpy.nan
    image[:,0:edge-1,:] = numpy.nan
    image[:,640-edge:640,:] = numpy.nan

    return image

def ground_truth_target_cyclopean_coords(target_coords,frame):
    '''grabs the ground truth cyclopean coordinates of a target
    '''

    target_point = target_coords[frame,:]
    #set (0,0,0) targets to nans
    if target_point[2] == 0:
        target_point = [numpy.nan , numpy.nan , numpy.nan]

    return target_point

def camera_target_camera_coords(target_point,target_camera_pixels,camera_coords,center_x,center_y,focal_length_x,focal_length_y,frame):
    '''if there is a ground truth target and it was seen by the camera, also grab camera pixel location of target
    from EPicker
    '''

    if numpy.isnan(target_point[2]) == 0 and bool(target_camera_pixels.get(frame)):
        #note that pixels are given in x,y here, so to index 3d coordinates we take y then the x (row then column)
        target_camera_pixel = target_camera_pixels.get(frame)
        target_camera_coords = camera_coords[target_camera_pixel[1],target_camera_pixel[0],:]
        #for points with inf distance, convert to just very very large distance
        #for points at 0,0,0 placeholder, replace with nans
        if target_camera_coords[2] == numpy.inf or target_camera_coords[2] == -numpy.inf:
            target_camera_coords = [ -1e100 * (-(center_x-target_camera_pixel[0])/focal_length_x) , -1e100 * (-(center_y-target_camera_pixel[1])/focal_length_y) , -1e100 ]
        elif target_camera_coords[2] == 0:
            target_camera_coords = [ numpy.nan , numpy.nan , numpy.nan ]

    else:
        target_camera_coords = [ numpy.nan , numpy.nan , numpy.nan ]

    return target_camera_coords

def initialize_images_and_mats():
    '''creates empty images and data matrices to fill up with data from each frame
    '''

    lefteye_image = numpy.zeros((550,550))
    righteye_image = numpy.zeros((550,550))
    cycloeye_image = numpy.zeros((550,550))
    disparity_image = numpy.zeros((550,550,3))
    vdisparity_image = numpy.zeros((550,550,3))
    disparity_mat = numpy.ones((550,550))*numpy.nan
    vdisparity_mat = numpy.ones((550,550))*numpy.nan
    depth_mat = numpy.ones((550,550,4))*numpy.nan
    disparity_cyclo_image = numpy.zeros((550,550,3))
    vdisparity_cyclo_image = numpy.zeros((550,550,3))
    disparity_cyclo_mat = numpy.ones((550,550))*numpy.nan
    vdisparity_cyclo_mat = numpy.ones((550,550))*numpy.nan

    return lefteye_image,righteye_image,cycloeye_image,disparity_image,vdisparity_image,disparity_mat,vdisparity_mat,depth_mat,disparity_cyclo_image,vdisparity_cyclo_image, disparity_cyclo_mat, vdisparity_cyclo_mat

def project_to_pixels(intrinsic,threedcoords,calc_radial_distance=0):
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

    #calculate radial distance to each point
    if calc_radial_distance == 1:
        threedcoords_X = threedcoords[0,:]
        threedcoords_Y = threedcoords[1,:]
        threedcoords_Z = threedcoords[2,:]
        radial_distance = numpy.sqrt(numpy.square(threedcoords_X) + numpy.square(threedcoords_Y) + numpy.square(threedcoords_Z))



    #reshape
    #if projecting a single point:
    if threedcoords.size == 3:
        pixels = numpy.append(threedcoords_x_pix,threedcoords_y_pix)
    #if projecting a camera image
    else:
        #deal with occlusions
        #create and sort matrix to find occlusions
        occlu_mat = numpy.array((numpy.round(threedcoords_x_pix),numpy.round(threedcoords_y_pix),threedcoords[2,:]))
        #sort matrix in ascending x, then ascending y, then ascending z
        occlu_mat_sort_ind = numpy.lexsort((threedcoords[2,:],numpy.round(threedcoords_y_pix),numpy.round(threedcoords_x_pix)))
        occlu_mat_sorted = occlu_mat[:,occlu_mat_sort_ind]
        #iterate through matrix,taking only first columns with any given pixel x,y (i.e., closest z)
        x0 = occlu_mat_sorted[0,0]
        y0 = occlu_mat_sorted[1,0]
        #if rounded pixel values are repeated, set pixels to nans
        for k in range(1,len(occlu_mat_sorted[0,:])):
            if occlu_mat_sorted[0,k] == x0 and occlu_mat_sorted[1,k] == y0:
                threedcoords_x_pix[occlu_mat_sort_ind[k]] = numpy.nan
                threedcoords_y_pix[occlu_mat_sort_ind[k]] = numpy.nan
                if calc_radial_distance == 1:
                    threedcoords_X[occlu_mat_sort_ind[k]] = numpy.nan
                    threedcoords_Y[occlu_mat_sort_ind[k]] = numpy.nan
                    threedcoords_Z[occlu_mat_sort_ind[k]] = numpy.nan
                    radial_distance[occlu_mat_sort_ind[k]] = numpy.nan
            #otherwise, move on the the next pixel values
            else:
                x0 = occlu_mat_sorted[0,k]
                y0 = occlu_mat_sorted[1,k]


        threedcoords_x_pix = numpy.reshape(threedcoords_x_pix,(480,640))
        threedcoords_y_pix = numpy.reshape(threedcoords_y_pix,(480,640))
        pixels = numpy.append(numpy.reshape(threedcoords_x_pix,(480,640,1)),numpy.reshape(threedcoords_y_pix,(480,640,1)),axis=2)
        if calc_radial_distance == 1:
            threedcoords_X = numpy.reshape(threedcoords_X,(480,640))
            threedcoords_Y = numpy.reshape(threedcoords_Y,(480,640))
            threedcoords_Z = numpy.reshape(threedcoords_Z,(480,640))
            radial_distance = numpy.reshape(radial_distance,(480,640))
            distance = numpy.append(numpy.reshape(threedcoords_X,(480,640,1)),numpy.reshape(threedcoords_Y,(480,640,1)),axis=2)
            distance = numpy.append(distance,numpy.reshape(threedcoords_Z,(480,640,1)),axis=2)
            distance = numpy.append(distance,numpy.reshape(radial_distance,(480,640,1)),axis=2)


    if calc_radial_distance == 1:
        return distance
    else:
        return pixels

def create_eye_image(eye_pixels,cam_image,eye_image):
    '''takes in look up table of camera to eye projections and camera image and creates eye images
    '''

    #vectorize x coords, y coords, and greyscale values
    eye_xs_all = numpy.reshape(eye_pixels[:,:,0],eye_pixels[:,:,0].size)
    eye_ys_all = numpy.reshape(eye_pixels[:,:,1],eye_pixels[:,:,1].size)
    eye_cam_all = numpy.reshape(cam_image,cam_image.size)

    #set selection criterian, ie, eye sensor size
    logic1 = numpy.logical_and(eye_xs_all > 0, eye_xs_all < 549)
    logic2 = numpy.logical_and(eye_ys_all > 0, eye_ys_all < 549)

    #apply criteria and convert to ints
    eye_xs = numpy.round(eye_xs_all[numpy.where(numpy.logical_and(logic1, logic2))])
    eye_ys = numpy.round(eye_ys_all[numpy.where(numpy.logical_and(logic1, logic2))])
    eye_cam = numpy.round(eye_cam_all[numpy.where(numpy.logical_and(logic1, logic2))])


    eye_image[-eye_ys.astype('int'),eye_xs.astype('int')] = eye_cam.astype('int')

    return eye_image

def create_depth_image(eye_pixels,cyclo_depth,eye_depth):
    '''takes in look up table of camera to eye projections, creates matrix with XYZ and radial distance of points from
    cyclopean eye to the left eye's view
    '''

    #vectorize x coords, y coords, and greyscale values
    eye_xs_all = numpy.reshape(eye_pixels[:,:,0],eye_pixels[:,:,0].size)
    eye_ys_all = numpy.reshape(eye_pixels[:,:,1],eye_pixels[:,:,1].size)
    eye_X = numpy.reshape(cyclo_depth[:,:,0],cyclo_depth[:,:,0].size)
    eye_Y = numpy.reshape(cyclo_depth[:,:,1],cyclo_depth[:,:,1].size)
    eye_Z = numpy.reshape(cyclo_depth[:,:,2],cyclo_depth[:,:,2].size)
    eye_D = numpy.reshape(cyclo_depth[:,:,3],cyclo_depth[:,:,3].size)

    #set selection criterian, ie, eye sensor size
    logic1 = numpy.logical_and(eye_xs_all > 0, eye_xs_all < 549)
    logic2 = numpy.logical_and(eye_ys_all > 0, eye_ys_all < 549)

    #apply criteria and convert to ints
    eye_xs = numpy.round(eye_xs_all[numpy.where(numpy.logical_and(logic1, logic2))])
    eye_ys = numpy.round(eye_ys_all[numpy.where(numpy.logical_and(logic1, logic2))])
    eye_Xs = eye_X[numpy.where(numpy.logical_and(logic1, logic2))]
    eye_Ys = eye_Y[numpy.where(numpy.logical_and(logic1, logic2))]
    eye_Zs = eye_Z[numpy.where(numpy.logical_and(logic1, logic2))]
    eye_Ds = eye_D[numpy.where(numpy.logical_and(logic1, logic2))]

    eye_depth[-eye_ys.astype('int'),eye_xs.astype('int'),0] = eye_Xs
    eye_depth[-eye_ys.astype('int'),eye_xs.astype('int'),1] = eye_Ys
    eye_depth[-eye_ys.astype('int'),eye_xs.astype('int'),2] = eye_Zs
    eye_depth[-eye_ys.astype('int'),eye_xs.astype('int'),3] = eye_Ds

    return eye_depth

def create_HH_table(intrinsic_eye):
    """
    Takes in eye's intrinsic matrix and creates a lookup table
    for converting each pixel location to Helmholtz elevation and azimuth
    in degrees.
    """
    """Takes in a 3x2 array the specifies the 3D coordinates of the left and right eye's
    fixation points respectively in cyclopean coordinates and outputs each eye's
    Helmholtz coordinates for that point.
    """

    if os.path.exists('HHcoords_of_eye_sensor_pixels.npy'):
        HH_table = numpy.load('HHcoords_of_eye_sensor_pixels.npy')
    else:

        x,y = numpy.meshgrid(range(0,550),range(0,550))
        eye_sensor_pixels = numpy.concatenate((numpy.reshape(x,(550,550,1)),numpy.reshape(y,(550,550,1)),intrinsic_eye[0,0]*numpy.ones((550,550,1))),axis=2)
        eye_sensor_cms    = eye_sensor_pixels*(6.0/10000.0)
        eye_sensor_cms = eye_sensor_cms - numpy.append([eye_sensor_cms[intrinsic_eye[0,2],intrinsic_eye[1,2],0:2]],0)
        #eye_sensor_cms_cols = numpy.reshape(eye_sensor_cms,(550*550,3),order='C').T

        HH_table = numpy.nan*numpy.ones((550,550,2))
        #phi = numpy.nan*numpy.ones((550,550,1))
        #theta = numpy.nan*numpy.ones((550,550,1))

        z_dir = numpy.array([0, 0, 1], dtype='double')

        #for point in eye_sensor_cms_cols.T:
        for xpoint in xrange(0,550):
            for ypoint in xrange(0,550):

                point = eye_sensor_cms[xpoint,ypoint,:]

                # Get the component of the  point that lies in the zy-plane
                zy_component = numpy.array(point)
                #zy_component = zy_component.squeeze()
                zy_component[0] = 0
                # Calculate the angle between the zy-component and the z-direction
                phi_cos_ratio = numpy.dot(z_dir, zy_component) / (numpy.linalg.norm(z_dir) * numpy.linalg.norm(zy_component))
                #if phi_cos_ratio > 1:
                #    print "Phi cos ratio is", phi_cos_ratio
                abs_phi = numpy.arccos(numpy.min([1, phi_cos_ratio]))
                # Correct the angle to be positive if looking down and negative if looking up
                phi_norm = numpy.cross(z_dir, zy_component)
                if phi_norm[0]:
                    phi_dir = phi_norm[0] / numpy.abs(phi_norm[0])
                else:
                    phi_dir = 1
                HH_table[xpoint,ypoint,0] = phi_dir * abs_phi
                #phi.append(phi_dir * abs_phi)

                # Calculate the angle between the zy-component and the fixation vector
                theta_cos_ratio = numpy.dot(zy_component.T, point) / (numpy.linalg.norm(zy_component) * numpy.linalg.norm(point))
                #if theta_cos_ratio > 1:
                #    print "Theta cos ratio is", theta_cos_ratio
                abs_theta = numpy.arccos(numpy.min([1, theta_cos_ratio]))
                # Correct the angle to positive if looking left and negative if looking right
                theta_norm = numpy.cross(point, zy_component)
                if theta_norm[1]:
                    theta_dir = theta_norm[1] / numpy.abs(theta_norm[1])
                else:
                    theta_dir = 1
                HH_table[xpoint,ypoint,1] = theta_dir * abs_theta

        HH_table = HH_table*(180/math.pi)
        numpy.save('HHcoords_of_eye_sensor_pixels.npy', HH_table)
                #theta.append(theta_dir * abs_theta)

    return HH_table


def convert_pixel_loc_to_HH(eye_pixels, intrinsic_eye):
    """
    Takes in the Helmholtz lookup table and the lookup table of
    camera to eye projections and creates a camera to eye helmholtz
    matrix.
    """

    #eye_pixels = numpy.concatenate((eye_pixels[:,:,0],eye_pixels[:,:,1],intrinsic_eye[0,0]*numpy.ones((480,640,1))),axis=2)
    eye_pixels = numpy.concatenate((eye_pixels,intrinsic_eye[0,0]*numpy.ones((480,640,1))),axis=2)
    eye_cms   = eye_pixels*(6.0/10000.0)

    eye_center_pixels = numpy.array([intrinsic_eye[0,2], intrinsic_eye[1,2], 0])
    eye_center_cms = eye_center_pixels*(6.0/10000.0)

    eye_cms = eye_cms - eye_center_cms
    #eye_cms = eye_cms - numpy.append([eye_cms[intrinsic_eye[0,2],intrinsic_eye[1,2],0:2]],0)
    #eye_sensor_cms_cols = numpy.reshape(eye_sensor_cms,(550*550,3),order='C').T

    eye_HH = numpy.nan*numpy.ones((480,640,2))
    #eye_HH = numpy.ones((480,640,2))
    #phi = numpy.nan*numpy.ones((550,550,1))
    #theta = numpy.nan*numpy.ones((550,550,1))

    z_dir = numpy.array([0, 0, 1], dtype='double')

    #for point in eye_sensor_cms_cols.T:
    for xpoint in xrange(0,640):
        for ypoint in xrange(0,480):

            point = eye_cms[ypoint,xpoint,:]

            # Get the component of the  point that lies in the zy-plane
            zy_component = numpy.array(point)
            #zy_component = zy_component.squeeze()
            zy_component[0] = 0
            # Calculate the angle between the zy-component and the z-direction
            phi_cos_ratio = numpy.dot(z_dir, zy_component) / (numpy.linalg.norm(z_dir) * numpy.linalg.norm(zy_component))
            #if phi_cos_ratio > 1:
            #    print "Phi cos ratio is", phi_cos_ratio
            abs_phi = numpy.arccos(numpy.min([1, phi_cos_ratio]))
            # Correct the angle to be positive if looking down and negative if looking up
            phi_norm = numpy.cross(z_dir, zy_component)
            if phi_norm[0]:
                phi_dir = phi_norm[0] / numpy.abs(phi_norm[0])
            else:
                phi_dir = 1
            eye_HH[ypoint,xpoint,1] = phi_dir * abs_phi
            #elevation

            # Calculate the angle between the zy-component and the fixation vector
            theta_cos_ratio = numpy.dot(zy_component.T, point) / (numpy.linalg.norm(zy_component) * numpy.linalg.norm(point))
            #if theta_cos_ratio > 1:
            #    print "Theta cos ratio is", theta_cos_ratio
            abs_theta = numpy.arccos(numpy.min([1, theta_cos_ratio]))
            # Correct the angle to positive if looking left and negative if looking right
            theta_norm = numpy.cross(point, zy_component)
            if theta_norm[1]:
                theta_dir = theta_norm[1] / numpy.abs(theta_norm[1])
            else:
                theta_dir = 1
            eye_HH[ypoint,xpoint,0] = theta_dir * abs_theta
            #azimuth

    eye_HH = eye_HH*(180/math.pi)

    return eye_HH

def convert_pixel_loc_to_HH_opt(eye_pixels, intrinsic_eye):
    """
    Takes in the Helmholtz lookup table and the lookup table of
    camera to eye projections and creates a camera to eye helmholtz
    matrix. FASTER VERSION
    """

    #convert from pixels to cycloeans cms
    eye_pixels = numpy.concatenate((eye_pixels,intrinsic_eye[0,0]*numpy.ones((480,640,1))),axis=2)
    eye_cms   = eye_pixels*(6.0/10000.0)
    eye_center_pixels = numpy.array([intrinsic_eye[0,2], intrinsic_eye[1,2], 0])
    eye_center_cms = eye_center_pixels*(6.0/10000.0)
    eye_cms = eye_cms - eye_center_cms
    #flip y axis sign to agree with right,up,forward positive cyclopean coords
    eye_cms[:,:,1] = eye_cms[:,:,1]*-1

    #allocate space for HH angles matrix
    eye_HH = numpy.nan*numpy.ones((480,640,2))

    #define z axis
    z_dir = numpy.array([0, 0, 1], dtype='double')
    # Get the component of the points that lie in the zy-plane
    zy_component = eye_cms.copy()
    zy_component[:,:,0] = 0

    #Elevation angle
    # Calculate the angle between the zy-component and the z-direction
    phi_cos_ratio = zy_component[:,:,2] / numpy.sqrt(numpy.sum(numpy.square(zy_component),axis=2))
    #make sure phi < 1
    phi_cos_ratio[phi_cos_ratio > 1] = 1
    #convert phi to angle
    abs_phi = numpy.arccos(phi_cos_ratio)
    #correct the angle to be positive if looking down and negative if looking up
    eye_HH[:,:,1] = -numpy.sign(zy_component[:,:,1]) * abs_phi

    #Azimuth angle
    # Calculate the angle between the zy-component and the fixation vector
    theta_cos_ratio = ((zy_component[:,:,1]*eye_cms[:,:,1]) + (zy_component[:,:,2]*eye_cms[:,:,2])) / (numpy.sqrt(numpy.sum(numpy.square(zy_component),axis=2)) * numpy.sqrt(numpy.sum(numpy.square(eye_cms),axis=2)))
    #make sure theta < 1 
    theta_cos_ratio[theta_cos_ratio > 1] = 1
    #convert theta to angle
    abs_theta = numpy.arccos(theta_cos_ratio)
    #correct the angle to be positive if looking left and negative if looking right
    eye_HH[:,:,0] = -numpy.sign(eye_cms[:,:,0]) * abs_theta

    eye_HH = eye_HH*(180/math.pi)

    return eye_HH

def create_angular_disparity_mats(lefteye_pixels,lefteye_HH,righteye_pixels,righteye_HH,cycloeye_pixels,disparity_mat,vdisparity_mat, disparity_cyclo_mat, vdisparity_cyclo_mat):
    '''takes in look up tables of camera to each eye projections and creates disparity matrix for left eye
    '''

    #vectorize x coords and y coords for each eye
    lefteye_xs_all = numpy.reshape(lefteye_pixels[:,:,0],lefteye_pixels[:,:,0].size)
    lefteye_ys_all = numpy.reshape(lefteye_pixels[:,:,1],lefteye_pixels[:,:,1].size)
    righteye_xs_all = numpy.reshape(righteye_pixels[:,:,0],righteye_pixels[:,:,0].size)
    righteye_ys_all = numpy.reshape(righteye_pixels[:,:,1],righteye_pixels[:,:,1].size)
    cycloeye_xs_all = numpy.reshape(cycloeye_pixels[:,:,0],cycloeye_pixels[:,:,0].size)
    cycloeye_ys_all = numpy.reshape(cycloeye_pixels[:,:,1],cycloeye_pixels[:,:,1].size)

    #vectorize HH azimuth and elevation for each eye
    lefteye_xs_all_HH = numpy.reshape(lefteye_HH[:,:,0],lefteye_HH[:,:,0].size)
    lefteye_ys_all_HH = numpy.reshape(lefteye_HH[:,:,1],lefteye_HH[:,:,1].size)
    righteye_xs_all_HH = numpy.reshape(righteye_HH[:,:,0],righteye_HH[:,:,0].size)
    righteye_ys_all_HH = numpy.reshape(righteye_HH[:,:,1],righteye_HH[:,:,1].size)

    #make disparity vectors
    disparity_vec_all = lefteye_xs_all_HH - righteye_xs_all_HH
    vdisparity_vec_all = lefteye_ys_all_HH - righteye_ys_all_HH

    #set selection criterian, ie, eye sensor size, for left eye
    left_logic1 = numpy.logical_and(lefteye_xs_all > 0, lefteye_xs_all < 549)
    left_logic2 = numpy.logical_and(lefteye_ys_all > 0, lefteye_ys_all < 549)
    cyclo_logic1 = numpy.logical_and(cycloeye_xs_all > 0, cycloeye_xs_all < 549)
    cyclo_logic2 = numpy.logical_and(cycloeye_ys_all > 0, cycloeye_ys_all < 549)

    #apply criteria and convert to ints
    lefteye_xs = numpy.round(lefteye_xs_all[numpy.where(numpy.logical_and(left_logic1, left_logic2))])
    lefteye_ys = numpy.round(lefteye_ys_all[numpy.where(numpy.logical_and(left_logic1, left_logic2))])
    cycloeye_xs = numpy.round(cycloeye_xs_all[numpy.where(numpy.logical_and(cyclo_logic1, cyclo_logic2))])
    cycloeye_ys = numpy.round(cycloeye_ys_all[numpy.where(numpy.logical_and(cyclo_logic1, cyclo_logic2))])
    disparity_vec = disparity_vec_all[numpy.where(numpy.logical_and(left_logic1, left_logic2))]
    vdisparity_vec = vdisparity_vec_all[numpy.where(numpy.logical_and(left_logic1, left_logic2))]

    disparity_cyclo_vec = disparity_vec_all[numpy.where(numpy.logical_and(cyclo_logic1, cyclo_logic2))]
    vdisparity_cyclo_vec = vdisparity_vec_all[numpy.where(numpy.logical_and(cyclo_logic1, cyclo_logic2))]


    disparity_mat[-lefteye_ys.astype('int'),lefteye_xs.astype('int')] = disparity_vec
    vdisparity_mat[-lefteye_ys.astype('int'),lefteye_xs.astype('int')] = vdisparity_vec

    disparity_cyclo_mat[-cycloeye_ys.astype('int'),cycloeye_xs.astype('int')] = disparity_cyclo_vec
    vdisparity_cyclo_mat[-cycloeye_ys.astype('int'),cycloeye_xs.astype('int')] = vdisparity_cyclo_vec

    return disparity_mat,vdisparity_mat, disparity_cyclo_mat, vdisparity_cyclo_mat

def create_disparity_mats(lefteye_pixels,righteye_pixels,disparity_mat,vdisparity_mat):
    '''takes in look up tables of camera to each eye projections and creates disparity matrix for left eye
    '''

    #vectorize x coords and y coords for each eye
    lefteye_xs_all = numpy.reshape(lefteye_pixels[:,:,0],lefteye_pixels[:,:,0].size)
    lefteye_ys_all = numpy.reshape(lefteye_pixels[:,:,1],lefteye_pixels[:,:,1].size)
    righteye_xs_all = numpy.reshape(righteye_pixels[:,:,0],righteye_pixels[:,:,0].size)
    righteye_ys_all = numpy.reshape(righteye_pixels[:,:,1],righteye_pixels[:,:,1].size)

    #make disparity vectors
    disparity_vec_all = lefteye_xs_all - righteye_xs_all
    vdisparity_vec_all = lefteye_ys_all - righteye_ys_all

    #set selection criterian, ie, eye sensor size, for left eye
    left_logic1 = numpy.logical_and(lefteye_xs_all > 0, lefteye_xs_all < 549)
    left_logic2 = numpy.logical_and(lefteye_ys_all > 0, lefteye_ys_all < 549)

    #apply criteria and convert to ints
    lefteye_xs = numpy.round(lefteye_xs_all[numpy.where(numpy.logical_and(left_logic1, left_logic2))])
    lefteye_ys = numpy.round(lefteye_ys_all[numpy.where(numpy.logical_and(left_logic1, left_logic2))])
    disparity_vec = disparity_vec_all[numpy.where(numpy.logical_and(left_logic1, left_logic2))]
    vdisparity_vec = vdisparity_vec_all[numpy.where(numpy.logical_and(left_logic1, left_logic2))]

    disparity_mat[-lefteye_ys.astype('int'),lefteye_xs.astype('int')] = disparity_vec
    vdisparity_mat[-lefteye_ys.astype('int'),lefteye_xs.astype('int')] = vdisparity_vec

    return disparity_mat,vdisparity_mat

def draw_fixation_point(image,fixation_pixels,rgb = (255,255,0)):
    '''add projection of 3d fixation point in to image as a cross
    '''

    image[int(round(-fixation_pixels[1]))-8:int(round(-fixation_pixels[1]))+8,int(round(fixation_pixels[0]))-1:int(round(fixation_pixels[0]))+2,0] = rgb[2]
    image[int(round(-fixation_pixels[1]))-8:int(round(-fixation_pixels[1]))+8,int(round(fixation_pixels[0]))-1:int(round(fixation_pixels[0]))+2,1] = rgb[1]
    image[int(round(-fixation_pixels[1]))-8:int(round(-fixation_pixels[1]))+8,int(round(fixation_pixels[0]))-1:int(round(fixation_pixels[0]))+2,2] = rgb[0]
    image[int(round(-fixation_pixels[1]))-1:int(round(-fixation_pixels[1]))+2,int(round(fixation_pixels[0]))-8:int(round(fixation_pixels[0]))+8,0] = rgb[2]
    image[int(round(-fixation_pixels[1]))-1:int(round(-fixation_pixels[1]))+2,int(round(fixation_pixels[0]))-8:int(round(fixation_pixels[0]))+8,1] = rgb[1]
    image[int(round(-fixation_pixels[1]))-1:int(round(-fixation_pixels[1]))+2,int(round(fixation_pixels[0]))-8:int(round(fixation_pixels[0]))+8,2] = rgb[0]

    return image

def make_colormap_disparity(disparity_image,disparity_mat):
    '''takes in disparity matrix and make color map of crossed and uncrossed disparities
    '''

    shift_vals = calc_shift_vals(disparity_mat)

    disparity_mat = disparity_mat.copy() - shift_vals

    if numpy.any(numpy.isnan(disparity_mat) == 0):

        #set max and min for crossed and uncross disparities
        cmax = numpy.percentile(disparity_mat[numpy.isnan(disparity_mat) == 0],99)
        cmin = numpy.percentile(disparity_mat[numpy.isnan(disparity_mat) == 0],1)
        #set invalid pixels to black
        disparity_image[numpy.isnan(disparity_mat),:] = numpy.array((0,0,0))
        #set crossed disparities to reds
        disparity_image[numpy.logical_and(disparity_mat < 0,numpy.isnan(disparity_mat) == 0),0] = 180 - (180*(disparity_mat[numpy.logical_and(disparity_mat < 0,numpy.isnan(disparity_mat) == 0)]/cmin))
        disparity_image[numpy.logical_and(disparity_mat < 0,numpy.isnan(disparity_mat) == 0),1] = 180 - (180*(disparity_mat[numpy.logical_and(disparity_mat < 0,numpy.isnan(disparity_mat) == 0)]/cmin))
        disparity_image[numpy.logical_and(disparity_mat < 0,numpy.isnan(disparity_mat) == 0),2] = 255
        #set uncrossed disparities to blues
        disparity_image[numpy.logical_and(disparity_mat > 0,numpy.isnan(disparity_mat) == 0),0] = 255
        disparity_image[numpy.logical_and(disparity_mat > 0,numpy.isnan(disparity_mat) == 0),1] = 180 - (180*(disparity_mat[numpy.logical_and(disparity_mat > 0,numpy.isnan(disparity_mat) == 0)]/cmax))
        disparity_image[numpy.logical_and(disparity_mat > 0,numpy.isnan(disparity_mat) == 0),2] = 180 - (180*(disparity_mat[numpy.logical_and(disparity_mat > 0,numpy.isnan(disparity_mat) == 0)]/cmax))
        #set zero disparity to yellow
        disparity_image[numpy.logical_and(disparity_mat == 0,numpy.isnan(disparity_mat) == 0),0] = 0
        disparity_image[numpy.logical_and(disparity_mat == 0,numpy.isnan(disparity_mat) == 0),1] = 255
        disparity_image[numpy.logical_and(disparity_mat == 0,numpy.isnan(disparity_mat) == 0),2] = 255


    return disparity_image

def calc_shift_vals(arr):

    shape = arr.shape
    center_row = shape[0]/2
    center_col = shape[1]/2
    box_size=1

    central_box = arr[center_row-box_size:center_row+box_size+1, center_col-box_size:center_col+box_size+1]
    if numpy.isnan(central_box).all():
        shift_vals = numpy.nan
    else:
        masked_center = numpy.ma.masked_invalid(central_box)
        center_avg = numpy.ma.mean(masked_center)
        # center_avg = numpy.ma.mean(center_avg, axis=0)
        # center_avg = center_avg.filled(numpy.nan)
        center = arr[center_row,center_col]
        shift_vals = numpy.where(numpy.isnan(center), center_avg, center)

    return shift_vals

def apply_circle_mask(img,circle_center):
    '''make 10 deg circle mask for eye images and disparity maps
    '''

    radius = 103. #10deg pixel radius with 583 fl

    start_pixel = int(circle_center-radius)
    end_pixel = int(circle_center+radius+1)

    #crop image to radius
    if len(img.shape) == 3:
        img = img[start_pixel:end_pixel,start_pixel:end_pixel,:]
    else:
        img = img[start_pixel:end_pixel,start_pixel:end_pixel]

    x,y = numpy.meshgrid(range(0,(2*int(radius))+1),range(0,(2*int(radius))+1))
    mask = (  (((x-radius)/(radius))**2 + ((y-radius)/(radius))**2) < 1 )
    img[mask == False] = numpy.nan

    return img.astype('float32') #convert to float32 here?

def find_eyes_Rt(fixation_pt, IPD):
    """Returns R_left, t_left, R_right, t_right, the rotation matrices and translation vectors
    for the left and right eyes given a fixation point in cyclopean eye coordinates and the IPD.

    fixation_pt can be 3x1 or 3x2 if the left and right eye's fixation points are different.
    In that case, left eye should go to the first column and right eye to the second.

    Units must match in all inputs! (e.g. both in cm)"""

    assert(fixation_pt.shape[1] < 3)

    t_left = numpy.matrix([[IPD/2.0],
                        [0],
                        [0]])

    t_right = numpy.matrix([[-IPD/2.0],
                        [0],
                        [0]])

    t_cyclopean = numpy.matrix([[0.],
                        [0],
                        [0]])

    # Format the fixation_pt input
    if fixation_pt.shape[1] != 2:
        fixation_pt = numpy.matrix(numpy.concatenate((fixation_pt, fixation_pt), axis=1))

    # Basically gluLookAt for a left handed coord system
    f_left = fixation_pt[:,0] + t_left
    f_left = f_left / numpy.linalg.norm(f_left)  # normalize
    f_right = fixation_pt[:,1] + t_right
    f_right = f_right / numpy.linalg.norm(f_right)  # normalize

    # add cyclopean eye
    f_cyclopean = fixation_pt[:,1] + t_right
    f_cyclopean = f_cyclopean / numpy.linalg.norm(f_cyclopean)

    up = numpy.matrix([[0], [1.0], [0]])  # up is positive y-axis

    R = list()
    for f in [f_left, f_right, f_cyclopean]:
        s = numpy.cross(up.squeeze(), f.squeeze())
        u = numpy.cross(f.squeeze(), s.squeeze())
        R.append(numpy.concatenate((s, u, f.T)))


    return R[0], t_left, R[1], t_right, R[2], t_cyclopean

def rotate_points(points_wc, R, t):
    """Returns the 3D coordinates of points in a new coordinate system given R and t from
    find_eyes_Rt(...)

    points_wc should be of shape 3xN and should have the same units as IPD and fixation_pt
    used to make R and t in find_eyes_Rt(...)"""

# R and t must be numpy.matrix types
    assert(type(R) == numpy.matrix)
    assert(type(t) == numpy.matrix)

    points_ec = R*(points_wc + t)

    return points_ec

def R_from_Helmholtz(phi, theta, psi=None):
    """Returns a rotation matrix R from Helmholtz coordinates phi, theta, and (optionally) psi.

    Phi represents the angle around the interocular axis and is positive for DOWNWARD movements.

    Theta represents the angle around the vertical axis and is positive for LEFTWARD movements.

    Psi optionally represents the torsion angle and is positive for clockwise torsion. If it
    is not explicitly set it is calculated from phi and theta based on Listing's Law.

    The returned matrix assumes a left handed coordinate system in which the x-axis is the
    interocular axis and increases to the right, the y-axis is the vertical axis and increases
    up, and the z-axis increases away from the viewer.
    """

    R = numpy.zeros((3,3))
    if psi is None:
        psi = -theta*phi / 2.0

    R[0,0] = numpy.cos(psi) * numpy.cos(theta)
    R[0,1] = -numpy.sin(psi) * numpy.cos(theta)
    R[0,2] = numpy.sin(theta)

    R[1,0] = -numpy.sin(phi)*numpy.sin(theta)*numpy.cos(psi) + numpy.sin(psi)*numpy.cos(phi)
    R[1,1] = numpy.sin(phi)*numpy.sin(psi)*numpy.sin(theta) + numpy.cos(phi)*numpy.cos(psi)
    R[1,2] = numpy.sin(phi)*numpy.cos(theta)

    R[2,0] = -numpy.sin(phi)*numpy.sin(psi) - numpy.sin(theta)*numpy.cos(phi)*numpy.cos(psi)
    R[2,1] = -numpy.sin(phi)*numpy.cos(psi) + numpy.sin(psi)*numpy.sin(theta)*numpy.cos(phi)
    R[2,2] = numpy.cos(phi)*numpy.cos(theta)

    R_phi = numpy.matrix([[1, 0, 0],
                       [0, numpy.cos(phi), -numpy.sin(phi)],
                       [0, numpy.sin(phi), numpy.cos(phi)]], dtype='double')
    R_theta = numpy.matrix([[numpy.cos(theta), 0, -numpy.sin(theta)],
                         [0, 1, 0],
                         [numpy.sin(theta), 0, numpy.cos(theta)]], dtype='double')
    R_psi = numpy.matrix([[numpy.cos(psi), numpy.sin(psi), 0],
                       [-numpy.sin(psi), numpy.cos(psi), 0],
                       [0, 0, 1]], dtype='double')

    R = R_phi*R_theta*R_psi

    return R

def extend_torsion(theta_L, theta_R, phi, G=.8):
    """Return the torsion angles for the left and right eye, psi_L and psi_R, using Listing's
    extended law with a gain, G, defaulting to .8.

    theta_L the Helmholtz horizontal coordinate for the left eye.
    theta_R the Helmholtz horizontal coordinate for the right eye.
    phi: the Helmholtz vertical component for both eyes.
    """

    D = theta_R - theta_L
    V = phi
    H_L = theta_L
    H_R = theta_R
    psi_L = -G*D*V/4 - H_L*V/2
    psi_R = G*D*V/4 - H_R*V/2

    return psi_L, psi_R

def Helmholtz_from_R(R):
    """Returns the Helmholtz coordinates for an arbitrary rotation matrix R.

    R should represent a rotation for a left-handed coordinate system in which the x-axis is
    the inter-ocular axis and increases to the right, the y-axis is the vertical axis and increases
    to upward, and the z-axis is the viewing axis and increases away from the viewer.

    The returned Helmholtz coordinates are phi, theta, and psi. Positive phi is a downward rotation
    about the interocular axis, positive theta is a leftward rotation about the vertical axis
    and positive psi is a clockwise torsional rotation. All angles returned in radians.
    """

    theta = numpy.arcsin(R[0,2])
    phi = numpy.arcsin(R[1,2] / numpy.cos(theta))
    psi = -numpy.arcsin(R[0,1] / numpy.cos(theta))

    return phi, theta, psi

def fixation_to_Helmholtz(fixation_points, IPD, G=0.8):
    """Takes in a 3x2 array the specifies the 3D coordinates of the left and right eye's
    fixation points respectively in cyclopean coordinates and outputs each eye's
    Helmholtz coordinates for that point.
    """

    fixation_points[0,0] += IPD/2.0
    fixation_points[0,1] -= IPD/2.0
    phi = []
    theta = []

    z_dir = numpy.array([0, 0, 1], dtype='double')

    for point in fixation_points.T:

        # Get the component of the fixation point that lies in the zy-plane
        zy_component = numpy.array(point)
        #zy_component = zy_component.squeeze()
        zy_component[0] = 0
        # Calculate the angle between the zy-component and the z-direction
        phi_cos_ratio = numpy.dot(z_dir, zy_component) / (numpy.linalg.norm(z_dir) * numpy.linalg.norm(zy_component))
        if phi_cos_ratio > 1:
            print "Phi cos ratio is", phi_cos_ratio
        abs_phi = numpy.arccos(numpy.min([1, phi_cos_ratio]))
        # Correct the angle to be positive if looking down and negative if looking up
        phi_norm = numpy.cross(z_dir, zy_component)
        if phi_norm[0]:
            phi_dir = phi_norm[0] / numpy.abs(phi_norm[0])
        else:
            phi_dir = 1
        phi.append(phi_dir * abs_phi)

        # Calculate the angle between the zy-component and the fixation vector
        theta_cos_ratio = numpy.dot(zy_component.T, point) / (numpy.linalg.norm(zy_component) * numpy.linalg.norm(point))
        if theta_cos_ratio > 1:
            print "Theta cos ratio is", theta_cos_ratio
        abs_theta = numpy.arccos(numpy.min([1, theta_cos_ratio]))
        # Correct the angle to positive if looking left and negative if looking right
        theta_norm = numpy.cross(point, zy_component)
        if theta_norm[1]:
            theta_dir = theta_norm[1] / numpy.abs(theta_norm[1])
        else:
            theta_dir = 1
        theta.append(theta_dir * abs_theta)

    assert numpy.round(phi[0],decimals=10) == numpy.round(phi[1],decimals=10), ('Elevation angles (phi) should come out equal', phi, fixation_points)

    # Calcualte torsion
    psi_L, psi_R = extend_torsion(theta[0], theta[1], phi[0], G)
    psi = [psi_L, psi_R]

    return numpy.array([phi, theta, psi])

# def fixation_to_Helmholtz_cyclopean(point, G=0.8):

#     # Convert cyclopean fixation point to Helmholtz coordinates, using same proceedure as fixation_to_Helmholtz
#     phi = []
#     theta = []

#     z_dir = numpy.array([0, 0, 1], dtype='double')

#     # for point in fixation_points.T:

#     # Get the component of the fixation point that lies in the zy-plane
#     zy_component = numpy.array(point)
#     #zy_component = zy_component.squeeze()
#     zy_component[0] = 0
#     # Calculate the angle between the zy-component and the z-direction
#     phi_cos_ratio = numpy.dot(z_dir, zy_component) / (numpy.linalg.norm(z_dir) * numpy.linalg.norm(zy_component))
#     if phi_cos_ratio > 1:
#         print "Phi cos ratio is", phi_cos_ratio
#     abs_phi = numpy.arccos(numpy.min([1, phi_cos_ratio]))
#     # Correct the angle to be positive if looking down and negative if looking up
#     phi_norm = numpy.cross(z_dir, zy_component)
#     if phi_norm[0]:
#         phi_dir = phi_norm[0] / numpy.abs(phi_norm[0])
#     else:
#         phi_dir = 1
#     phi.append(phi_dir * abs_phi)

#     # Calculate the angle between the zy-component and the fixation vector
#     theta_cos_ratio = numpy.dot(zy_component.T, point) / (numpy.linalg.norm(zy_component) * numpy.linalg.norm(point))
#     if theta_cos_ratio > 1:
#         print "Theta cos ratio is", theta_cos_ratio
#     abs_theta = numpy.arccos(numpy.min([1, theta_cos_ratio]))
#     # Correct the angle to positive if looking left and negative if looking right
#     theta_norm = numpy.cross(point, zy_component)
#     if theta_norm[1]:
#         theta_dir = theta_norm[1] / numpy.abs(theta_norm[1])
#     else:
#         theta_dir = 1
#     theta.append(theta_dir * abs_theta)

#     # assert numpy.round(phi[0],decimals=10) == numpy.round(phi[1],decimals=10), ('Elevation angles (phi) should come out equal', phi, fixation_points)

#     # Calcualte torsion
#     psi_L, psi_R = extend_torsion(theta[0], theta[1], phi[0], G)
#     psi = [psi_L, psi_R]

#     return numpy.array([phi, theta, psi])
