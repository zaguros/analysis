"""
A collection of functions to manipulate and analyze scans or camera images
Requires a working installation of opencv for python (a.k.a. cv2)
NK 2017
"""

import numpy as np
import cv2
from matplotlib import pyplot as plt
from itertools import compress
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common


##############################
# general image manipulation #
##############################



### apply a brightness filter we only want to detect bright features.
def apply_brightness_filter(array2d,threshold_min,threshold_max = 256):
    return array2d*(array2d>threshold_min)*(array2d<threshold_max)

### stamp out the central region of the image a.k.a. what is our field of view?
def stamp_out_relevant_field_of_view(array2d,xsize=300,ysize=300,xoffset = 10, yoffset = 60):
    """
    Slices an array according to given sizes.
    Input parameters with offset give an offset from the CENTRE of the array!
    """
    
    (xlen,ylen) =np.shape(array2d)

    xmin = int(xlen/2.-xsize/2.)+xoffset
    xmax = int(xlen/2.+xsize/2.)+xoffset
    ymin = int(ylen/2.-ysize/2.)+yoffset
    ymax = int(ylen/2.+ysize/2.)+yoffset

    return array2d[xmin:xmax,ymin:ymax] 

def make_brightness_hist(array2d, rng = (0,255)):
    """
    Makes a histogram according to the grayscale values of an array.
    assumes that the value range of the array is 0 to 255
    Excludes 0!
    """
    
    n,bins = np.histogram(array2d,bins=255,range=rng)
    fig = plt.figure()
    ax = plt.subplot()
    x = (bins[1:]-(bins[1]-bins[0])/2.)[1:] ### exclude 0
    plt.plot(x,n[1:]) ### exclude 0
    ax.set_xlabel('Pixel brightness')
    ax.set_ylabel('Occurences #')
    
    return x,n[1:]

def show_image(array2d,size = 7,no_col_bar = False):
    fig = plt.figure(figsize = (size,size))
    ax = plt.subplot()
    im = plt.imshow(array2d,cmap='gray',interpolation='none')
    if not no_col_bar:
        fig.colorbar(im) 
    plt.show()
    plt.close('all')

def make_binary(array2d,threshold_min):
    return (array2d>threshold_min)*np.ones(np.shape(array2d))


def rescale_to_8_bit_grayscale(img):
    """
    scales an image from 0 to 255.
    """
    img = np.abs(img)-(np.amin(np.abs(img))) ## set min to zero and scale to 255
    img = img*255/np.amax(img)
    return img.astype(np.uint8)

def rotate_image(img,angle):
    """
    rotates a given image by a certain angle. Watch out: applies interpolation and might smear out features.
    can also be useful for automated scanning. when movements of the piezos are expressed in lithography coordinates
    """
    image_center = tuple(np.array(img.shape)/2)
    print image_center
    rot_mat = cv2.getRotationMatrix2D(image_center,angle,1.0)
    print rot_mat
    result = cv2.warpAffine(img, rot_mat, img.shape)#,flags=cv2.INTER_LINEAR)
    return result



#####################################
# functions for grid identification #
#####################################

def detect_1D_clusters(array1d,distance = 13,verbose = False):
    """
    finds clustering in a 1d array. returns points within the clusters.
    This function is mainly used for grid detection and assumes a few parameters 
    Only clusters of a length of more than 6 values are taken in.
    Clusters are identified by differences in value of more than 10 (kwarg distance). 
    Distance is given in units of array indices of the original image.
    """
    x_cluster_distance = np.round(np.abs(np.diff(np.sort(array1d))))    
    ### try to identify chains of markers with a length > 6
    x_cluster_diff = np.diff(np.array(range(len(x_cluster_distance)))[x_cluster_distance>distance]) 
    
    clusters_to_be_fitted = x_cluster_diff>=6
    cluster_indices = np.cumsum(x_cluster_diff)[clusters_to_be_fitted]
    if verbose:
        print x_cluster_distance
        print x_cluster_diff
        print cluster_indices
    xcentres = np.sort(array1d)[cluster_indices] ### get some indication for the desired x value where we start the fit.
    
    return xcentres
    
def fit_1D_lines_to_clusters(x,y,cluster_positions):
    """
    Input: x coordinates of all keypoints
            y coordinates of all keypoints
            cluster positions along y
            Note that x and y might be interchanged to fit vertical or horizontal lines in the original image
    Output:
        fit_offsets; list containing the offsets
        fit_slopes; list containing the slopes
        fit_results; list containing the fit dictionaries
    """
    fit_offsets = [] 
    fit_slopes = []
    fit_results = []
    for x_centre in cluster_positions:
        x_filt = (x > x_centre-10) & (x<x_centre+10) ## get all relevant values
        key_xs_filt = x[x_filt]
        key_ys_filt = y[x_filt]
        
        
        ### fit a line to the estimated keypoints we took a slice along x so y becomes the effective x-coordinate for this fit.
        p0,fitfunc,fitfunc_str = common.fit_line(np.average(key_xs_filt),0.0)
        fit_result = fit.fit1d(key_ys_filt,key_xs_filt,None,p0= p0,fitfunc=fitfunc,do_print=False,ret=True,fixed=[])
        
        fit_offsets.append(fit_result['params_dict']['a']);fit_slopes.append(fit_result['params_dict']['b'])
        fit_results.append(fit_result)
        
    return fit_offsets,fit_slopes,fit_results
def find_marker_grid_via_blob_finder(filtered_laplacian_img,plot_histogram = False):
    laplace_8bit = rescale_to_8_bit_grayscale(filtered_laplacian_img)
    laplace_8bit  = apply_brightness_filter(laplace_8bit,54)

    ####
    # blob detector params
    ####
    if plot_histogram:
        x,y = make_brightness_hist(laplace_8bit)


    blob_params = cv2.SimpleBlobDetector_Params()
    blob_params.filterByColor = False

    blob_params.minThreshold =55
    blob_params.maxThreshold = 256
    # Filter by area
    blob_params.filterByArea = True
    blob_params.minArea = 0
    blob_params.maxArea = 30
    #Filter by Circularity
    blob_params.filterByCircularity = True
    blob_params.minCircularity = 0.5

    #Filter by Converxity
    blob_params.filterByConvexity = True
    blob_params.minConvexity = 1.0

    # Filter by Inertia
    blob_params.filterByInertia = True
    blob_params.minInertiaRatio = 0.3


    detector = cv2.SimpleBlobDetector_create(blob_params)
    keypoints = detector.detect(laplace_8bit)
    print len(keypoints)
    im_with_keypoints = cv2.drawKeypoints(laplace_8bit,keypoints,np.array([]),
                                          (255,0,0),cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)

    return im_with_keypoints, keypoints,laplace_8bit

def estimate_grid_from_keypts(im,keypoints,plot_fitted_lines = False,verbose = False):
    """
    tries to fit a rectangular grid to the estimated keypoints
    currently only works for relatively small angular deviations modulo 90 degrees
    A grid at 45 degree rotation is the wrost case scenario
    returns angle, pitch and offset from the origin x0
    """
    ## extract all x and y values
    key_xs = []
    key_ys = []
    for k in keypoints:
        (x,y) = k.pt
        key_xs.append(x);key_ys.append(y)
    key_xs = np.array(key_xs); key_ys = np.array(key_ys)
    
    ### now take several x-slices and fit a line to the given keypoints
    
    ### we first cluster the keypoint_xcoords by sorting them and looking for the largest gaps. 
    ### (after all the markers are expected to be almost on a straightline
    xcentres = detect_1D_clusters(key_xs,verbose = verbose)
    ycentres = detect_1D_clusters(key_ys,verbose = verbose)
    x_offset,x_slope,x_fits = fit_1D_lines_to_clusters(key_xs,key_ys,xcentres)
    y_offset,y_slope,y_fits = fit_1D_lines_to_clusters(key_ys,key_xs,ycentres)
    
    if plot_fitted_lines:
        for xcent,fit_result in zip(xcentres,x_fits):
            ymin = np.amin(key_ys).astype(np.uint16)
            ymax = np.amax(key_ys).astype(np.uint16)
            cv2.line(im,(fit_result['fitfunc'](ymin).astype(np.uint16),ymin),
                               (fit_result['fitfunc'](ymax).astype(np.uint16),ymax),(255,0,0),1)
            
        for ycent,fit_result in zip(ycentres,y_fits):
            ymin = np.amin(key_xs).astype(np.uint16)
            ymax = np.amax(key_xs).astype(np.uint16)
            cv2.line(im,(ymin,fit_result['fitfunc'](ymin).astype(np.uint16)),
                               (ymax,fit_result['fitfunc'](ymax).astype(np.uint16)),(255,0,0),1)
    
    
    ### average over all the slopes:
    print "this is the average slope: ", np.average(np.append(x_slope,-1*y_slope))
    avg_slope = np.average(np.append(x_slope,-1*y_slope))
    grid_distances = np.append(np.diff(y_offset),np.diff(x_offset))
    avg_distance = np.average(grid_distances[grid_distances < 30]) ### one could do something smarter here...
    print "this is the average grid distance", avg_distance ## assume a typical distance smaller than 30 indices
    
    
    ### estimate origin by choosing the best fit (smallest error) and from there the point which is closest to the fitted grid.
    err_dict_offset = []
    err_dict_slope = []
    print len(x_fits),len(y_fits)
    all_fits = x_fits +y_fits
    for fit in all_fits:
        err_dict_offset.append(fit['error_dict']['a'])
        err_dict_slope.append(fit['error_dict']['b'])
    
    err_dict_offset = np.array(err_dict_offset);err_dict_slope = np.array(err_dict_slope)
    total_error = err_dict_offset + err_dict_slope # simple addition of the uncertainties to gain some fit fitness
    best_fit = all_fits[np.argmin(total_error)]
    best_fit_is_x_fit = False
    if np.argmin(total_error) < len(x_fits):
        best_fit_is_x_fit = True
    x0 = best_fit['x'][0]; y0 = best_fit['y'][0]
    
    return avg_slope,avg_distance,x0,y0

def get_grid_crossing_pts(x0,y0,avg_slope,avg_distance):
    grid_positions = np.array(range(-10,10,1))
    grid_angle = np.arctan(avg_slope)
    
    ## generate crossing points for 1 horizontal line
    x_pos1 = grid_positions*avg_distance*np.cos(grid_angle)+x0
    y_pos1 = grid_positions*avg_distance*np.sin(grid_angle)+y0
    xy_pos = []
    ### we now loop over the crossing points of this horizontal line and generate all crossing points over the entire image
    for x,y in zip(x_pos1,y_pos1):
        xs = -grid_positions*avg_distance*np.sin(grid_angle)+x
        ys = grid_positions*avg_distance*np.cos(grid_angle)+y
        xy = xs+1j*ys ### make this numbers imaginary for ease of handling.
        xy_pos.extend(xy)
        
    x_filt = (np.array(xy_pos).real > 0) & (np.array(xy_pos).real < 300)
    y_filt = (np.array(xy_pos).imag > 0) & (np.array(xy_pos).imag < 300)
    xy_filt = np.logical_and(x_filt,y_filt)
    xy_pos = np.array(xy_pos)[xy_filt]
    return xy_pos


def stamp_out_marker_grid_onto_im(im,all_marker_pts,marker_size,show_images = False):
    """
    Sets a specific window (marker_size) at a certain position (all_marker_pts) to 0.
    """
    marker_pos = np.zeros(np.shape(im))
    [indices_x,indices_y] = np.indices(np.shape(im))
    some_ones = np.ones(np.shape(im))
    for mx,my in zip(all_marker_pts.real,all_marker_pts.imag):
        marker_pos = marker_pos + some_ones*(indices_x < mx + marker_size)*(indices_x > mx - marker_size)*(indices_y > my - marker_size)*(indices_y < my + marker_size)
    if show_images:
        show_image(im)
        show_image(im*(marker_pos-1)*(-1))
    return im*(marker_pos-1)*(-1)


####################################
# Bit marker detection and readout # (a.k.a. where am i currently on the sample?)
####################################

def find_bit_marker_in_image(filtered_img,do_plot = False):
    """
    expects a edged version of the original image. then returns the coordinates of the bit marker in the image    
    edged versions of images can be obtained with e.g.
    cv2.Canny(ct.rescale_to_8_bit_grayscale(original_image),100,185)
    """    
    filtered_img_8bit = rescale_to_8_bit_grayscale(filtered_img)
    ideal_bitmarker_identifier = np.array([[[0,0]],[[0,15]],[[3,15]],[[3,0]]])

    thresh,filtered_img_8bit = cv2.threshold(filtered_img_8bit,80,255,0)


    edges = cv2.Canny(filtered_img_8bit,210,255)
    
    

    ideal_bitmarker_identifier = np.array([[[0,0]],[[0,15]],[[3,15]],[[3,0]]])
    
    im2,contours,hierarchy = cv2.findContours(edges, 1, cv2.CHAIN_APPROX_SIMPLE)

    areas = np.zeros(len(contours))
    cxs = areas.copy()
    cys =  areas.copy()
    c_filt =  np.array([False]*len(contours))
    shape_matching =  areas.copy()
    for ii,c in enumerate(contours):
        M = cv2.moments(c)
        if M['m00'] != 0: ## area equals 0??
            areas[ii] = M['m00']
            cxs[ii] = int(M['m10']/M['m00'])
            cys[ii] = int(M['m01']/M['m00'])
            c_filt[ii] = True
            epsilon = 0.01*cv2.arcLength(c,True)
            approx = cv2.approxPolyDP(c,epsilon,True)
            x,y,w,h = cv2.boundingRect(c)
            c_filt[ii] = c_filt[ii] and w>h*2 ## we only want elongated rectangles

            shape_matching[ii] = cv2.matchShapes(approx,ideal_bitmarker_identifier,1,0.0)

    shape_fltr = np.logical_and(np.logical_and(c_filt,areas>30),areas<100)
    # print shape_matching[shape_fltr]
    best_shape_fit = np.amin(shape_matching[shape_fltr])
    if best_shape_fit > 1.:
        print 'could not find bit marker in image'
        return False
    contour_index = np.argmin(shape_matching[shape_fltr])
    real_contour_index = np.where(np.cumsum(shape_fltr) == contour_index+1)[0]

    bit_marker_identifier = contours[real_contour_index[0]]

    #### now want to find the position of the identifier in the image
    cv2.drawContours(edges, [bit_marker_identifier], -1, (110, 0, 0), 3)
    if do_plot:
        show_image(edges,size = 7,no_col_bar = True)

    bit_y = np.amin(bit_marker_identifier.T[0]) + 10
    bit_x = np.amin(bit_marker_identifier.T[1]) + 2
    
    return bit_x,bit_y

def distance_bitm_and_laserspot(bit_x,bit_y,spot_x,spot_y):
    """
    Assumes a resolution of the camera image of 200 nm per pixel.
    returns the distance in microns between detected marker shape (bit marker) and a calibrated position on the camera (ideally this is the laser spot)
    """
    return int((bit_x-spot_x/2.)/5),int((bit_y-spot_y/2.)/5)


def zoom_in_on_bitmarker(img,bit_x,bit_y):
    """
    takes cam image and cuts out the bitmarker. Assumes a resoultion of 200 nm per pixel and that the image is prerotated.
    """
    xlen,ylen = np.shape(img)
    xoff = int(bit_x - xlen/2.)-16
    yoff = int(bit_y - ylen/2.)


    bit_marker_img = stamp_out_relevant_field_of_view(img,xsize = 25,ysize = 25,xoffset = xoff,yoffset = yoff)
    bit_marker_img = rescale_to_8_bit_grayscale(bit_marker_img)
    return make_binary(bit_marker_img,100)


def generate_marker_shape(img,xpos,ypos,size):
    """
    generates a rectangular series of 1s in an array that is the size of the input img (img)
    """
    xlen,ylen = np.shape(img)
    new_img = np.ones((xlen,ylen))
    xinds,yinds = np.indices((xlen,ylen))
    x_filt = np.logical_and(xinds < (xpos + size),xinds > (xpos-size))
    y_filt = np.logical_and(yinds < (ypos + size),yinds > (ypos-size))
    
    return new_img*x_filt*y_filt
    
def get_bit_marker_array(binary_bit_marker_img):
    """
    compares a given BINARY bit marker image with marker shapes at several different positions
    """
    bit_mrkr_string = np.zeros(16)
    xcoords,ycoords = np.where(binary_bit_marker_img == 1)
    x_start = np.amax(xcoords)
    y_start = np.amax(ycoords)
    marker_size = 2.5
    ii = 0
    for x_pos in range(4):
        for y_pos in range(4):
            marker_shape = generate_marker_shape(binary_bit_marker_img.copy(),x_start-marker_size*2*(0.5+x_pos),
                                                                 y_start-marker_size*2*(0.5+y_pos),marker_size)
            bit_mrkr_string[ii] = np.sum(marker_shape*binary_bit_marker_img)
            ii+=1
    return np.round(bit_mrkr_string/np.amax(bit_mrkr_string))

def get_bitm_xy_from_array(bitm_array):
    """converts a bitmarker pattern (that is ideally received from camera) into x,y relative position of the bitm"""
    bity = bitm_array[-2,0]*2**0 + bitm_array[-3,0]*2**1 + bitm_array[-4,0]*2**2 + bitm_array[-4,1]*2**3 + bitm_array[-4,2]*2**4 + bitm_array[-3,1]*2**5 
    bitx = bitm_array[-1,1]*2**0 + bitm_array[-1,2]*2**1 + bitm_array[-1,3]*2**2 + bitm_array[-2,3]*2**3 + bitm_array[-3,3]*2**4 + bitm_array[-2,2]*2**5 
    return bitx,bity


######
# drawing & manioulating the ebeam pattern
######

def generate_bitm_array(bit_x,bit_y):
    bitx_string = '{0:06b}'.format(bit_x)
    bity_string = '{0:06b}'.format(bit_y)

    return np.array([[bity_string[-3],bity_string[-4],bity_string[-5],1],[bity_string[-2],bity_string[-6],1,bitx_string[-5]],
                            [bity_string[-1],1,bitx_string[-6],bitx_string[-4]],[1,bitx_string[-1],bitx_string[-2],bitx_string[-3]],
                            [0,0,0,0],
                            [1,1,1,1]],dtype=int)


def generate_marker_pattern(small_marker_pitch=5,bit_x_max = 36,bit_y_max = 15,pitch_bitm = 60):
    """
    generates the marker pattern from a few simple codewords
    Assumes an array resolution of 1 um.
    """
    x_dim = bit_x_max*pitch_bitm + 6 ## 6 is the size of a bitmarker in micron.
    y_dim = bit_y_max*pitch_bitm + 6 
    full_img = np.ones((y_dim,x_dim))
    x_inds,y_inds = np.indices((y_dim,x_dim))
    marker_image = full_img*(np.mod(np.flipud(x_inds)-2,small_marker_pitch) == 0) * (np.mod(y_inds,small_marker_pitch) == 0)
    
    current_bit_m_number = 0
    for y in range(bit_y_max):
        for x in range(bit_x_max):
            x_start = x*pitch_bitm
            y_start = y_dim - y*pitch_bitm
            
            marker_image[y_start-6:y_start,x_start:x_start+4] = generate_bitm_array(x,y)
            current_bit_m_number +=1
            
    return marker_image

def add_striplines_to_img(img,stripline_width = 30,stripline_centre = [240+4,244+474]):
    """always assumes that the stripline is horizontal."""
    x_ind,y_ind = np.indices(np.shape(img))
    for c in stripline_centre:
        line = np.ones(np.shape(img))*(x_ind > c-stripline_width/2.)*(x_ind < c+stripline_width/2.)
        img = np.logical_or(img,np.flipud(line))
    return img*np.ones(np.shape(img))


def correct_faulty_lift_off(bitm_x,bitm_y,bitm_pitch,x_dim,y_dim):
    """
    some times bitmarkers might not make sense! Then we need to correct the detected bitmarker pattern
    """
    if bitm_x*bitm_pitch > x_dim:
        bitm_x = bitm_x - 32 
    if bitm_y*bitm_pitch > y_dim : ### note that the output image will always be 90 degrees rotated....
        bitm_y = bitm_y - 32
    return bitm_x,bitm_y

def pattern_zoom_on_bitm(ebeam_pattern,bitm_x,bitm_y,bitm_pitch = 60,rel_size = 1,rel_shift_x = 0,rel_shift_y = 0,correct_lift_off_problems = True):
    """ returns a cut out of the size of the bitm_pitch """


    rel_shift_y = -rel_shift_y

    y_dim,x_dim = np.shape(ebeam_pattern)
    if correct_lift_off_problems:
        bitm_x,bitm_y = correct_faulty_lift_off(bitm_x,bitm_y,bitm_pitch,x_dim,y_dim)

    if bitm_y == 0: ## got to shift the image because the lowest row is 
        bitm_y = 1
    if bitm_x == 0: ## got to shift the image because the lowest row is 
        bitm_x = 1

    ymin = y_dim-int((0.5*rel_size+bitm_y)*bitm_pitch) + rel_shift_y
    ymax = y_dim+int((1/2.*rel_size-bitm_y)*bitm_pitch) + rel_shift_y
    xmin = int((bitm_x-0.5*rel_size)*bitm_pitch) + rel_shift_x
    xmax = int((1/2.*rel_size+bitm_x)*bitm_pitch) + rel_shift_x
    return ebeam_pattern[ymin:ymax,xmin:xmax] ### note that arrays are flipped by 90 degrees when plotting with imshow

def draw_spot_onto_pattern(pattern,brightness,bitm_x,bitm_y,rel_x,rel_y,bitm_pitch = 60,correct_lift_off_problems=True):
    """
    does this need a resolution? should be handled by a different function that forwards indices. We assume a micron
    """
    
    x,y = find_current_ebeam_coordinates(pattern,bitm_x,bitm_y,rel_x,rel_y,bitm_pitch = bitm_pitch,correct_lift_off_problems=correct_lift_off_problems)
    pattern[y,x] = brightness
    return pattern

def find_current_ebeam_coordinates(pattern,bitm_x,bitm_y,rel_x,rel_y,bitm_pitch = 60,correct_lift_off_problems=True):
    rel_y = -rel_y
    y_dim,x_dim = np.shape(pattern)

    if correct_lift_off_problems:
        bitm_x,bitm_y = correct_faulty_lift_off(bitm_x,bitm_y,bitm_pitch,x_dim,y_dim)

    y = int(y_dim-bitm_y*bitm_pitch + rel_y)
    x = int(bitm_x*bitm_pitch + rel_x)
    
    return x,y