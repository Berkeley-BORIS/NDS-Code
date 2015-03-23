"""
This script plots histograms of the disparity data, either as a 2D plot of
eccentrity vs disparity or as a 1D histogram of disparity collapsed
across eccentricities.
"""
import os
import sys

import matplotlib.pyplot as plt
from matplotlib import cm, colors
import scipy.io as spio
import numpy as np
import pandas
from scipy.interpolate import interp1d
import weighted_percentile


EXPECTED_SHAPE = (207, 207)
MERIDIAN_PIXEL = 103
FOCAL_LENGTH = 583
WEIGHT_DICT = {'sandwich': .29, 'ordering_coffee': .51, 'nature_walk_1':.20}

def plot_disp_hexbin(disp_data, grouping_name=()):
    """
    Plots a 2D hexbin of the disparity data with eccentricity on the abscissa
    and disparity on the ordinate.
    """
    xvals = []
    disparities = []

    disp_range = (-.75, .75)
    ecc_range = (-5, 5)

    ndisp_data = disp_data - disp_data.ix[MERIDIAN_PIXEL]

    # Build up the list of eccentricity/disparity pairs for hexbin
    for row, sers in ndisp_data.iterrows():
        nona_sers = sers.dropna()
        temp_x = np.r_[[row]*len(nona_sers)]
        xvals.append(temp_x)
        disparities.append(nona_sers.values)

    xvals = np.concatenate(xvals)
    disparities = np.concatenate(disparities)

    fig, ax = plt.subplots(1,1)
    mincnt = 50 #min(25, ndisp_data.shape[1] / 1000)
    # NOTE The default number of bins/bin sizes causes aliasing with our
    # discrete data points. The gridsize should be an even divisor of the length
    # of the data (ie 2, 4, 6...)
    hb = ax.hexbin(xvals, disparities, gridsize=disp_data.shape[0]/2, mincnt=mincnt,
                   extent=(-10, 10, -2, 2), bins='log')
    plt.colorbar(hb)
    medians = ndisp_data.median(axis=1)
    line, = ax.plot(ndisp_data.index.tolist(), medians.values)
    ax.set_xlim(ecc_range)
    ax.set_ylim(disp_range)

    ax.set_xlabel('Eccentricity (degrees)')
    ax.set_ylabel('Disparity (degrees)')

    return fig

def plot_disp_image(name_str, meridian, data, grouping_name=()):
    disp_range = (-.75, .75)
    ecc_range = (-5, 5)

    bins_ecc = get_gridded_binedges(data, samples_per_bin=2)
    bins_disparity = np.linspace(-2, 2, num=(len(bins_ecc)*2))
    bins_disparity = np.r_[-1e10, bins_disparity, 1e10]

    normed_data = data-data.ix[MERIDIAN_PIXEL]
    normed_data[MERIDIAN_PIXEL] = np.nan

    args = normed_data, (bins_ecc, bins_disparity)
    if grouping_name == 'all':
        H, xedges, yedges, medians = calc_weighted_2d_hist(*args)

    else:
        H, xedges, yedges = calc_2d_hist(*args)
        medians = normed_data.median(axis=1)

    fig, ax = plt.subplots(1,1)

    im = ax.imshow(np.log(H.T[1:-1, :]), extent=[xedges[0], xedges[-1], yedges[1], yedges[-2]],
                   origin='lower', interpolation='none', aspect='auto',
                   cmap=cm.autumn, vmin=-9.5, vmax=-4.5)
    
    #im = ax.imshow(H.T[1:-1, :], extent=[xedges[0], xedges[-1], yedges[1], yedges[-2]],
    #                              origin='lower', interpolation='none', aspect='auto',
    #                              cmap=cm.autumn, vmin=.0001, vmax=.008)
                                                 
    im.cmap.set_under('w')
    cb = plt.colorbar(im)

    #line, = ax.plot(normed_data.index.tolist(), medians.values, lw=2, color='b')
    
    #plot horopter    
    hor = spio.loadmat('horopter_medians_and_cis.mat')
    
    if meridian[0:2] == 'VM':
        line2, = ax.plot(hor['vh_ecc'].ravel(),-hor['vh_med'].ravel(), lw=2,color='k')
        line5, = ax.plot(hor['vh_ecc'].ravel(),-hor['vh_cis'][0].ravel(), 'k--', lw=2)
        line6, = ax.plot(hor['vh_ecc'].ravel(),-hor['vh_cis'][1].ravel(), 'k--', lw=2)
    else:
        line3, = ax.plot(hor['hh_ecc'].ravel(),-hor['hh_med'].ravel(), lw=2,color='g')
        line4, = ax.plot(hor['hh_ecc'].ravel(),-hor['hh_cis'][0].ravel(), 'g--',lw=2)
        line7, = ax.plot(hor['hh_ecc'].ravel(),-hor['hh_cis'][1].ravel(), 'g--',lw=2)
    
    #plot weighted medians
    if name_str[4:] == 'corrected' or name_str[4:] == 'random':
        print 'Loading weighted medians', name_str
        vm = spio.loadmat('../VM_dispdata_' + name_str[4:] + '_weighted_median_and_ci.mat')
        hm = spio.loadmat('../HM_dispdata_' + name_str[4:] + '_weighted_median_and_ci.mat')
    
        if meridian[0:2] == 'VM':
            line8, = ax.plot(normed_data.index.tolist(),vm['weighted_medians'].ravel(), lw=2,color='m')
            line9, = ax.plot(normed_data.index.tolist(),vm['cis'][0].ravel(), 'm--', lw=2)
            line10, = ax.plot(normed_data.index.tolist(),vm['cis'][1].ravel(), 'm--', lw=2)
            
            chi_square(hor['vh_ecc'],-hor['vh_med'],-hor['vh_cis'],normed_data.index.tolist(),vm['weighted_medians'],vm['cis'])
        else:
            line11, = ax.plot(normed_data.index.tolist(),hm['weighted_medians'].ravel(), lw=2,color='c')
            line12, = ax.plot(normed_data.index.tolist(),hm['cis'][0].ravel(), 'c--',lw=2)
            line13, = ax.plot(normed_data.index.tolist(),hm['cis'][1].ravel(), 'c--',lw=2)
            
            chi_square(hor['hh_ecc'],-hor['hh_med'],-hor['hh_cis'],normed_data.index.tolist(),hm['weighted_medians'],hm['cis'])
    
    
    #import pdb; pdb.set_trace()
    
    ax.set_xlabel('Eccentricity (degrees)')
    ax.set_ylabel('Disparity (degrees)')
    ax.set_ylim(disp_range)
    ax.set_xlim(ecc_range)

    return fig

def get_gridded_binedges(data, samples_per_bin):

    # TODO These bins aren't exactly centered, so when plotted they are a little
    # off. Note that this makes no difference to the analysis. It's
    # purely aesthetic.
    row_inds = data.index.tolist()
    row_inds.sort()
    row_inds = np.array(row_inds)

    shift_amount = .02
    row_inds -= shift_amount
    bins = np.r_[row_inds, row_inds[-1] + shift_amount]

    return bins[::samples_per_bin]

def calc_weighted_2d_hist(data, bins):
    
    groups = data.groupby(level='task_id', axis=1)
    H = {}
    medians = {}
    varis = {}
    mat_dict = {}
    for i, group in groups:
        print 'i:', i
        task_id = group.columns[0][2]
        H[task_id], xedges, yedges = calc_2d_hist(group, bins)
        medians[task_id] = group.median(axis=1)
        
        # build mat dict for saving
        if task_id:
            mat_dict[task_id] = group.values
        
    
    fixation_id = group.columns[0][1]
    #print 'Saving mat', fixation_id
    #spio.savemat('dispdata_' + fixation_id + '.mat', mat_dict)


    H_weighted = np.zeros(H[task_id].shape)
    medians_weighted = pandas.Series(np.zeros(len(medians[task_id])), index=medians[task_id].index)
    for task in WEIGHT_DICT:
        H_weighted += H[task] * WEIGHT_DICT[task]
        medians_weighted += medians[task] * WEIGHT_DICT[task]

    return H_weighted, xedges, yedges, medians_weighted

def calc_bootstrap_and_median(data):
    
    groups = data.groupby(level='task_id', axis=1)
    H = {}
    medians = {}
    vars = {}
    for i, group in groups:
        task_id = group.columns[0][2]
        for s in xrange:
            medians[task_id] = group.median(axis=1)
    
    return medians_weighted

def calc_2d_hist(data, bins):
    xvals = []
    disparities = []
    # Build up the eccentricity/disparity pairs for histogram2d
    for row, sers in data.iterrows():
        nona_sers = sers.dropna()
        temp_x = np.r_[[row]*len(nona_sers)]
        xvals.append(temp_x)
        disparities.append(nona_sers.values)

    xvals = np.concatenate(xvals)
    disparities = np.concatenate(disparities)

    H, xedges, yedges = np.histogram2d(xvals, disparities, bins=bins)

    H = H / np.float(np.sum(H))

    return H, xedges, yedges

def get_name_string(grouping_name, name):
    if type(name) == str:
        name = [name]
    else:
        name = list(name)

    if grouping_name == 'each_experiment':
        grouping_name = ''

    for i, n in enumerate(name):
        if grouping_name and not n:
            name[i] = 'raw'

    name.insert(0, grouping_name)
    name_str = '_'.join([n for n in name if n])

    return name_str

def printnow(s):
    sys.stdout.write(s+'\n')
    sys.stdout.flush()

def calc_hist(data, drange, num_bins):
    srs = pandas.Series(data.values.flatten())
    nona_srs = srs.dropna()


    # subdivide the mins into num_bins evenly spaced bins in disp_range, and
    # put everything else into bins before and after that
    binedges = np.linspace(*drange, num=num_bins)#.tolist()
    binedges = np.r_[-1e20, binedges, 1e20]
    # binedges.insert(0, -1e20)
    # binedges.append(1e20)
    # low_binedge = nona_srs.min()
    # high_binedge = nona_srs.max()
    # if low_binedge < min(binedges): binedges.insert(0, low_binedge)
    # if high_binedge > max(binedges): binedges.append(high_binedge)
    n, binedges = np.histogram(nona_srs, bins=binedges)
    n = n/np.float(np.sum(n))

    return n, binedges

def calc_weighted_hist(data, drange, num_bins):
    weight_dict = {'sandwich': .29, 'ordering_coffee': .51, 'nature_walk_1':.20}

    groups = data.groupby(level='task_id', axis=1)
    n = {}
    for i, group in groups:
        task_id = group.columns[0][2]
        n[task_id], binedges = calc_hist(group, drange, num_bins)

    n_weighted = np.zeros(n[task_id].shape)
    for task in weight_dict:
        n_weighted += n[task] * weight_dict[task]

    return n_weighted, binedges

def plot_linear_disp_hist(group, grouping_name):
    normed_data = group - group.ix[MERIDIAN_PIXEL]
    disp_range = (-2,2)
    num_bins = 76

    if grouping_name == 'all':
        n, binedges = calc_weighted_hist(normed_data, disp_range, num_bins)
    else:
        n, binedges = calc_hist(normed_data, disp_range, num_bins)

    fig, ax = plt.subplots(1,1)
    y = interp1d(binedges[:-1] + np.mean(np.diff(binedges[1:-1]))/2, n, kind='cubic')
    x = np.linspace(*disp_range, num=1000)
    lines, = ax.plot(x, y(x), lw=2)
    # rects = ax.bar(binedges[:-1], n, width=np.diff(binedges))

    ax.set_ylabel('Probability')
    ax.set_xlabel('Disparity (degrees)')
    ax.set_xlim(disp_range)
    ax.set_ylim(0,.45)

    return fig

#-------------------EDIT THESE SETTINGS-----------------------------------------

# SELECT ONE
# hist_type = '2dhists'
hist_type = 'linear_hists'

groupings = {'all':['fix_fixation_id']}
             #'each_experiment':['subject', 'fix_fixation_id', 'task_id'],
             #'subjects':['subject', 'fix_fixation_id'],
             #'tasks':['fix_fixation_id', 'task_id']}


#-------------------END EDIT SECTION--------------------------------------------

plots_dpath = '../plots_disparities'
if not os.path.exists(plots_dpath):
    os.makedirs(plots_dpath)

dispstore = pandas.io.pytables.HDFStore('disparities.h5')



if hist_type == '2dhists':
    # plot_func = plot_disp_hexbin
    plot_func = plot_disp_image
elif hist_type == 'linear_hists':
    plot_func = plot_linear_disp_hist


print "Plotting {0}...\n".format(hist_type)
for meridian in ['HMeridian', 'VMeridian']:
    for disparity_type in ['horizontal']:#, 'vertical']:
        print "Loading", meridian, disparity_type, "disparities into memory..."
        dispdata = dispstore[meridian][disparity_type].drop('Unk', axis=1, level='subject')
        dispdata = dispdata.sortlevel(0, axis=1)
        dispdata = dispdata.select(lambda x: x[2] == 'ordering_coffee'
                                                 or x[2] == 'sandwich'
                                                 or x[2] == 'nature_walk_1', axis=1)

        # for each of the ways we want to group things...
        for grouping_name, grouping in groupings.iteritems():
            subgroup_dpath = os.path.join(plots_dpath, grouping_name)
            print "Grouping by", grouping_name, ":", grouping
            groups = dispdata.groupby(level=grouping, axis=1)
            for name, group in groups:
                subgroup_dpath = os.path.join(plots_dpath, hist_type, grouping_name,
                                              meridian, disparity_type + '_disparity')
                if not os.path.exists(subgroup_dpath):
                    os.makedirs(subgroup_dpath)
                name_str = get_name_string(grouping_name, name)
                printnow("Plotting {0} for {1}".format(hist_type, name_str))

                fig = plot_func(name_str,meridian,group, grouping_name)
                plot_fname = '_'.join([meridian, disparity_type, "disparity", name_str]) + ".eps"
                fig.axes[0].set_title(plot_fname[:-4])
                fig.savefig(os.path.join(subgroup_dpath, plot_fname))

        del dispdata


# Make the 'linear' disparity histogram as in Liu et al
# print "Plotting linear hists...\n"
# for meridian in ['HMeridian', 'VMeridian']:

#     for disparity_type in ['horizontal', 'vertical']:
#         print "Loading", meridian, disparity_type, "disparities into memory..."
#         dispdata = dispstore[meridian][disparity_type].drop('Unk', axis=1, level='subject')
#         dispdata = dispdata.sortlevel(0, axis=1)
#         for grouping_name, grouping in groupings.iteritems():
#             print "Grouping by", grouping_name, ":", grouping
#             groups = dispdata.groupby(level=grouping, axis=1)
#             for name, group in groups:
#                 linear_hist_dpath = os.path.join(plots_dpath, 'linear_hists', grouping_name,
#                                                  meridian, disparity_type+"_disparity")
#                 if not os.path.exists(linear_hist_dpath):
#                     os.makedirs(linear_hist_dpath)

#                 name_str = get_name_string(grouping_name, name)
#                 printnow("Plotting linear histogram for {0}".format(name_str))

#                 fig = plot_linear_disp_hist(group)
#                 plot_fname = '_'.join([meridian, disparity_type, "disparity", name_str]) + '.eps'
#                 fig.axes[0].set_title(plot_fname[:-4])
#                 fig.savefig(os.path.join(linear_hist_dpath, plot_fname))

#         del dispdata

dispstore.close()


