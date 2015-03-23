import os
import re
from glob import glob
import sys

import numpy as np
import pandas
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import scipy.io as spio

MERIDIAN = 103  # line of pixels to pull out
FOCAL_LENGTH = 583.  # focal length of "eye" in pixels

def calc_hist(data, drange, num_bins):
    # srs = pandas.Series(data.flatten())
    # nona_srs = srs.dropna()
    if np.isnan(data).any():
        raise ValueError('NaNs in histogram data!')


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
    n, binedges = np.histogram(data, bins=binedges)
    n = n/np.float(np.sum(n))

    return n, binedges

def calc_weighted_hist(data, drange, num_bins):
    weight_dict = {'sandwich': .29, 'ordering_coffee': .51, 'nature_walk_1':.20}

    n = {}
    for task_id, disps in data.iteritems():
        n[task_id], binedges = calc_hist(disps, drange, num_bins)

    n_weighted = np.zeros(n[task_id].shape)
    for task in weight_dict:
        n_weighted += n[task] * weight_dict[task]

    return n_weighted, binedges

def plot_linear_disp_hist(normed_data, grouping_name='all'):
    # normed_data = {}
    # for task_id, disps in data.iteritems():
    #     normed_data[task_id] = disps - disps[:,:,MERIDIAN]
    print min(normed_data), max(normed_data)
    disp_range = (-.05,.05)
    num_bins = 76

    if grouping_name == 'all':
        n, binedges = calc_weighted_hist(normed_data, disp_range, num_bins)
    else:
        n, binedges = calc_hist(normed_data, disp_range, num_bins)

    fig, ax = plt.subplots(1,1)
    # y = interp1d(binedges[:-1] + np.mean(np.diff(binedges[1:-1]))/2, n, kind='cubic')
    x = np.linspace(*disp_range, num=1000)
    # lines, = ax.plot(x, y(x), lw=2)
    rects = ax.bar(binedges[:-1], n, width=np.diff(binedges))

    ax.set_ylabel('Probability')
    ax.set_xlabel('Disparity (degrees)')
    ax.set_xlim(disp_range)
    ax.set_ylim(0,.45)

    return fig

def save_hist_data(fig, name, save_dpath='.'):

    ax = fig.axes[0]
    patches = ax.patches
    left_binedges = []
    binwidths = []
    vergs = []
    for patch in patches:
        left_binedges.append(patch.get_x())
        vergs.append(patch.get_height())
        binwidths.append(patch.get_width())

    data = np.c_[left_binedges, vergs, binwidths]
    save_fpath = os.path.join(save_dpath, ax.get_title()+'.mat')
    print 'Saving hist plot data to', save_fpath
    spio.savemat(save_fpath, {name:data})

if __name__ == '__main__':

    if 'disparities' in dir():
        print 'Data already detected, reload not necessary'
        if not os.path.exists(save_dpath):
            print 'Making directories...'
            os.makedirs(save_dpath)
        print 'Calculating histogram...'
        sys.stdout.flush()
        linear_hist_fig = plot_linear_disp_hist(disparities)
        linear_hist_fig.axes[0].set_title(save_fname[:-4])
        print 'Saving histogram to', save_fpath
        linear_hist_fig.savefig(save_fpath)
        save_hist_data(linear_hist_fig, 'disp_data', save_dpath)
        sys.exit()


    disparity_type = 'vertical_disparity'
    if disparity_type == 'horizontal_disparity':
        disparity_type_fname = 'disparity'
    else:
        disparity_type_fname = disparity_type

    fix_fixation_id = 'corrected'
    tasks_of_interest = ('ordering_coffee', 'nature_walk_1', 'sandwich')

    fpath_regexp = re.compile(r"\.\./data/[a-z]{3}[a-z]{3,4}?\d" + \
                              r"_" + fix_fixation_id + \
                              r"/(?P<exp_id>[a-z]{3}[a-z]{3,4}?\d)" + \
                              r"_task_(?P<task_id>\w*?)" + \
                              r"_" + disparity_type_fname + r"\.npy")

    save_dpath = os.path.join('..', 'plots_disparities', 'linear_hists', 'all', '7deg', disparity_type)
    if not os.path.exists(save_dpath):
        print 'Making directories...'
        os.makedirs(save_dpath)
    save_fname = '_'.join(['7deg', disparity_type, 'all', fix_fixation_id]) + '.eps'
    save_fpath = os.path.join(save_dpath, save_fname)

    disparity_fpaths = glob("../data/*/*disparity.npy")

    disparities = {}
    usable_fpaths = []
    for fpath in disparity_fpaths:
        fpath_match = fpath_regexp.match(fpath)
        if fpath_match is None:
            print "Match failed with string:"
            print fpath
            print ''
            continue

        task_id = fpath_match.group('task_id')
        if task_id not in tasks_of_interest:
            print 'Skipping task_id', task_id
            print ''
            continue

        usable_fpaths.append((fpath, task_id))

    print 'Using paths'
    print usable_fpaths
    print ''
    sys.stdout.flush()

    for fpath, task_id in usable_fpaths:
        print 'Using path:\n', fpath
        print ''
        sys.stdout.flush()
        disp_mat = np.load(fpath)
        disp_mat = np.array(disp_mat, dtype=np.float32)
        disp_mat = disp_mat - disp_mat[MERIDIAN, MERIDIAN, :]
        # create mask
        n = disp_mat.shape[0]
        r = .7*103
        y,x = np.ogrid[-MERIDIAN:n-MERIDIAN, -MERIDIAN:n-MERIDIAN]
        mask = x*x + y*y <= r*r
        mask = np.reshape(mask, (207,207,1))
        mask = np.repeat(mask, disp_mat.shape[2], axis=2)
        disp_mat = disp_mat[mask]
        disp_mat = pandas.Series(disp_mat.flatten())
        disp_mat = disp_mat.dropna().values
        current_disparities = disparities.get(task_id)
        if current_disparities is not None:
            disparities[task_id] = np.concatenate((current_disparities, disp_mat), axis=2)
        else:
            disparities[task_id] = disp_mat

        del disp_mat

    print 'Calculating histogram...'
    sys.stdout.flush()
    linear_hist_fig = plot_linear_disp_hist(disparities)
    linear_hist_fig.axes[0].set_title(save_fname[:-4])
    print 'Saving histogram to', save_fpath
    linear_hist_fig.savefig(save_fpath)
    save_hist_data(linear_hist_fig, 'disp_data', save_dpath)