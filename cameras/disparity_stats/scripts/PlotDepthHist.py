import os
import sys

import pandas
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

def plot_depth_hist(group):

    vert_merid_length = 163
    grid = -1*(np.r_[0:vert_merid_length] - vert_merid_length/2)
    hh_coords = np.rad2deg(np.arctan(grid/583.))
    depth_bin_num = 100
    ecc_bin_num = len(hh_coords)
    ecc_bins = np.linspace(-8,8, ecc_bin_num)
    depth_bins = np.linspace(0, 1000, depth_bin_num)

    count = np.zeros((ecc_bin_num-1, depth_bin_num-1))

    for path in group:
        printnow("Loading {path}".format(path=path))
        depth_data = np.load(path)
        vert_meridian = depth_data[:,np.int(depth_data.shape[1]/2.),:].copy()
        del depth_data

        assert(vert_meridian.shape[0] == 163)

        depths = pandas.DataFrame(vert_meridian, index=hh_coords)
        depths = depths.stack().dropna()
        eccentricities = [i[0] for i in depths.index.values]

        printnow("Counting...")
        H, xedges, yedges = np.histogram2d(eccentricities, depths.values, (ecc_bins, depth_bins))
        count += H

    fig, ax = plt.subplots(1,1)
    imh = ax.imshow(H,
                    extent=(min(depth_bins), max(depth_bins), min(ecc_bins), max(ecc_bins)),
                    aspect='auto', origin='lower')
    plt.colorbar(imh)
    ax.set_xlabel('Depth')
    ax.set_ylabel('Eccenticity (deg)')

    return fig


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

def plot_linear_hist(group, grouping_by, vis_field=''):


    depth_bin_num = 76
    depth_bin_range = (0,3)  # diopters
    #depth_bin_num = 38
    #depth_bin_range = (0,1.5)  # diopters
    depth_bins = np.linspace(*depth_bin_range, num=depth_bin_num)
    depth_bins = np.r_[-1, depth_bins, 2e100]

    if grouping_by == 'all':
        weight_dict = {'sandwich': .29, 'ordering_coffee': .51, 'nature_walk_1':.20}
        counts = {'sandwich':np.zeros(depth_bins.shape[0]-1),
                  'ordering_coffee':np.zeros(depth_bins.shape[0]-1),
                  'nature_walk_1':np.zeros(depth_bins.shape[0]-1)}
    else:
        weight_dict = {'sandwich': 1, 'ordering_coffee': 1, 'nature_walk_1':1}
        counts = {group.index[0][-1]:np.zeros(depth_bins.shape[0]-1)}

    for name, path in group.iteritems():
        printnow("Loading {path}".format(path=path))
        depth_data = np.load(path)
        if vis_field == 'upper':
            depth_data = depth_data[:(depth_data.shape[0]/2) - 21, 103]
        elif vis_field == 'lower':
            depth_data = depth_data[(depth_data.shape[0]/2) + 21:, 103]
        elif vis_field == 'left':
                depth_data = depth_data[103, :(depth_data.shape[0]/2) + 21]
        elif vis_field == 'right':
                depth_data = depth_data[103, (depth_data.shape[0]/2) + 21:]

        depth_data = depth_data.ravel()
        depth_data = depth_data[np.logical_not(np.isnan(depth_data))]
        depth_data = 100./depth_data  # convert cm to diopters

        c, binedges = np.histogram(depth_data, bins=depth_bins)

        counts[name[2]] += c
        del depth_data

    freqs = np.zeros(counts[name[2]].shape)
    for k, v in counts.iteritems():
        freqs += v/np.sum(v) * weight_dict[k]

    fig, ax = plt.subplots(1,1)
    y = interp1d(binedges[:-1] + np.mean(np.diff(binedges[1:-1]))/2, freqs, kind='cubic')
    x = np.linspace(*depth_bin_range, num=1000)
    lines, = ax.plot(x, y(x), lw=2)
    # rects = ax.bar(binedges[:-1], freqs, width=np.diff(binedges))

    ax.set_xlim(depth_bin_range)
    ax.set_ylim(0,.28)
    ax.set_xlabel('Diopters')
    ax.set_ylabel('Probability')

    return fig

def printnow(s):
    sys.stdout.write(s+'\n')
    sys.stdout.flush()

plots_dpath = '../plots_depths'
if not os.path.exists(plots_dpath):
    os.makedirs(plots_dpath)

depthstore = pandas.io.pytables.HDFStore('depths.h5')

depth_paths = depthstore['paths']
depth_paths = depth_paths.unstack()[['sandwich', 'ordering_coffee', 'nature_walk_1']]
depth_paths = depth_paths.stack()

#groupings = {'all':['fix_fixation_id']}
             #'each_experiment':['subject', 'fix_fixation_id', 'task_id'],
             #'subjects':['subject', 'fix_fixation_id'],
             #'tasks':['fix_fixation_id', 'task_id']}

groupings = {'tasks':['fix_fixation_id', 'task_id']}

vis_field = 'lower'  # NOTE Change this to 'upper' or 'lower' to grab only that part of the visual field

if not vis_field: vis_field = 'whole'

for grouping_name, grouping in groupings.iteritems():
    groups = depth_paths.groupby(level=grouping)

    printnow("Grouping {0}:{1}".format(grouping_name, grouping))

    for name, group in groups:
        depth_hist_dpath = os.path.join(plots_dpath, 'linear_hist', vis_field, grouping_name)
        if not os.path.exists(depth_hist_dpath):
            os.makedirs(depth_hist_dpath)

        name_str = get_name_string(grouping_name, name)

        printnow("Plotting hist for {0}".format(name_str))
        fig = plot_linear_hist(group, grouping_name, vis_field)
        plot_fname = '_'.join(["depth", name_str]) + '.pdf'
        fig.axes[0].set_title(plot_fname[:-4])
        fig.savefig(os.path.join(depth_hist_dpath, plot_fname))

depthstore.close()