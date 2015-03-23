#!/usr/bin/env ipython-2.7
import itertools
import numpy as np
import numpy.random as nprand
import os
import scipy.stats

#data type
data_type = 'corrected'
vis_field = 'upper'
task = 'nature'

if task == 'nature':
    # where to look for depth files
    dpath1 = '../data/bwsnat5_' + data_type + '/bwsnat5_task_nature_walk_1_depth.npy'
    dpath2 = '../data/gesnat1_' + data_type + '/gesnat1_task_nature_walk_1_depth.npy'
    dpath3 = '../data/tkinat1_' + data_type + '/tkinat1_task_nature_walk_1_depth.npy'

    # where to look for fixation files
    fpath1 = '../data/bwsnat5_' + data_type + '/bwsnat5_' + data_type + '_task_nature_walk_1_foveal_depth_vec.npy'
    fpath2 = '../data/gesnat1_' + data_type + '/gesnat1_' + data_type + '_task_nature_walk_1_foveal_depth_vec.npy'
    fpath3 = '../data/tkinat1_' + data_type + '/tkinat1_' + data_type + '_task_nature_walk_1_foveal_depth_vec.npy'
elif task == 'sandwich':
    # where to look for depth files
    dpath1 = '../data/bwssand2_' + data_type + '/bwssand2_task_sandwich_depth.npy'
    dpath2 = '../data/gessand1_' + data_type + '/gessand1_task_sandwich_depth.npy'
    dpath3 = '../data/tkisand2_' + data_type + '/tkisand2_task_sandwich_depth.npy'

    # where to look for fixation files
    fpath1 = '../data/bwssand2_' + data_type + '/bwssand2_' + data_type + '_task_sandwich_foveal_depth_vec.npy'
    fpath2 = '../data/gessand1_' + data_type + '/gessand1_' + data_type + '_task_sandwich_foveal_depth_vec.npy'
    fpath3 = '../data/tkisand2_' + data_type + '/tkisand2_' + data_type + '_task_sandwich_foveal_depth_vec.npy'


print "loading depth"
depth_data1tmp = np.load(dpath1)
print "loaded depth"

if vis_field == 'upper':
    depth_data1 = depth_data1tmp[:(depth_data1tmp.shape[0]/2) - 21, :]
elif vis_field == 'lower':
    depth_data1 = depth_data1tmp[(depth_data1tmp.shape[0]/2) + 21:, :]
elif vis_field == 'left':
    depth_data1 = depth_data1tmp[:, :(depth_data1tmp.shape[0]/2) + 21]
elif vis_field == 'right':
    depth_data1 = depth_data1tmp[:, (depth_data1tmp.shape[0]/2) + 21:]

del depth_data1tmp

depth_data1 = depth_data1.ravel()
depth_data1 = depth_data1[np.logical_not(np.isnan(depth_data1))]

print "loading depth"
depth_data2tmp = np.load(dpath2)
print "loaded depth"

if vis_field == 'upper':
    depth_data2 = depth_data2tmp[:(depth_data2tmp.shape[0]/2) - 21, :]
elif vis_field == 'lower':
    depth_data2 = depth_data2tmp[(depth_data2tmp.shape[0]/2) + 21:, :]
elif vis_field == 'left':
    depth_data2 = depth_data2tmp[:, :(depth_data2tmp.shape[0]/2) + 21]
elif vis_field == 'right':
    depth_data2 = depth_data2tmp[:, (depth_data2tmp.shape[0]/2) + 21:]

del depth_data2tmp

depth_data2 = depth_data2.ravel()
depth_data2 = depth_data2[np.logical_not(np.isnan(depth_data2))]

depth_data = np.concatenate([depth_data1,depth_data2])

del depth_data1
del depth_data2

print "loading depth"
depth_data3tmp = np.load(dpath3)
print "loaded depth"

if vis_field == 'upper':
    depth_data3 = depth_data3tmp[:(depth_data3tmp.shape[0]/2) - 21, :]
elif vis_field == 'lower':
    depth_data3 = depth_data3tmp[(depth_data3tmp.shape[0]/2) + 21:, :]
elif vis_field == 'left':
    depth_data3 = depth_data3tmp[:, :(depth_data3tmp.shape[0]/2) + 21]
elif vis_field == 'right':
    depth_data3 = depth_data3tmp[:, (depth_data3tmp.shape[0]/2) + 21:]

del depth_data3tmp

depth_data3 = depth_data3.ravel()
depth_data3 = depth_data3[np.logical_not(np.isnan(depth_data3))]

depth_data = np.concatenate([depth_data,depth_data3])

del depth_data3

print "loading fix" 
fix_data1 = np.load(fpath1)
print "loaded fix"

fix_data1 = fix_data1.ravel()

print "loading fix" 
fix_data2 = np.load(fpath2)
print "loaded fix"

fix_data2 = fix_data2.ravel()

print "loading fix" 
fix_data3 = np.load(fpath3)
print "loaded fix"

fix_data3 = fix_data3.ravel()

fix_data = np.concatenate([fix_data1,fix_data2,fix_data3])

#bootstrap median distributions
samples = 1000;
size = 1000;
depth_medians = np.zeros(samples)
fix_medians = np.zeros(samples)
for sample in range(samples):
    depth_data_tmp = nprand.choice(depth_data,size=size,replace=True)
    depth_medians[sample] = np.median(depth_data_tmp)
    
    fix_data_tmp = nprand.choice(fix_data,size=size,replace=True)
    fix_medians[sample] = np.median(fix_data_tmp)
    #print '%.2f' % depth_medians[sample]    
    
print 'depth', data_type, vis_field, task, '%.2f' % np.median(depth_medians)
print 'fixation', data_type, vis_field, task, '%.2f' % np.median(fix_medians)