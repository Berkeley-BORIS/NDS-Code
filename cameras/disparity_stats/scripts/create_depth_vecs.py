import sys
import os
import re

import numpy as np

exp_dname_regexp = re.compile(r"(?P<exp_id>^[a-z]{3}[a-z]{3,4}?\d)_?(?P<fix_fixation_id>(corrected|random|shuffle)?)$")
depth_matrix_fname_regexp = re.compile(
    r"(?P<exp_id>^[a-z]{3}[a-z]{3,4}?\d)_task_(?P<task_id>\w*?)_depth\.npy")

if __name__ == '__main__':
    # setup paths to depth npy arrays
    data_dpath = os.path.join('..', 'data')

    exp_list = [exp_dname for exp_dname in os.listdir(data_dpath)
                if exp_dname_regexp.match(exp_dname)]
    for experiment in exp_list:
        exp_dpath = os.path.join(data_dpath, experiment)
        exp_name_match = exp_dname_regexp.match(experiment)
        exp_id = exp_name_match.group('exp_id')
        fix_fixation_id = exp_name_match.group('fix_fixation_id')
        try:
            data_mat_list = os.listdir(exp_dpath)
            # filter mat list for depth npy arrays
            data_mat_names = [mat_name for mat_name in data_mat_list if mat_name.endswith('depth.npy')]
        except OSError as e:
            print 'OSError'
            print e
            print 'continuing...\n'
            continue

        print 'Creating depth vec for', experiment, '...'
        for mat_name in data_mat_names:
            task_id = depth_matrix_fname_regexp.match(mat_name).group('task_id')

            depth_data_fpath = os.path.join(exp_dpath, mat_name)
            used_frames_fname = '_'.join([exp_id, 'task', task_id, 'frames_used']) + '.npy'

            depth_mat = np.load(depth_data_fpath)
            used_frames = np.load(os.path.join(exp_dpath, used_frames_fname))

            center_y = depth_mat.shape[0]/2
            center_x = depth_mat.shape[1]/2
            foveal_depth_vec = depth_mat[center_y, center_x, :]
            foveal_depth_vec = np.c_[used_frames, foveal_depth_vec]

            depth_vec_fname_parts = filter(None,
                                    [exp_id, fix_fixation_id, 'task', task_id, 'foveal_depth_vec.npy'])

            np.save(os.path.join(exp_dpath, '_'.join(depth_vec_fname_parts)), foveal_depth_vec)
        print 'Done!\n'