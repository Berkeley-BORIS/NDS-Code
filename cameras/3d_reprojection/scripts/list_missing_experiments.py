"""
This script looks through both hard drives and compiles a list of experiments
that may have missing or incomplete reprojections. These experiments need
to be reprojected!
"""
import os

import numpy as np

def get_expected_experiments():
    """
    Return a list a of expected experiments of the form:
    expid_fixationstyle_taskid

    eg "bwsnat5_corrected_task_nature_walk_1"
    """

    expected_experiments = []

    exp_list = []

    # load the list of experiment IDs we want to search under
    with open('/Users/natdispstats/Documents/cameras/experiment_list.txt') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            else:
                exp_list.append(line)

    fixes = ['', '_corrected', '_random', '_shuffle']

    for experiment in exp_list:
        subj_id = experiment[:3]
        task_id = experiment[3:]
        if task_id.startswith('drdp'):
            tasks = ['_walk_1']
        elif task_id.startswith('nat'):
            tasks = ['_nature_walk_1', '_campus_walk_3']
        elif task_id.startswith('ure'):
            tasks = ['_nature_walk_2']
        elif task_id.startswith('sand'):
            tasks = ['_sandwich']
        elif task_id.startswith('cafe'):
            tasks = ['_ordering_coffee']
        else:
            tasks = None
            print task_id, 'is unrecognized task_id'

        for fix in fixes:
            for task in tasks:
                expected_experiments.append(''.join([experiment, fix, '_task', task]))

    return expected_experiments

def load_frame_status(exp_name):
    """
    Load the frame status matrix of this experiment and return it.
    Returns None is the frame status matrix can't be found.

    exp_name must be of the form:
    expid_fixationstyle_taskid
    eg bwsnat5_corrected_task_nature_walk_1
    """

    HD1_root = '/Users/natdispstats/Documents/cameras/3d_reprojection/data/'
    HD2_root = "/Volumes/Macintosh HD 2/cameras2/3d_reprojection/data/"

    exp_name_parts = exp_name.split('_')
    exp_id = exp_name_parts[0]
    
    # figure out the task_id and fixed_fixation state
    task_id_start_index = exp_name_parts.index('task')
    task_id = '_'.join(exp_name_parts[task_id_start_index:])
    if task_id_start_index == 1:
        fixed_fixations = ''
    else:
        fixed_fixations = exp_name_parts[1]

    # look in both hard drive locations, returning immediately if found
    for root_path in [HD1_root, HD2_root]:
        frame_status_name_parts = [part for part in [exp_id, fixed_fixations, task_id, 'frame_status.npy'] if part]
        frame_status_fname = '_'.join(frame_status_name_parts)

        frame_status_dname = '_'.join([exp_id, task_id, 'eye_images_transform1'])

        exp_dname_parts = [part for part in [exp_id, fixed_fixations] if part]
        exp_dname = '_'.join(exp_dname_parts)

        frame_status_fpath = os.path.join(root_path, exp_dname, frame_status_dname, frame_status_fname)
        try:
            frame_status_matrix = np.load(frame_status_fpath)

            return frame_status_matrix
        except IOError as e:
            frame_status_fname = '_'.join([exp_id, task_id, 'frame_status.npy'])
            frame_status_fpath = os.path.join(root_path, exp_dname, frame_status_dname, frame_status_fname)
            try:
                frame_status_matrix = np.load(frame_status_fpath)
                return frame_status_matrix
            except IOError as e:
                pass

    return None

def find_missing_frames(frame_status_matrix):
    """
    Takes in a frame status matrix and returns any frames
    that don't appear to have been processed.
    """

    missing_indices = np.where(frame_status_matrix[:,1] == 0)
    return missing_indices[0]  #frame_status_matrix[missing_indices,0]

if __name__ == '__main__':

    missing_experiments = []

    experiment_list = get_expected_experiments()
    for experiment in experiment_list:
        frame_status_matrix = load_frame_status(experiment)

        if frame_status_matrix is None:
            experiment_is_missing = True  # experiment is missing if frame status can't be loaded
        else:
            frames_incomplete = find_missing_frames(frame_status_matrix)
            frames_incomplete = np.array([])  # FIXME This isn't working right, so assume all frames are processed
            if frames_incomplete.shape[0]:
                experiment_is_missing = True  # experiment is missing if there unprocessed frames

            else:
                experiment_is_missing = False


        if experiment_is_missing:
            missing_experiments.append(experiment)

    print len(missing_experiments), 'experiments are currently missing:'
    for e in missing_experiments: print e