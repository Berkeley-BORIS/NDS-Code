"""
This script accumulates disparity data into two DataFrames for the horizontal
and vertical meridians of the eye. The DataFrames are indexed by the Helmholtz
coordinates of each pixel with left an down in the visual field assigned as
negative.

The columns of the DataFrames are heirarchically indexed as:

(disparity_type, subject, fix_fixation_id, task_id, frame)

where disparity_type is either 'horizontal' or 'vertical' disparity, subject
is the 3 initials of the subject (e.g. 'bws'), fix_fixation_id is how the
fixation points were perturbed, task_id is the name of the task (e.g.
'ordering coffee'), and frame is the frame number.
"""
import os.path
from glob import glob
import re

import numpy as np
import pandas
import matplotlib.pyplot as plt

MERIDIAN = 103  # line of pixels to pull out
FOCAL_LENGTH = 583.  # focal length of "eye" in pixels

fpath_regexp = re.compile(r"\.\./data/[a-z]{3}[a-z]{3,4}?\d_?" + \
                          r"(?P<fix_fixation_id>(corrected|random|shuffle)?)" + \
                          r"/(?P<exp_id>[a-z]{3}[a-z]{3,4}?\d)" + \
                          r"_task_(?P<task_id>\w*?)" + \
                          r"(_vertical)?_disparity\.npy")

def extract_experiment_identity(fname, regexp):

    fname_match = regexp.match(fname)
    if fname_match is None:
      print "Match failed with string:\n", fname
      return ('Unknown', 'Unknown', 'Unknown')
    exp_id = fname_match.group('exp_id')
    fix_fixation_id = fname_match.group('fix_fixation_id')
    task_id = fname_match.group('task_id')

    return exp_id, fix_fixation_id, task_id

def make_dataframe(srs_dict):
    dataframe = pandas.concat(srs_dict, axis=1)
    col_index = pandas.MultiIndex.from_tuples(dataframe.columns.tolist(),
                                              names=["disparity_type",
                                                     "subject",
                                                     "fix_fixation_id",
                                                     "task_id"])
    return dataframe.reindex(columns=col_index)

disparity_fpaths = glob("../data/*/*disparity.npy")

hmeridian = []
vmeridian = []

for disparity_fpath in disparity_fpaths:
    print "Loading", disparity_fpath
    if disparity_fpath.endswith('vertical_disparity.npy'):
      disparity_type = 'vertical'
    else:
      disparity_type = 'horizontal'

    full_dispdata = np.load(disparity_fpath)
    hm_dispdata = full_dispdata[MERIDIAN,:,:].copy()
    vm_dispdata = full_dispdata[:,MERIDIAN,:].copy()
    del full_dispdata
    (exp_id, fix_fixation_id, task_id) = extract_experiment_identity(disparity_fpath, fpath_regexp)
    subject = exp_id[:3]

    cols = pandas.MultiIndex.from_tuples([(disparity_type, subject, fix_fixation_id, task_id, i)
                                          for i in xrange(hm_dispdata.shape[1])],
                                         names=['disparity_type', 'subject',
                                                'fix_fixation_id', 'task_id', 'frame'])

    # TODO Egregious violoation of DRY follows...clean this up
    vert_merid_length = vm_dispdata.shape[0]
    # recenter pixel indices by making an array from 0 to len(meridian),
    # subtract the meridian pixel, and invert so positive is at the beginning,
    # meaning up in the visual field is positive
    shifted_pixels = -1*(np.arange(vert_merid_length) - MERIDIAN)
    hh_coords = np.rad2deg(np.arctan(shifted_pixels/FOCAL_LENGTH))
    vmeridian_index = pandas.Index(hh_coords)

    horz_merid_length = hm_dispdata.shape[0]
    # recenter pixel indices by making an array from 0 to len(meridian),
    # subtract the meridian pixel. NO INVERSION so left in the visual field is
    # negative
    shifted_pixels = np.arange(horz_merid_length) - MERIDIAN
    hh_coords = np.rad2deg(np.arctan(shifted_pixels/FOCAL_LENGTH))
    hmeridian_index = pandas.Index(hh_coords)

    hmeridian_dispdata = pandas.DataFrame(hm_dispdata, columns=cols, index=hmeridian_index)
    vmeridian_dispdata = pandas.DataFrame(vm_dispdata, columns=cols, index=vmeridian_index)

    hmeridian.append(hmeridian_dispdata)
    vmeridian.append(vmeridian_dispdata)



print "Concatenating into dataframes..."
hmeridian = pandas.concat(hmeridian, axis=1)
vmeridian = pandas.concat(vmeridian, axis=1)

dispstore = pandas.io.pytables.HDFStore('disparities.h5')
dispstore['HMeridian'] = hmeridian
dispstore['VMeridian'] = vmeridian
dispstore.close()
