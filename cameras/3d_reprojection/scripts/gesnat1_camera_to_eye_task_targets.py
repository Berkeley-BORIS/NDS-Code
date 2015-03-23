from camera_image_to_eye_image_parallel import *

subj = 'gesnat1'
fixed_fixations = ''

if fixed_fixations:  # if fixed_fixations isn't empty, we might need to add an underscore
	if not fixed_fixations.startswith('_'):  # so if it doesn't start with an underscore...
		fixed_fixations = '_' + fixed_fixations  # add one to the beginning

fixation_src_file = '../../../eyes/data_processing/data/' + subj  + '/' + subj + '_fixation_points' + fixed_fixations + '.mat'
transform_src_file = '../../camera_registration/data/' + subj + '/' + subj + '_transform_1.npz'
ipdfile = '../../../eyes/ipds/' + subj[0:3] + '.txt'
ipd_cm = numpy.loadtxt(ipdfile)
frames = (0,0)
targets = 1
drift = 0

inputs = ['radials_50_1','radials_50_2','radials_100_1','radials_100_2','radials_450_1','radials_450_2']

ppservers = ()
ncpus = 6
job_server = pp.Server(ncpus,ppservers=ppservers)
print "Starting pp with", job_server.get_ncpus(), "workers"

jobs = [(input,job_server.submit(camera_image_to_eye_image, (subj,transform_src_file,fixation_src_file,ipd_cm,frames,targets,drift,input, fixed_fixations), (), ("cv","os.path","cv2","numpy","scipy.io","os","math","string","camera_to_eye_utils_parallel", ))) for input in inputs]

for input,job in jobs:
	print "ran",job()
	
job_server.print_stats()
