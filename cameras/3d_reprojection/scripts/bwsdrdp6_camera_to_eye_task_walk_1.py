from camera_image_to_eye_image_parallel import *

subj = 'bwsdrdp6'
fixation_src_file = '../../../eyes/data_processing/data/' + subj  + '/' + subj + '_fixation_points.mat'
if len(fixation_src_file.strip('.mat').split('_')) == 5:
	fixed_fixations = fixation_src_file.strip('.mat').split('_')[-1]
else:
	fixed_fixations = ''
transform_src_file = '../../camera_registration/data/' + subj + '/' + subj + '_transform_1.npz'
ipdfile = '../../../eyes/ipds/' + subj[0:3] + '.txt'
ipd_cm = numpy.loadtxt(ipdfile)
frames = (0,0)
targets = 0
drift = 0


inputs = ['task_walk_1']


ppservers = ()
ncpus = 1
job_server = pp.Server(ncpus,ppservers=ppservers)
print "Starting pp with", job_server.get_ncpus(), "workers"

jobs = [(input,job_server.submit(camera_image_to_eye_image, (subj,transform_src_file,fixation_src_file,ipd_cm,frames,targets,drift,input, fixed_fixations ), (), ("cv","os.path","cv2","numpy","scipy.io","os","math","string","camera_to_eye_utils_parallel", ))) for input in inputs]

for input,job in jobs:
	print "ran",job()
	
job_server.print_stats()
