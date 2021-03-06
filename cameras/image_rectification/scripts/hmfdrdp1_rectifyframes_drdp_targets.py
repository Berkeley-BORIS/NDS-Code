from rectify_frames_parallel import *


subj = 'hmfdrdp1'
calibration_date = '2012-07-18'

cam_dir = '../../stereo_calibration/data/' + subj + '/calibration_frames_' + calibration_date + '/'

inputs = ['radials_50_1','radials_50_2','radials_50_3','radials_50_4','radials_100_1','radials_100_2','radials_100_3','radials_100_4','radials_200_3','radials_450_1','radials_450_2','radials_450_3','radials_450_4']


ppservers = ()
ncpus = 6
job_server = pp.Server(ncpus,ppservers=ppservers)
print "Starting pp with", job_server.get_ncpus(), "workers"

jobs = [(input,job_server.submit(rectify_frames, (subj,cam_dir,input, ), (), ("cv","os","string" ))) for input in inputs]

for input,job in jobs:
	print "ran",job()
	
job_server.print_stats()