from disparity_sgbm_parallel import *

##EDIT
subj = 'gescafe2'

inputs = []
##MANUALLY SELECT RELEVANT FRAMES FROM EACH TASK
inputs.append(["task_ordering_coffee",(0,0)])


min_distance_mm=100
SAD_window_size=(3,5,17)
disp_12_maxdiff=2
uniqueness_ratio=10
speckle_filter = False
smoothing = True
speckle_window_size = 100
speckle_range = 32
prefilter_cap = 0
full_DP=True
denoise = True
time_smooth = False
three_d = True
save_images=True

ppservers = ()
ncpus = 1
job_server = pp.Server(ncpus,ppservers=ppservers)
print "Starting pp with", job_server.get_ncpus(), "workers"

jobs = [(input,job_server.submit(disparity_sgbm, (min_distance_mm,SAD_window_size,disp_12_maxdiff,uniqueness_ratio,speckle_filter,smoothing,speckle_window_size,speckle_range,prefilter_cap,full_DP,save_images,denoise,time_smooth,three_d,subj,input, ), (), ("cv","cv2","bottleneck","numpy","os","math","string","disparity_sgbm_parallel_utils","fnmatch" ))) for input in inputs]

for input,job in jobs:
	print "ran",job()
	
job_server.print_stats()



