from rectify_frames_parallel import *

subj = 'gesdrdp1'
inputs = ['walk_2']

ppservers = ()
ncpus = 1
job_server = pp.Server(ncpus,ppservers=ppservers)
print "Starting pp with", job_server.get_ncpus(), "workers"

jobs = [(input,job_server.submit(rectify_frames, (subj,input, ), (), ("cv","os","string","fnmatch" ))) for input in inputs]

for input,job in jobs:
	print "ran",job()
	
job_server.print_stats()