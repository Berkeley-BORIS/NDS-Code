from rectify_frames_parallel import *

subj = 'tkidrdp1'
inputs = ['walk_1','walk_2','walk_3']

ppservers = ()
ncpus = 3
job_server = pp.Server(ncpus,ppservers=ppservers)
print "Starting pp with", job_server.get_ncpus(), "workers"

jobs = [(input,job_server.submit(rectify_frames, (subj,input, ), (), ("cv","os","string","fnmatch" ))) for input in inputs]

for input,job in jobs:
	print "ran",job()
	
job_server.print_stats()