from rectify_frames_parallel import *

subj = 'bwssand2'
inputs = ['radials_50_1','radials_50_2','radials_100_1','radials_100_2','radials_450_1','radials_450_2']

ppservers = ()
ncpus = 3
job_server = pp.Server(ncpus,ppservers=ppservers)
print "Starting pp with", job_server.get_ncpus(), "workers"

jobs = [(input,job_server.submit(rectify_frames, (subj,input, ), (), ("cv","os","string","fnmatch" ))) for input in inputs]

for input,job in jobs:
	print "ran",job()
	
job_server.print_stats()