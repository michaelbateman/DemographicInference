import sys

from mdb import *

num_haps = int(sys.argv[1])
num_samples = num_haps / 2
#num_cohorts = int(sys.argv[2])
num_times = int(sys.argv[2])

#call_string = 'python managedata.py'
	#print call_string
	#subprocess.call(call_string, shell=True)
	
pop_list = ['FIN', 'GBR', 'CEU','YRI',  'PEL', 'ESN']
IC_list = ['flat50000']
for pop in pop_list:
    for starting_guess in IC_list:
	ic_file = starting_guess +'.txt'
	make_initial_conditions(ic_file,  num_haps, num_times)
	root = pop + '.' + str(num_haps/2) + '.composite.max.binned'
	s = './maxfiles/random/' + root
	
	
	#ic_file = 'CEU.long.ic.txt'
	t = root +'.'+str(num_times) + '.'+ ic_file + '.fire'
	call_string = 'gcc ../fire.c -o fire -lm -lgsl -lgslcblas'
	print call_string
	subprocess.call(call_string, shell=True)	

	fire(s, ic_file, num_haps, num_times, t)

	
	


