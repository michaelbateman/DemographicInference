import sys
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import time

from scipy import stats


num_haps = int(sys.argv[1])
num_times = int(sys.argv[2])
#my_pop = sys.argv[3]

pop_list = ['FIN', 'GBR', 'CEU','YRI',  'PEL', 'ESN']
#pop_list = ['PEL', 'ESN']
	
def fire( s, ic_file, num_haps, num_times, t): # runs fire on the file with name s and returns a fire file called t
	
	call_string = 'time ../fire ' + str(num_haps) + ' ' + str(num_times) + ' ' + s + ' ' + ic_file + ' > ' + t
	print call_string
	subprocess.call(call_string, shell=True)	
	
	print ' Fire is done.'
	
	return

	
def readfire(s): # reads the output of a fire file and 
		 # returns the time vector and the population vector
	
	# input file should have two columns: first column is time in generations
	# second column is population
	
	time = []
	pop = []
	
	with open(s, 'r') as input_file:
		for line in input_file:
		    temp = line.strip()
		    if 'f()' in temp:
			time = []
			pop = []
			temp = line.strip()
			a,b = temp.split()
			time.append(float(a))
			pop.append(float(b))
	
	#with open(s, 'r') as input_file:
		#throwaway = input_file.readline()
		#while throwaway.strip() != 'START HERE':
			#throwaway = input_file.readline()
		#for line in input_file:
			#print 'hello'
			#temp = line.strip()
			#a,b = temp.split()
			#time.append(float(a))
			#pop.append(float(b))
			#print a, b	
	
	print 'readfire is done'
	
	return [time, pop]

def make_initial_conditions(s, num_haps, num_times): # writes two columns, first is times, second is initial population size
	time = [0] * num_times
	init_pop = [0] * num_times
	pairs = num_haps * (num_haps -1) / 2.0
	gap = (1.0  / num_times) * 4.0 *  (20000 / pairs)
	time[1] = 5.0 * gap
	time[2] = time[1] + 3.0 * gap
	for j in range(3, num_times - 2):
	    time[j] = time[j-1] + gap

	time[num_times - 2] = time[num_times - 3] + 5.0 * gap
	time[num_times - 1] = time[num_times - 2] + 10.0 * gap

	for j in range(0, num_times):
	    init_pop[j] = 50000
	    
	with open(s, 'w') as out_file:
	    for j in range(0, num_times):
		out_file.write( str(time[j]) + '\t' + str(init_pop[j]) + '\n' )
		
	
IC_list = ['CEU']
for pop in pop_list:
    for starting_guess in IC_list:
	ic_file = 'neutral.txt'
	make_initial_conditions(ic_file,  num_haps, num_times)
	root = pop + '.' + str(num_haps/2) + '.composite.max.binned'
	s = './maxfiles/random/' + root
	
	
	#ic_file = 'CEU.long.ic.txt'
	t = root +'.'+ ic_file + '.fire'
	call_string = 'gcc ../fire.c -o fire -lm -lgsl -lgslcblas'
	print call_string
	subprocess.call(call_string, shell=True)	

	fire(s, ic_file, num_haps, num_times, t)
