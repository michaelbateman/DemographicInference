import sys
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import time
import pylab
from scipy import stats

num_haps = int(sys.argv[1])
num_times = int(sys.argv[2])


	
def readfire(s): # reads the output of a fire file and 
		 # returns the time vector and the population vector
	
	# input file should have two columns: first column is time in generations
	# second column is population
	
	time = []
	pop = []
	
	with open(s, 'r') as input_file:
		throwaway = input_file.readline()
		while throwaway.strip() != 'START HERE':
			throwaway = input_file.readline()
		for line in input_file:
		    temp = line.strip()
		    L = temp.split()
		    if 'START' in temp:
			print 'START'
			time = []
			pop = []
		    elif 'f' in temp:
			pass
		    elif len(L) >= 2:
			print 'no'
			print temp
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


pop_list = ['GBR', 'CEU','YRI', 'FIN', 'PEL', 'ESN']
IC_list = ['FIN']

for pop in pop_list:
    for start in IC_list:
	root = pop + '.' + str(num_haps/2) + '.composite.max.binned'
	t = root + '.fire'
	print t
	#t = pop + '.fire'
	if pop == 'ESN':
	  [T, P] = readfire(t)
	else:
	  ic_file = start + '.ic.txt'
	  t = root +'.'+ ic_file + '.fire'
	  [T, P] = readfire(t)
	
	plt.plot(np.multiply(28,T) ,P, '-o', label = pop )
	plt.yscale('log')
	plt.xlabel('years')
	#plt.title('re-started from ' + start +' curve, everyone *except* ESN improves')
	plt.title('best overall')
	#fig = plt.figure()
	#fig.set_yscale('log')
plt.legend(loc = 'lower right')
title = 'ic.' + start + '.' +str(num_haps/2) 
for pop in pop_list:
      title+= '.' + pop 


pylab.savefig(title + '.png', bbox_inches='tight')

plt.show()


