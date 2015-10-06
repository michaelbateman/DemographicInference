
import sys
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import time
import pylab as P

from scipy import stats	
def make_initial_conditions(s, num_haps, num_times): # writes two columns, first is times, second is initial population size
	
	
	time = [0] * num_times
	init_pop = [0] * num_times
	pairs = num_haps * (num_haps -1) / 2.0
	gap = (1.0  / num_times) * .20 *  (40000 / pairs)
	if (num_times >1): time[1] = 5.0 * gap
	if (num_times >2): time[2] = time[1] + 3.0 * gap
	for j in range(3, num_times - 2):
	    time[j] = time[j-1] + gap

	if (num_times >5):time[num_times - 2] = time[num_times - 3] + 5.0 * gap
	if (num_times >5): time[num_times - 1] = time[num_times - 2] + 10.0 * gap

	for j in range(0, num_times):
	    init_pop[j] = 50000
	    
	with open(s, 'w') as out_file:
	    for j in range(0, num_times):
		out_file.write( str(time[j]) + '\t' + str(init_pop[j]) + '\n' )
		
		
	
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
      
      
def fire( s, ic_file, num_haps, num_times, t): # runs fire on the file with name s and returns a fire file called t
	
	call_string = 'time ../fire ' + str(num_haps) + ' ' + str(num_times) + ' ' + s + ' ' + ic_file + ' > ' + t
	print call_string
	subprocess.call(call_string, shell=True)	
	
	print ' Fire is done.'
	
	return
      
      
      

def bin_data_2(s, bins): # takes data in file root.match, puts it in bins with boundaries described by bins
			# then writes to file root.binned
			#
			# bins need to be in increasing order, largest value should be roughly infinity, b[0] = 0
			# we sort the matches first so that putting into bins is (hopefully) faster
			
			
	L = []
	filename = s
	with open(filename, 'r') as input_file:
		for line in input_file:
			temp = line.strip()
			if int(temp) >=0: L.append(int(temp))
		
	L = np.sort(L)
	bin_counter = [0] * len(bins)		
	i = 0  # i = 1 so that we start with b[1] as first upper bound.  
	j = 0  # start with first match length
	while j < len(L) and i + 1 < len(bins):
		while(L[j] < bins[i+1]):
			bin_counter[i] +=1
			if j== len(L)- 1: break
			j+=1
			#print i,j, len(bins),len(L), bin_counter[i]
		i +=1
		if j== len(L)- 1: break
		
	#i = 0
	#j = 0 
	#for j in range(0, len(A)):
		#if bins[i] <= L[j]  and L[j] < bins[i+1]:
			#bin_counter[i] +=1
	filename = s + '.binned'
	with open(filename, 'w') as out_file:
		for i in range(0, len(bins)):
			out_string = str(bins[i]) + '\t' + str(bin_counter[i]) + '\n'
			out_file.write( out_string )
		
	weighted_counter = [0] * len(bins)			
	weighted_counter[0] = bin_counter[0]	
	for i in range(1, len(bins)):
		weighted_counter[i] = round( bin_counter[i])  * round(bins[i])
	
	return bin_counter, weighted_counter
		
		
		
		

def make_bins():
	bins = [0] * 150
	bins[0] = 0;
	transition = 20
	for i in range(1, transition):
	  bins[i] = 100 * i;
	for i in range(transition, len(bins)):
		bins[i] = bins[i-1] * 1.05;
	
	bins[-1] = 1e9
	
	return bins




def make_quantiles(s, quant_list):
	
	with open(s, 'r') as input_file:
		vec = []
		for line in input_file:
			temp = line.strip()
			temp = int(temp)
			vec.append(temp)
			#print temp
			
		quantiles = []
		print quant_list
		print len(vec)
		for x in quant_list:
			#np.percentile(vec, x)
			#num = stats.scoreatpercentile(vec,x)
			quantiles.append(stats.scoreatpercentile(vec, x))
	return quantiles


def read_data(s): # reads data in file called s and returns it as an array
        L = []
	filename = s
	with open(filename, 'r') as input_file:
		for line in input_file:
			temp = line.strip()
			if int(temp) >=0: L.append(int(temp))
		
	L = np.sort(L)

	return L

def write_data(s, vec): # writes the elements of vec to s, one element per line
	n = len(vec)
	with open(s, 'w') as out_file:
	  for i in range(0, n):
	    out_file.write(str(vec[i]) + '\n')
	
	return True
	      
def compile_data(pop, chrom_list, num_samples):
	    data = []
	    for chrom in chrom_list:
		
		for ind in range(0,20):
			if chrom == 'chr10000000' and pop not in ['PEL', 'ESN']:
				max_file = './maxfiles/random/' + pop +  '.' + str(num_samples) +'.'+ str(ind) + '.max'
			else:
				max_file = './maxfiles/random/' + pop + '.' + chrom + '.' + str(num_samples) +'.'+ str(ind) + '.max'
			temp = read_data(max_file)
			data.extend(temp)
			print len(temp)
	    composite_max_file  = './maxfiles/random/' + pop + '.' + str(num_samples) + '.' + 'composite' + '.max'
	    print len(data)
	    write_data(composite_max_file, data)
	    print 'data written for', pop


def bin_and_plot_data(pop, num_samples):
	max_file = './maxfiles/random/' + pop + '.' + str(num_samples) + '.' + str(ind) + '.max'
	composite_max_file  = './maxfiles/random/' + pop + '.' + str(num_samples) + '.' + 'composite' + '.max'
	#make_quantiles(match_file, quant_list)

	data = read_data(composite_max_file)


	my_bins = make_bins()
	
	start =   time.clock() 
	counter_2, weighted_counter = bin_data_2(composite_max_file, my_bins)
	print 'my sum is', sum(counter_2)
	end =  time.clock()

	print  'The SORTED binning process took', end - start, 'seconds.'


	#(n, bins, patches) = P.hist(data, bins = 100, normed = True, cumulative = False, alpha = .5, label = pop)


	#center = (bins[:-1] + bins[1:]) / 2
	plt.plot(my_bins, counter_2, label = pop)
	#plt.plot(center,n)
	P.xlim([0,1e6])
	#for i in range(0, len(my_bins)):
	#  print my_bins[i], counter_2[i]
	print pop, len(data)
	print min(data)
	print stats.scoreatpercentile(data, .5)
	print max(data)
