import sys
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import time

from scipy import stats

def open_file(s):
	try:
            input_file = open(s, "r")
            print "input file is", input_file.name
    	except IOError:
            print "Oh no, we couldn't open our file!"

	return

def pbwt_macs_Current(root): # runs pbwt on a macs file and returns a max file
	
	macs_file_name = root + '.macs'
	match_file_name = root + '.match'
	#call_string = './pbwtCurrent/pbwt/pbwt -readMacs ' + root + '.macs -writeAll ' + root + ' -maxWithin > ' + root + '.match'
	call_string = './pbwtCurrent/pbwt/pbwt -readMacs ' + root + '.macs -writeAll ' + root + ' -maxWithin > ' + root + '.match'  
   	print call_string
    	subprocess.call(call_string, shell=True)
	
	return


def pbwt_macs_May(root): # runs pbwt on a macs file and returns a max file
	
	macs_file_name = root + '.macs'
	match_file_name = root + '.match'
	#call_string = './pbwtCurrent/pbwt/pbwt -readMacs ' + root + '.macs -writeAll ' + root + ' -maxWithin > ' + root + '.match'
	call_string = './pbwtMay/pbwt/pbwt -readMacs ' + root + '.macs -writeAll ' + root + ' -maxWithin > ' + root + '.match'  
   	print call_string
    	subprocess.call(call_string, shell=True)
	
	return


def sim_macs(s, length, num_haps): 
	
	#pop_size = 2e4
	chrom_length = length
	t = .001
	r = .001
	
	macs_file_name = s + '.macs'
	
	call_string = './pbwtMay/pbwt/macs/macs ' + str(num_haps) + ' ' + str(chrom_length) + ' -t ' + str(t) + ' -r ' + str(r) + ' > ' + macs_file_name + ' 2>/dev/null'
	print call_string
	print 'Running macs...' 
	subprocess.call(call_string, shell=True)	
	print call_string
	return 

	
	
def fire( s, num_times, num_haps, t): # runs fire on the file with name s and returns a fire file called t
	
	call_string = 'fire ' + str(num_haps) + ' ' + str(num_times) + ' ' + s+'.match' #+ ' > ' + t
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
		throwaway = input_file.readline()
		while throwaway.strip() != 'START HERE':
			throwaway = input_file.readline()
		for line in input_file:
			print 'hello'
			temp = line.strip()
			a,b = temp.split()
			time.append(float(a))
			pop.append(float(b))
			print a, b	
	
	print 'readfire is done'
	
	return [time, pop]
		
		

def convert_data(vec,s): 
	
	(n, bins) = np.histogram(vec, bins = 100, normed = True)
        Y = np.zeros((100,))
        X = np.zeros((100,))
        delta = np.zeros((100,))
        mass = np.zeros((100,))

    
        total = 0
    
        for i in range(0, 100):
       	    Y[i] = np.exp(bins[i]) * n[i]
	    X[i] = bins[i]
	    delta[i] = bins[i+1] - bins[i]
            mass[i] = Y[i] * delta[i]
	    total +=mass[i]
	
	with open(s, 'w') as out_file:
		for i in range(0, len(Y)):
			out_file.write( str(Y[i]) + '\n')

	return X,Y
        #print 'Total = ', total
        #plt.plot(X, np.divide(Y, total*1.0), label = pop_dict[pop])
	


def bin_data(root, bins): # takes data in file root.match, puts it in bins with boundaries described by bins
			# then writes to file root.binned
			#
			# bins need to be in increasing order, largest value should be roughly infinity
			# we sort the matches first so that putting into bins is (hopefully) faster
			
			
	L = []
	filename = root + '.match'
	with open(filename, 'r') as input_file:
		for line in input_file:
			temp = line.strip()
			L.append(int(temp))
			
	#L = np.sort(L)
	bin_counter = [0] * len(bins)
	for j in range(0, len(L)):
		
		for i in range(0, len(bins) - 1):
			if bins[i] <= L[j]  and L[j] < bins[i+1]:
				bin_counter[i] +=1
				break


	filename = root + '.binned'
	with open(filename, 'w') as out_file:
		for i in range(0, len(bins)):
			if i == 0:
				out_string = str(round(bins[i])) + '\t' + str(round(bin_counter[i])) + '\n'
			else:
				#out_string = str(round(bins[i])) +  '\t' + str(round(bin_counter[i]) * round(bins[i])) + '\n'
				out_string = str(round(bins[i])) +  '\t' + str(round(bin_counter[i])) + '\n' #* round(bins[i])) 
			out_file.write( out_string )
		

	return bin_counter
		
		
		


def bin_data_2(root, bins): # takes data in file root.match, puts it in bins with boundaries described by bins
			# then writes to file root.binned
			#
			# bins need to be in increasing order, largest value should be roughly infinity, b[0] = 0
			# we sort the matches first so that putting into bins is (hopefully) faster
			
			
	L = []
	filename = root + '.match'
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
	filename = root + '.binned'
	with open(filename, 'w') as out_file:
		for i in range(0, len(bins)):
			if i == 0:
				out_string = str(round(bins[i])) + '\t' + str(round(bin_counter[i])) + '\n'
			else:
				out_string = str(round(bins[i])) +  '\t' + str( round(bin_counter[i]) * round(bins[i]) ) + '\n'
				#out_string = str(round(bins[i])) +  '\t' + str(round(bin_counter[i])) + '\n' #* round(bins[i])) 
			out_file.write( out_string )
		
	weighted_counter = [0] * len(bins)			
	weighted_counter[0] = bin_counter[0]	
	for i in range(1, len(bins)):
		weighted_counter[i] = round( bin_counter[i])  * round(bins[i])
	
	return bin_counter, weighted_counter
		
		
		
		

def make_bins():
	bins = [0] * 1000
	for i in range(0, 100):
		bins[i] = i;
	for i in range(100, len(bins) - 1):
		bins[i] = bins[i-1] * 1.01
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




num_haps = int(sys.argv[1])

chrom_length = 1e9

file_root = 'temp'

macs_file = file_root + '.macs'
match_file = file_root + '.match'
num_times = 3
t = file_root + '.fire'

sim_macs(file_root, chrom_length, num_haps)




pbwt_macs_May(file_root)

bins = make_bins()
start =   time.clock() 
counter_2, weighted_counter = bin_data_2(file_root, bins)
end =  time.clock()

print  'The SORTED binning process took', end - start, 'seconds.'


#fire(file_root, num_times, t)
#time, pop = readfire(t)

#total = np.sum(counter_2)
##print 'total number of '
#plt.plot( np.log(bins),np.divide(counter_2, float(total) ), 'r' )
##plt.plot(time,pop)

#pbwt_macs_Current(file_root)

#bins = make_bins()
#start =   time.clock() 
#counter_2, weighted_counter = bin_data_2(file_root, bins)
#end =  time.clock()

#print  'The SORTED binning process took', end - start, 'seconds.'


#quant_list = [25, 50, 75]
quant_list = [ 50 ]
#quant_list = []
#for i in range(1, 1000):
	#quant_list.append(90 + .01*i)


quant = make_quantiles(match_file, quant_list)

quant_file = file_root+'.quant'
with open(quant_file, 'w') as output_file:
	for i in range(0, len(quant)):
		output_file.write(str(quant_list[i]) + '\t' + str(quant[i])+ '\n')
		

print 'The quantiles are', quant

#fire(file_root, num_times, t)
#time, pop = readfire(t)

total = np.sum(counter_2)
weighted_total = np.sum(weighted_counter)
plt.plot(np.log(bins),np.divide(weighted_counter, float(weighted_total) ), 'b' )
#plt.plot(time,pop)


#print X
#print Y

#A = []
#B = []
#for i in range(0, 100):
	#A.append(i)
	#B.append(i*i)

#plt.plot(time,pop)
#plt.show()
