import sys
import matplotlib as plt
import numpy as np


def open_file(s):
	try:
            input_file = open(s, "r")
            print "input file is", input_file.name
    	except IOError:
            print "Oh no, we couldn't open our file!"


def pbwt_max(): # runs pbwt and returns a maxfile
	

def sim_max(t): 
	
	pop_size = 1e5
	chrom_length = 1e7
	t = .001
	r = .001
	
	max_file_name = t
	
	call_string = 'macs 2 ' + chrom_length + ' -t ' + t + ' -r ' + r + ' > ' + max_file + '2>/dev/null'
	print call_string
	print 'Running macs...' 
	subprocess.call(call_string, shell=True)	
	
	return 

	
	
def fire( s, num_times, t): # runs fire on the file with name s and returns a fire file called t
	
	call_string = 'fire 2 ' + num_times + ' ' + s + ' > ' + t
	print call_string
	subprocess.call(call_string, shell=True)	
	
	return

	
def readfire(s): # reads the output of a fire file and 
		 # returns the time vector and the population vector
	
	# input file should have two columns: first column is time in generations
	# second column is population
	
	X = []
	Y = []
	
	with open(s, 'r') as input_file:
		for line in input_file:
			temp = line.strip()
			a,b = temp.split()
			X.append(float(a))
			Y.append(float(b))
			print a, b	
	
	return X, Y
		
		


s = 'temp.max'
num_times = 5
t = 'output.fire'

sim_max(s)
fire(s, num_times, t)