	
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
