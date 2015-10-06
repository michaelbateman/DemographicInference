#try:
    #input_file = open("1kgPOPSsorted.list", "r")
    #print "Input file is", input_file.name
#except IOError:
    #print "Oh no, we couldn't open our file!"

#pop_list = []

#for line in input_file:
    #pop_list.append(line.strip() )

#pop_list.remove('ASW')
#pop_list.remove('MXL')



region_size = 5000000
num_regions = 50e6 / region_size

print num_regions
print int(num_regions)


pop_list = [  'ABS','ABM','AC', 'AJ','AK','AN', 'AR', 'AS', 'AV'  ] 
#pop_list = [ 'ABS','ABM'  ] 
print pop_list

import numpy as np
from scipy import stats
import subprocess
#call_string = 'cd individuals'
#subprocess.call(call_string, shell=True)


import random
import sys


for i in range(0,1):   #1,int(num_regions)):
    #i = int(num_regions) -1 - k
    left = i * region_size
    right = (i+1) * region_size
    window = str(left) + '-' + str(right)
    file = 'ag1000g.phase1.AR2.3L.PASS.' + window + '.phased.vcf.gz'
    file2 = 'ag1000g.phase1.AR2.3L.PASS.' + window + '.phased.vcf'
    
    #call_string = 'time gunzip ' + file
    #print call_string
    #subprocess.call(call_string, shell=True)

    #call_string = 'time ../tabix/bgzip ' + file2
    #print call_string
    #subprocess.call(call_string, shell=True)

    #call_string = 'time ../bcftools/bcftools index -t ' + file
    #print call_string
    #subprocess.call(call_string, shell=True)

    #call_string = 'time ../tabix/bgzip -d ' + file
    #print call_string
    #subprocess.call(call_string, shell=True)

    #call_string = 'time gzip ' + file2
    #print call_string
    #subprocess.call(call_string, shell=True)

    




    
	
    for pop in pop_list:
	
	s = pop + '.samples'
	print s
	try:
		input_file = open(s, "r")
		print "input file is", input_file.name
	except IOError:
		print "Oh no, we couldn't open our file!"
	
	print "hello"
	
	
	pop_file = 'ag1000g.phase1.AR2.3L.PASS.'+ pop + '.'+ window + '.vcf.gz'
	pop_file2 = 'ag1000g.phase1.AR2.3L.PASS.'+ pop + '.'+ window + '.vcf'
	call_string = 'gunzip ' + pop_file
        print call_string
        subprocess.call(call_string, shell=True)
	
	
	#call_string = 'time ../bcftools/bcftools view ' + file + ' -o ' + 'temp.'+pop_file + ' -S ' + s
	#print call_string
	#subprocess.call(call_string, shell=True)
	
	#call_string = 'time ../bcftools/bcftools view ' + 'temp.'+pop_file + ' -o ' + pop_file + ' -c 1'
	#print call_string
	#subprocess.call(call_string, shell=True)
	
	#call_string = 'rm ' + 'temp.'+pop_file
	#print call_string
	#subprocess.call(call_string, shell=True)
	
		
		
		
	
	input_file.close()
#call_string = 'cd ..'
#subprocess.call(call_string, shell=True)
