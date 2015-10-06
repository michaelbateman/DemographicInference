try:
    input_file = open("1kgPOPSsorted.list", "r")
    print "Input file is", input_file.name
except IOError:
    print "Oh no, we couldn't open our file!"

pop_list = []

for line in input_file:
    pop_list.append(line.strip() )




#pop_list = [ 'YRI','GWD', 'CHS']#,'MXL','TSI', 'IBS', 'GBR', 'FIN'  ] #

print pop_list

import numpy as np
from scipy import stats


deciles_dict = {}
pop_dict = {}

len_dict = {}
start_dict = {}
midpoint_list = []



import subprocess
#call_string = 'cd individuals'
#subprocess.call(call_string, shell=True)

import sys
num_samples = int(sys.argv[1])
num_cohorts = int(sys.argv[2])

f = 'temp_file' + str(num_samples)



import random

k = 0

for pop in pop_list:
    
    s = pop + '.samples'
    print s
    try:
        input_file = open(s, "r")
        print "input file is", input_file.name
    except IOError:
    	print "Oh no, we couldn't open our file!"
    
    print "hello"
    
   
    for i in range(0, num_cohorts):
	line = random.choice(open(s).readlines())
        sam1 = line.strip()
        total_string = sam1
    
        j=1
        while (j < num_samples):	
	    line = random.choice(open(s).readlines())
	    u = line.strip()
	    if (u in total_string):
		pass
	    else:
	        total_string = total_string + ',' + u
        	j+=1
		
		
	pop_string = pop + '.hets.chr1.vcf.gz'
	call_string = 'time ../bcftools/bcftools view ' + pop_string + ' -o ' + f + '.vcf'+ ' -s ' + total_string
	print call_string
	subprocess.call(call_string, shell=True)
	#subprocess.call('../bcftools/bcftools view', s1, '-h')
	max_file = str(num_samples)+'.'+str(k) +'.'+ pop + '.max'
	print max_file
	call_string = '../manyhaps/pbwt/pbwt -readVcfGT ' + f + '.vcf -write ' + f + '.pbwt -writeSites ' + f + '.sites -maxWithin > ' + max_file
	print call_string
	subprocess.call(call_string, shell=True)
	call_string = 'rm ' + f + '.pbwt'
	subprocess.call(call_string, shell=True)
	call_string = 'rm ' + f + '.sites'
	subprocess.call(call_string, shell=True)
	call_string = 'rm ' + f + '.vcf'
	subprocess.call(call_string, shell=True)

	
	call_string = 'mv ' + max_file + ' ./maxfiles/random/'+str(num_samples)
	subprocess.call(call_string, shell=True)
	
	k+=1
	
   
    input_file.close()
#call_string = 'cd ..'
#subprocess.call(call_string, shell=True)
