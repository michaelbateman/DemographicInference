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


#pop_list = [ 'AB','AC', 'AJ','AK','AN', 'AR', 'AS', 'AV'  ] 
pop_list = [ 'ABM','ABS','AC', 'AJ','AK','AN', 'AR', 'AS', 'AV'  ] 
print pop_list


pop_dict = {'ABM':'Burkina Faso M','ABS':'Burkina Faso S','AC':'Uganda', 'AJ':'Guinea-Bissau','AK':'Kenya','AN':'Cameroon', 'AR':'Angola', 'AS':'Gabon', 'AV':'Guinea'}


import numpy as np
from scipy import stats
import subprocess
#call_string = 'cd individuals'
#subprocess.call(call_string, shell=True)

import sys
num_samples = int(sys.argv[1])
num_cohorts = int(sys.argv[2])


import random

for i in range(0,int(num_regions)):
    #i = int(num_regions) -1 - k
    left = i * region_size
    right = (i+1) * region_size
    window = str(left) + '-' + str(right)
    file = 'ag1000g.phase1.AR2.3L.PASS.' + window + '.vcf.gz'







    f = 'temp_file' + str(num_samples)


    k = 0

    for pop in pop_list:
    	k_pop = 0    
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
		
		
	    pop_file = 'ag1000g.phase1.AR2.3L.PASS.'+ pop + '.'+ window + '.vcf.gz'
	    call_string = 'time ../bcftools/bcftools view ' + pop_file + ' -o ' + f + '.vcf'+ ' -s ' + total_string
	    print call_string
	    subprocess.call(call_string, shell=True)
	    #subprocess.call('../bcftools/bcftools view', s1, '-h')
	    max_file = str(num_samples)+'.'+str(k_pop) +'.'+ pop + '.' + window + '.max'
	    print max_file
	    call_string = '../pbwtAllmatches/pbwt/pbwt -readVcfGT ' + f + '.vcf -write ' + f + '.pbwt -writeSites ' + f + '.sites -maxWithin > ' + max_file
	    print call_string
	    subprocess.call(call_string, shell=True)
	    call_string = 'rm ' + f + '.pbwt'
	    subprocess.call(call_string, shell=True)
	    call_string = 'rm ' + f + '.sites'
	    subprocess.call(call_string, shell=True)
	    call_string = 'rm ' + f + '.vcf'
	    subprocess.call(call_string, shell=True)

	
	    call_string = 'mv ' + max_file + ' ./maxfiles/all/'+str(num_samples)
	    subprocess.call(call_string, shell=True)
	    k_pop += 1
	    k+=1
	
    
        input_file.close()

