try:
    input_file = open("1kgPOPSsorted.list", "r")
    print "Input file is", input_file.name
except IOError:
    print "Oh no, we couldn't open our file!"

pop_list = []

for line in input_file:
    pop_list.append(line.strip() )




pop_list = [ 'CEU', 'FIN', 'YRI']#,'MXL','TSI', 'IBS', 'GBR', 'FIN'  ] #

pop_list = ['GBR', 'CEU', 'FIN', 'YRI']
pop_list = ['PEL', 'ESN']
print pop_list

import numpy as np
from scipy import stats


deciles_dict = {}
pop_dict = {}

len_dict = {}
start_dict = {}
midpoint_list = []



import subprocess
import random
import sys
num_samples = int(sys.argv[1])
num_cohorts = int(sys.argv[2])

f = 'temp_file' + str(num_samples)




chrom_list =  ['chr1','chr22', 'chr21', 'chr20', 'chr19', 'chr18', 'chr17', 'chr16']#, 'chr15', 'chr14', 'chr13']


k = 0
for chrom in chrom_list:
  for pop in pop_list:
    
    s = pop + '.samples'
    print s
    try:
        input_file = open(s, "r")
        print "input file is", input_file.name
    except IOError:
    	print "Oh no, we couldn't open our file!"
    
    print "hello"
    
   
    for i in range(10, 10 + num_cohorts):
	sample_file = 'temp.'+ str(num_samples) +'.samples'
	with open(sample_file, 'w') as outfile:
		line = random.choice(open(s).readlines())
        	sam1 = line.strip()
        	total_string = sam1
    		outfile.write(sam1+'\n')
        	j=1
        	while (j < num_samples):	
	    		line = random.choice(open(s).readlines())
	    		u = line.strip()
	    		if (u in total_string):
				pass
	    		else:
			        total_string = total_string + ',' + u
				outfile.write(u+'\n')
        			j+=1
		
	
	#sample_file = 'temp.'+ str(num_samples) +'.samples'
		#with open(sample_file, 'w') as outfile:
			#line = random.choice(open(s).readlines())
			#sam1 = line.strip()
			#total_string = sam1
			#outfile.write(sam1+'\n')
			#j=1
			#while (j < num_samples):	
			    #line = random.choice(open(s).readlines())
			    #u = line.strip()
			    #if (u in total_string):
				#pass
			    #else:
				#total_string = total_string + ',' + u
				#outfile.write(u+'\n')
				#j+=1
	
	max_file = pop + '.' + chrom + '.' + str(num_samples) +'.'+ str(i) + '.max'
	#pop_string = pop + '.hets.chr1.vcf.gz'
	pop_string = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf'
	pbwt_pop_string = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.pbwt'
	#call_string = 'time ../bcftools/bcftools view ' + pop_string + ' -o ' + f + '.vcf'+ ' -s ' + total_string
	
	
	pop_string = pop + '.' + chrom
	
	call_string = 'time ../pbwtDec/pbwt/pbwt -readAll ' + pop_string  + ' -selectSamples ' + sample_file + ' -maxWithin > ' + max_file 
	print call_string
	subprocess.call(call_string, shell=True)
	
	#call_string = 'time ../new/pbwt/pbwt -readAll ' + pop_string  + ' -selectSamples ' + sample_file + ' -maxWithin > ' + max_file 
	#print call_string
	#subprocess.call(call_string, shell=True)
	
	#subprocess.call('../bcftools/bcftools view', s1, '-h')
	#max_file = str(num_samples)+'.'+str(k) +'.'+ pop + '.max'
	#print max_file
	#call_string = 'time ../pbwtMay/pbwt/pbwt -readVcfGT ' + f + '.vcf -write ' + f + '.pbwt -writeSites ' + f + '.sites -maxWithin > ' + max_file
	#print call_string
	#subprocess.call(call_string, shell=True)
	
	#call_string = 'rm ' + f + '.pbwt'
	#subprocess.call(call_string, shell=True)
	#call_string = 'rm ' + f + '.sites'
	#subprocess.call(call_string, shell=True)
	#call_string = 'rm ' + f + '.vcf'
	#subprocess.call(call_string, shell=True)

	
	call_string = 'mv ' + max_file + ' ./maxfiles/random/'#+str(num_samples)
	subprocess.call(call_string, shell=True)
	
	k+=1
	
   
    input_file.close()
#call_string = 'cd ..'
#subprocess.call(call_string, shell=True)
