pop_list = [  'AV' ,'ABS','AC','AK','ABM', 'AJ','AN', 'AR', 'AS' ] 
#pop_list = [ 'AB','AC', 'AJ','AK','AN', 'AR', 'AS', 'AV'  ] 


region_size = 5000000
num_regions = 50e6 / region_size

#print num_regions
#print int(num_regions)


#pop_list = [ 'AB','AC', 'AJ','AK','AN', 'AR', 'AS', 'AV'  ] 
#pop_list = [ 'ABM','ABS','AC', 'AJ','AK','AN', 'AR', 'AS', 'AV'  ] 
#print pop_list


#pop_dict = {'ABM':'Burkina Faso M','ABS':'Burkina Faso S','AC':'Uganda', 'AJ':'Guinea-Bissau','AK':'Kenya','AN':'Cameroon', 'AR':'Angola', 'AS':'Gabon', 'AV':'Guinea'}

#pop_list = [ 'YRI'] #,'GWD', 'CHS']#,'MXL','TSI', 'IBS', 'GBR', 'FIN'  ] #

import numpy as np
from scipy import stats
import subprocess
#call_string = 'cd individuals'
#subprocess.call(call_string, shell=True)

import sys
num_samples = int(sys.argv[1])
num_cohorts = int(sys.argv[2])



chrom = str(3) + 'L'
import random


base_name = 'ag1000g.AC.phase1.AR2.3L.PASS.vcf.gz'

#for i in range(0,int(num_regions)):
for i in range(0,1):   # Extra loop here in case we want a loop over regions of the genome
    #i = int(num_regions) -1 - k
    #left = i * region_size
    #right = (i+1) * region_size
    #window = str(left) + '-' + str(right)
    #file = 'ag1000g.phase1.AR2.3L.PASS.' + window + '.vcf.gz'







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
         
	 
	
	pop_file = 'ag1000g.' + pop + '.phase1.AR2.' + chrom + '.PASS.vcf' 


	if (num_samples == 1):
	    i = 0
	    for line in input_file:
	        sam = line.strip()
	        sample_file = 'temp.'+ str(num_samples) +'.samples'
		with open(sample_file, 'w') as outfile:
		    print sam
		    outfile.write(sam+'\n')
	    #samples_already_used = ''
	    #for i in range(0, num_cohorts):
		#sample_file = 'temp.'+ str(num_samples) +'.samples'
		#with open(sample_file, 'w') as outfile:
			#j = 0
			#while ( j == 0 ):    
			    #line = random.choice(open(s).readlines())
			    #sam = line.strip()	
			    #if (sam in samples_already_used):
				#print 'Sample already used, drawing another randomly...'
				#print ' '
				#print ' '
				#print ' '
				#print ' '
				#print ' '
			    #else:
				#print samples_already_used
				#print sam
				#samples_already_used = samples_already_used + ',' + sam
				#outfile.write(sam+'\n')
				#j+=1
				
		max_file = pop_file + '.' + str(num_samples) +'.'+ str(i) + '.max'
	    	#call_string = '../../pbwtMay/pbwt/pbwt -read ' + pop_file + '.pbwt -readSites ' + pop_file + '.sites' + ' -selectSamples ' + sample_file + ' -maxWithin > ' + max_file
	    
	    	call_string = '../../new/pbwt/pbwt -readAll ' + pop_file  + ' -selectSamples ' + sample_file + ' -sfs'
	    	print call_string
	    	subprocess.call(call_string, shell=True)	
			
	    	call_string = 'mv ' + max_file + ' ../maxfiles/all/'+str(num_samples)
	    	subprocess.call(call_string, shell=True)
	    	k_pop += 1
	    	k+=1		
		i+=1		
				
	    
	else:
	    for i in range(0, num_cohorts):
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
	    
	    
	    
		
		
		
	    
	    
	    
	    	max_file = pop_file + '.' + str(num_samples) +'.'+ str(i) + '.max'
	    	#call_string = '../../pbwtMay/pbwt/pbwt -read ' + pop_file + '.pbwt -readSites ' + pop_file + '.sites' + ' -selectSamples ' + sample_file + ' -maxWithin > ' + max_file
	    
	    	call_string = '../../new/pbwt/pbwt -readAll ' + pop_file  + ' -selectSamples ' + sample_file + ' -maxWithin > ' + max_file 
	    	print call_string
	    	subprocess.call(call_string, shell=True)	
			
	    	call_string = 'mv ' + max_file + ' ../maxfiles/all/'+str(num_samples)
	    	subprocess.call(call_string, shell=True)
	    	k_pop += 1
	    	k+=1
	
    
        input_file.close()

