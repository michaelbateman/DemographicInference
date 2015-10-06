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


num_samples = 24

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
    for line in input_file:
        
        sam1 = line.strip()
	total_string = sam1
	try:
	    for x in range(2,num_samples+1):	
                temp = input_file.next()
	        u = temp.strip()
	        total_string = total_string + ',' + u
	    #temp3 = input_file.next()
	    #temp4 = input_file.next()
	    #temp5 = input_file.next()
	    #temp6 = input_file.next()
	    #temp7 = input_file.next()
	    #temp8 = input_file.next()
	    #temp9 = input_file.next()
	    #temp10 = input_file.next()
	    #temp11 = input_file.next()
	    #temp12 = input_file.next()
	    #temp13 = input_file.next()
	    #temp14 = input_file.next()
	    #temp15 = input_file.next()
	    #temp16 = input_file.next()
	    #sam2 = temp2.strip()
	    #sam3 = temp3.strip()
	    #sam4 = temp4.strip()
	    #sam5 = temp5.strip()
	    #sam6 = temp6.strip()
	    #sam7 = temp7.strip()
	    #sam8 = temp8.strip()
	    
	    #sam16 = temp2.strip()
	    #sam9 = temp3.strip()
	    #sam10 = temp4.strip()
	    #sam11 = temp5.strip()
	    #sam12 = temp6.strip()
	    #sam13 = temp7.strip()
	    #sam14 = temp8.strip()
	except StopIteration:
	    pass
        print k
        #s1 = sam1+'.'+pop+'.hets.chr1.vcf.gz'
        #s2 = sam2+'.'+pop+'.hets.chr1.vcf.gz'
	#s3 = sam3+'.'+pop+'.hets.chr1.vcf.gz'
	#s4 = sam4+'.'+pop+'.hets.chr1.vcf.gz'
	#s5 = sam5+'.'+pop+'.hets.chr1.vcf.gz'
        #s6 = sam6+'.'+pop+'.hets.chr1.vcf.gz'
	#s7 = sam7+'.'+pop+'.hets.chr1.vcf.gz'
	#s8 = sam8+'.'+pop+'.hets.chr1.vcf.gz'
	#print s1
	#print s2
	#print s3
	#print s4
	#print s5
	#print s6
	#print s7
	#print s8
	pop_string = pop + '.hets.chr1.vcf.gz'
	#call_string = 'time ../../bcftools/bcftools merge ' + s1 + ' ' + s2 +' '+ s3 +' -o combo.vcf'
	#call_string = 'time ../bcftools/bcftools view ' + pop_string + ' -o combo8.vcf'+ ' -s ' + sam1+','+sam2+','+sam3+','+sam4+','+sam5+','+sam6+','+sam7+','+sam8
	call_string = 'time ../bcftools/bcftools view ' + pop_string + ' -o combo24.vcf'+ ' -s ' + total_string
	print call_string
	subprocess.call(call_string, shell=True)
	#subprocess.call('../bcftools/bcftools view', s1, '-h')
	max_file = str(num_samples)+'.'+str(k) +'.'+ pop + '.max'
	print max_file
	call_string = '../manyhaps/pbwt/pbwt -readVcfGT combo24.vcf -write combo24.pbwt -writeSites combo24.sites -maxWithin > ' + max_file
	print call_string
	subprocess.call(call_string, shell=True)
	call_string = 'rm combo24.pbwt'
	subprocess.call(call_string, shell=True)
	call_string = 'rm combo24.sites'
	subprocess.call(call_string, shell=True)
	call_string = 'rm combo24.vcf'
	subprocess.call(call_string, shell=True)

	
	call_string = 'mv ' + max_file + ' ./maxfiles/'+str(num_samples)
	subprocess.call(call_string, shell=True)
	
	k+=1
	
   
    input_file.close()
#call_string = 'cd ..'
#subprocess.call(call_string, shell=True)
