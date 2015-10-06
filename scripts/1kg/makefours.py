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
	try:
            temp2 = input_file.next()
	    temp3 = input_file.next()
	    temp4 = input_file.next()
	    temp5 = input_file.next()
	    temp6 = input_file.next()
	    temp7 = input_file.next()
	    temp8 = input_file.next()
	    sam2 = temp2.strip()
	    sam3 = temp3.strip()
	    sam4 = temp4.strip()
	    sam5 = temp2.strip()
	    sam6 = temp3.strip()
	    sam7 = temp4.strip()
	    sam8 = temp4.strip()
	except StopIteration:
	    pass
        print k
        s1 = sam1+'.'+pop+'.hets.chr1.vcf.gz'
        s2 = sam2+'.'+pop+'.hets.chr1.vcf.gz'
	s3 = sam3+'.'+pop+'.hets.chr1.vcf.gz'
	s4 = sam4+'.'+pop+'.hets.chr1.vcf.gz'
	s5 = sam5+'.'+pop+'.hets.chr1.vcf.gz'
        s6 = sam6+'.'+pop+'.hets.chr1.vcf.gz'
	s7 = sam7+'.'+pop+'.hets.chr1.vcf.gz'
	s8 = sam8+'.'+pop+'.hets.chr1.vcf.gz'
	print s1
	print s2
	print s3
	print s4
	print s5
	print s6
	print s7
	print s8
	pop_string = pop + '.hets.chr1.vcf.gz'
	#call_string = 'time ../../bcftools/bcftools merge ' + s1 + ' ' + s2 +' '+ s3 +' -o combo.vcf'
	call_string = 'time ../bcftools/bcftools view ' + pop_string + ' -o combo8.vcf'+ ' -s ' + sam1+','+sam2+','+sam3+','+sam4+','+sam5+','+sam6+','+sam7+','sam8
	print call_string
	subprocess.call(call_string, shell=True)
	#subprocess.call('../bcftools/bcftools view', s1, '-h')
	max_file = 'pair.'+str(k) +'.'+ pop + '.max'
	print max_file
	call_string = '../manyhaps/pbwt/pbwt -readVcfGT combo8.vcf -write combo8.pbwt -writeSites combo8.sites -maxWithin > ' + max_file
	print call_string
	subprocess.call(call_string, shell=True)
	call_string = 'rm combo8.pbwt'
	subprocess.call(call_string, shell=True)
	call_string = 'rm combo8.sites'
	subprocess.call(call_string, shell=True)
	call_string = 'rm combo8.vcf'
	subprocess.call(call_string, shell=True)

	
	call_string = 'mv ' + max_file + ' ./maxfiles/8'
	subprocess.call(call_string, shell=True)
	
	k+=1
	
   
    input_file.close()
#call_string = 'cd ..'
#subprocess.call(call_string, shell=True)
