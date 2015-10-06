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

for pop in pop_list:
    k = 0
    s = pop + '.samples'
    try:
        input_file = open(s, "r")
        print "input file is", input_file.name
    except IOError:
    	print "Oh no, we couldn't open our file!"
    
 #   print "hello"
    for line in input_file:
        k+=1
        sam1 = line.strip()
	try:
            temp = input_file.next()
	    sam2 = temp.strip()
	except StopIteration:
	    pass
        print k
        s1 = sam1+'.'+pop+'.hets.chr1.vcf.gz'
        s2 = sam2+'.'+pop+'.hets.chr1.vcf.gz'
	print s1
	print s2
	pop_string = pop + '.hets.chr1.vcf.gz'
#	call_string = '../bcftools/bcftools merge ' + s1 + ' ' + s2 +' -o combo.vcf'
	call_string = 'time ../bcftools/bcftools view ' + pop_string + ' -o combo.vcf'+ ' -s ' + sam1+','+sam2 
	print call_string
	subprocess.call(call_string, shell=True)
	#subprocess.call('../bcftools/bcftools view', s1, '-h')
	max_file = 'pair.'+str(k) +'.'+ pop + '.max'
	print max_file
	call_string = '../manyhaps/pbwt/pbwt -readVcfGT combo.vcf -write combo.pbwt -writeSites combo.sites -maxWithin > ' + max_file
	print call_string
	subprocess.call(call_string, shell=True)
	call_string = 'rm combo.pbwt'
	subprocess.call(call_string, shell=True)
	call_string = 'rm combo.sites'
	subprocess.call(call_string, shell=True)
	call_string = 'rm combo.vcf'
	subprocess.call(call_string, shell=True)

	
	call_string = 'mv ' + max_file + ' ./maxfiles/2'
	subprocess.call(call_string, shell=True)
	
	
    
    input_file.close()

