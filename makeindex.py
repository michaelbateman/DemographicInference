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


import subprocess
#call_string = 'cd individuals'
#subprocess.call(call_string, shell=True)



for pop in pop_list:
    
    s = pop + '.samples'
    print s
    try:
        input_file = open(s, "r")
        print "input file is", input_file.name
    except IOError:
    	print "Oh no, we couldn't open our file!"
    
 #   print "hello"
    for line in input_file:
        
        sam1 = line.strip()
	s1 = sam1+'.'+pop+'.hets.chr1.vcf.gz'
	v1 = sam1+'.'+pop+'.hets.chr1.vcf'
	
	call_string = 'gunzip ' + s1
	subprocess.call(call_string, shell=True)
	print call_string
	
	call_string = '../../tabix/bgzip ' + v1
	subprocess.call(call_string, shell=True)
	print call_string
	
	call_string = '../../bcftools/bcftools index ' + s1
	subprocess.call(call_string, shell=True)
   	print call_string 
	
    input_file.close()
#call_string = 'cd ..'
#subprocess.call(call_string, shell=True)
