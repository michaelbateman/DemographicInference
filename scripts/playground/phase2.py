import subprocess

in_file = 'ag1000g.phase1.AR2.3L.PASS.vcf.gz'


region_size = 1000000
num_regions = 50e6 / region_size

print num_regions
print int(num_regions)



for i in range(7,int(num_regions)):
    if (i % 3 == 0) or (i % 3 == 1):
	print '	passing...'
    else:
	left = i * region_size
	right = (i+1) * region_size
	window = str(left) + '-' + str(right)
	out_file = 'ag1000g.phase1.AR2.3L.PASS.' + window + '.vcf'
	
	call_string = 'time ../bcftools/bcftools view -o ' + out_file + ' ' + in_file + ' -r 3L:' + window
	print call_string
	print 'Now creating the file:   ', out_file
	print '.....'
	subprocess.call(call_string, shell=True)
	
	
	call_string = '../pbwtCurrent/pbwt/pbwt -readVcfGT ' + out_file + ' -writeAll temp_file_name'
	print call_string
	print 'Now preparing site file...'
	subprocess.call(call_string, shell=True)
	
	
		# The 1530 just below is the number of haplotypes in 765 samples
		# Should change in different situation
	
	phased_name = 'ag1000g.phase1.AR2.3L.PASS.' + window + '.phased.vcf'
	call_string = '../pbwtCurrent/pbwt/pbwt -readVcfGT ' + out_file + ' -phase 1530 -readAll temp_file_name -writeVcf ' + phased_name
	print call_string
	print 'Now phasing...'
	subprocess.call(call_string, shell=True)
	
	call_string = 'time gzip ' + phased_name
	print call_string
	subprocess.call(call_string, shell=True)
	
	call_string = 'rm ' + out_file
	print call_string
	subprocess.call(call_string, shell=True)
	
	
	print 'Progress:  %d out of %d regions complete.' %(i+1, num_regions)

print call_string