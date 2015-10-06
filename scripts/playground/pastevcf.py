import subprocess

in_file = 'ag1000g.phase1.AR2.3L.PASS.vcf.gz'


region_size = 1000000
num_regions = 50e6 / region_size

print num_regions
print int(num_regions)

list_string = ''


for i in range(0,10):
    	left = i * region_size
	right = (i+1) * region_size
	window = str(left) + '-' + str(right)
		
	
	
	
	phased_name = 'ag1000g.phase1.AR2.3L.PASS.' + window + '.phased.vcf.gz'
	
	list_string = list_string + ' ' + phased_name
	
	
	

call_string = 'time ../bcftools/bcftools concat -O z -o ag1000g.phase1.AR2.3L.PASS.10Mb.0.phased.vcf.gz ' + list_string
print call_string
subprocess.call(call_string, shell=True)
