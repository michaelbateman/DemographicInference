import sys
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import time
import pylab as P

from scipy import stats		


try:
    input_file = open("1kgPOPSsorted.list", "r")
    print "Input file is", input_file.name
except IOError:
    print "Oh no, we couldn't open our file!"

pop_list = []

for line in input_file:
    pop_list.append(line.strip() )




#pop_list = [ 'CEU', 'FIN', 'YRI']#,'MXL','TSI', 'IBS', 'GBR', 'FIN'  ] #

pop_list = ['GBR', 'PEL', 'ESN', 'CEU', 'FIN', 'YRI']

#pop_list = ['PEL', 'ESN']

print pop_list


#chrom_list = ['chr1', 'chr22', 'chr21', 'chr20', 'chr19', 'chr18', 'chr17', 'chr16', 'chr15', 'chr14', 'chr13', 'chr12']
chrom_list = [ 'chr11', 'chr10', 'chr9', 'chr8', 'chr7', 'chr6', 'chr5', 'chr4', 'chr3', 'chr2']

for chrom in chrom_list:
	fileroot = 'ALL.'+ chrom + '.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes'
	filename = 'ALL.'+ chrom + '.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf'
	
	special = 'ALL.chr1'
	call_string = 'wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/' + filename + '.gz'
	print call_string
	subprocess.call(call_string, shell=True)
	
	call_string = 'wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/' + filename + '.gz.tbi'
	print call_string
	subprocess.call(call_string, shell=True)
	
	call_string = 'time gunzip ' + filename + '.gz'
	print call_string
	subprocess.call(call_string, shell=True)
	
	
	
	call_string = '../pbwtDec/pbwt/pbwt -readVcfGT ' + filename + ' -writeAll '  + fileroot
	print call_string
	subprocess.call(call_string, shell=True)
	
	call_string = 'time gzip ' + filename
	print call_string
	subprocess.call(call_string, shell=True)
	
	for pop in pop_list:
		pop_root = pop + '.' + chrom
		if chrom == 'chr1':
		    call_string = '../pbwtDec/pbwt/pbwt -readAll ' + special + ' -selectSamples ' + pop + '.samples' + ' -writeAll ' + pop_root
		else:
		    call_string = '../pbwtDec/pbwt/pbwt -readAll ' + fileroot + ' -selectSamples ' + pop + '.samples' + ' -writeAll ' + pop_root
		print call_string
		subprocess.call(call_string, shell=True)


	