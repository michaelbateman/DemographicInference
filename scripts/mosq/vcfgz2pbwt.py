import subprocess

#pop_list = [ 'ABM','ABS','AC', 'AJ','AK','AN', 'AR', 'AS', 'AV'  ] 
#pop_list = [ 'AB','AC', 'AJ','AK',
#pop_list = ['AN', 'AR', 'AS', 'AV'  ] 

pop_list = [ 'ABM','ABS']

chrom = str(3) + 'L'


# pop_file is outside the pop loop to fix ABM and ABS

big_pop_file = 'ag1000g.AB.phase1.AR2.' + chrom + '.PASS.vcf' 
#call_string = 'time gunzip ' + big_pop_file + '.gz'
#print call_string
#subprocess.call(call_string, shell=True)
    
for pop in pop_list:

    
    pop_file = 'ag1000g.' + pop + '.phase1.AR2.' + chrom + '.PASS.vcf' 

    sample_file = pop + '.samples'
    call_string = '../../pbwtCurrent/pbwt/pbwt -readAll ' + big_pop_file + ' -selectSamples ' + sample_file + ' -writeAll ' + pop_file
    print call_string
    subprocess.call(call_string, shell=True)


	

	
    #call_string = '../../pbwtCurrent/pbwt/pbwt -readVcfGT ' + pop_file +' -writeAll ' + pop_file
    #print call_string
    #subprocess.call(call_string, shell=True)
    

