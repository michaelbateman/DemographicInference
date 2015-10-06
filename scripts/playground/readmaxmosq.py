#pop_list = [ 'ABM','ABS','AC', 'AJ','AK','AN', 'AR', 'AS', 'AV'  ] 
pop_list = [ 'ABM','ABS','AC', 'AJ','AN', 'AR', 'AS', 'AV'  ] 
print pop_list

import numpy as np
from scipy import stats

pop_dict = {'ABM':'Burkina Faso M','ABS':'Burkina Faso S','AC':'Uganda', 'AJ':'Guinea-Bissau','AK':'Kenya','AN':'Cameroon', 'AR':'Angola', 'AS':'Gabon', 'AV':'Guinea'}

deciles_dict = {}
pop_dict = {}
local_index_dict = {}
len_dict = {}
#start_dict = {}
midpoint_list = []

import sys
num_samples = int(sys.argv[1])
num_cohorts = int(sys.argv[2])

#starting_tile = int(sys.argv[3])

tile_range = []

#for j in range(0,9):
    #tile_range.append(j*1 + starting_tile)
#print tile_range

k_total = 0  # k_total is never reset

for pop in pop_list:
    k_pop = 0 #k_pop is reset for each new population
    s = pop + '.samples'
    
   # local_index_dict.update({k_total:k_pop})
    try:
        input_file = open(s, "r")
        print "input file is", input_file.name
    except IOError:
    	print "Oh no, we couldn't open our file!"
    
 #   print "hello"
    for m in range(0, num_cohorts):
	k_pop+=1  # k_pop is used to count how many samples in each population
	#local_index_dict.update({k_total:k_pop})
	
	#sam1 = line.strip()
	#try:
	    #for x in range(2,num_samples+1):	
                #temp = input_file.next()
	        #u = temp.strip()
            ##temp2 = input_file.next()
	    ##temp3 = input_file.next()
	    ##temp4 = input_file.next()
	    ##temp5 = input_file.next()
	    ##temp6 = input_file.next()
	    ##temp7 = input_file.next()
	    ##temp8 = input_file.next()
	    ##sam2 = temp2.strip()
	    ##sam3 = temp3.strip()
	    ##sam4 = temp4.strip()
	    ##sam5 = temp5.strip()
	    ##sam6 = temp6.strip()
	    ##sam7 = temp7.strip()
	    ##sam8 = temp8.strip()
	#except StopIteration:
	    #pass
	
	#pop_dict.update({sam: pop})
	#print "hello2"
	#print sam.rstrip(), 2, 3
	
	region_size = 5000000
	num_regions = 50e6 / region_size
	filenames = []
	for i in range(0,int(num_regions)):
	    left = i * region_size
    	    right = (i+1) * region_size
            window = str(left) + '-' + str(right)
	    temp_max_file = './maxfiles/all/'+str(num_samples)+'/'+str(num_samples)+'.' + str(m) +'.'+pop + '.' + window+ '.max'
	    filenames.append(temp_max_file)
	
	t = './maxfiles/all/'+str(num_samples)+'/'+str(num_samples)+'.' + str(m) +'.'+pop + '.max'
	
	
	with open(t, 'w') as outfile:
    	    for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
	
	
    	try:
            max_file = open(t, "r")
            #print "input file is", max_file.name
        except IOError:
	    print "Oh no, we couldn't open our file:", max_file.name
	#print "hello4"
	vec=[]	
	for line0 in max_file:
	    tempnum = line0.rstrip()
	    num = int(tempnum)
	    vec.append(num)
	    
	#print 'There were', len(vec), 'matches collected'
#	print sam1,pop, k_pop, ' 25th = ', stats.scoreatpercentile(vec, 25), ' 50th = ', stats.scoreatpercentile(vec, 50),' 75th = ', stats.scoreatpercentile(vec, 75)
	if len(vec)>0:
	    
	    deciles = []	
	    for i in [10,20,30,40,50,60,70,80,90]:
	   # for i in tile_range:
	        hello = 3
		#print '%d th percentile = ' %i, stats.scoreatpercentile(vec, i)
		weight = 1
		deciles.append(weight *weight * stats.scoreatpercentile(vec, i))
	    deciles_dict.update({k_total: deciles})	
	    #print '25th = ', stats.scoreatpercentile(vec, 25)
	    #print '50th = ', stats.scoreatpercentile(vec, 50)
	    #print '75th = ', stats.scoreatpercentile(vec, 75)
	    #print 'The mean match length for', sam.strip(), 'from', pop, 'is', np.mean(vec)
	    #a = np.mean(vec)
	else:
	    pass
	
	pop_dict.update({k_total:pop})
	k_total += 1  #k_total is the index for all samples
    
    len_dict.update({pop:k_pop})
    midpoint_list.append(k_total - (k_pop) / 2 )
    input_file.close()
##print len_dict
#print midpoint_list
#print pop_list

#print deciles_dict.keys()
	
#for sam1 in deciles_dict.keys():
    #for sam2 in deciles_dict.keys():
	#diff = np.sqrt( np.sum( np.square( np.subtract( deciles_dict[sam1], deciles_dict[sam2] ) ) ) )
	#print "diff for", pop_dict[sam1], pop_dict[sam2], "=", diff

total_length = k_total

print 'Total length is ', total_length
#print len_dict

A = len(deciles_dict.keys())
M = np.zeros((A,A))
X = np.zeros((A,1))
print np.shape(M)
#print deciles_dict

i = 0
j = 0
#k_total = 0  # k_total is never reset

bad = []

print pop_dict #local_index_dict[0]
print pop_dict[1] #local_index_dict[1]
for i in range(0, total_length):
   #print i
    for j in range(0, total_length):
	#print j
	if( i <= j):
		
	    diff =  np.sqrt( np.sum( np.square( np.subtract( deciles_dict[i], deciles_dict[j] ) ) ) )
		        ##print i, j, A
	    M[i][j] = diff
	                ##print "diff for", pop1, pop2, "=", diff
			
	else:
	    M[i][j] = M[j][i]
	    
    X[i] = np.sum(M[i][:]) / A
    #print pop_dict[i], local_index_dict[i],'avg = ', X[i]
    print pop_dict[i], 'avg = ', X[i]
    if ( X[i] > 10000):# and pop_dict[i] !='FIN' and pop_dict[i] !='PEL'):
	print pop_dict[i], i, 'avg = ', X[i]
	bad.append(i)
	asdf = 1
    #print i
bad
temp_new_M = np.delete(M, bad, axis =0)
new_M = np.delete(temp_new_M, bad, axis =1)


num_cohorts


midpoint_list = []

i = 0
for pop in pop_list:
    midpoint_list.append(num_cohorts * (i + .5) )
    i += 1

import matplotlib.pyplot as plt
plt.matshow(new_M)
plt.summer()
#x = midpoint_dict.values()
#names = midpoint_dict.keys()
plt.xticks(midpoint_list, pop_list)
plt.yticks(midpoint_list, pop_list)
title = ' L2 distance between deciles of length dist, '+str(2*num_samples)+' haps = '+str(num_samples)+' ind'
plt.title(title)
#plt.subplot(121)

import pylab
#figure_name = 'figure.ag1000g.'+str(num_samples)+'.'+str(num_cohorts)+'.png'
#pylab.savefig(figure_name)
#figure_name = 'figure.ag1000g.'+str(num_samples)+'.'+str(num_cohorts)+'.pdf'
#pylab.savefig(figure_name)
plt.show()




