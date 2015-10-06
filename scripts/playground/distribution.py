try:
    input_file = open("1kgPOPSsortednew.list", "r")
    print "Input file is", input_file.name
except IOError:
    print "Oh no, we couldn't open our file!"

pop_list = []

for line in input_file:
    pop_list.append(line.strip() )

#pop_list.remove('ASW')
#pop_list.remove('MXL')

#pop_list = [ 'GWD', 'ESN', 'YRI', 'MSL','LWK', 'ACB', 'ASW'  ] #

africa = [ 'GWD', 'ESN', 'YRI', 'MSL','LWK', 'ACB', 'ASW'  ] 
eastasia = ['CHS', 'CHB', 'CDX', 'KHV', 'JPT']
southasia = ['ITU', 'STU', 'BEB', 'PJL','GIH']
europe = ['CEU', 'GBR', 'FIN', 'TSI', 'IBS']
america = ['MXL', 'PEL', 'CLM', 'PUR']

my_list = ['GBR', 'FIN', 'PEL', 'ESN', 'GWD', 'CEU']

my_list = [ 'ESN', 'GWD']





print pop_list
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


dist_dict = {}
pop_dict = {}
local_index_dict = {}
len_dict = {}
#start_dict = {}
midpoint_list = []

import sys
#num_samples = int(sys.argv[1])
num_cohorts = int(sys.argv[2])


num_samples_list = [int(sys.argv[1])]

num_samples_list = [.5, 1, 2, 4, 8, 16, 30]

#starting_tile = int(sys.argv[3])

tile_range = [5,10,12, 15, 17,20,21]

#tile_range = []

#for j in range(0,9):
    #tile_range.append(j*1 + starting_tile)
#print tile_range

i=1
plotting_list = ['PEL']#, 'PEL']
print num_samples_list
print my_list
for num_samples in num_samples_list:
	#for pop in pop_list:
    for pop in my_list:
	print pop
	my_s = pop + '.' + pop
	vec=[]	
	for m in range(0, num_cohorts):
	
	
		t = './maxfiles/all/'+str(int(2*num_samples))+'/'+str(int(2*num_samples))+'.' + str(m) +'.'+pop + '.max'
		print t
		try:
		    max_file = open(t, "r")
		    #print "input file is", max_file.name
		except IOError:
		    print "Oh no, we couldn't open our file:", max_file.name
		
		
		for line0 in max_file:
		    tempnum = line0.rstrip()
		    num = int(tempnum)
		    if (num >0):
			vec.append(np.log(num))
		    else:
			pass		
		if len(vec)>0:
		    dist_dict.update({my_s: vec})	
		    #print len(vec)
		else:
		    pass
	#plt.subplot(3, 1, i)
	if (5==5):
	    plt.hist(dist_dict[my_s], bins = 100, alpha = .25, label =my_s, normed = True)
	else:
	    pass
	i+=1
	groupA = ['PEL', 'FIN']
	groupB = ['PEL', 'FIN']
	
	
	
    for pop in my_list:
	    print pop
	    for pop2 in my_list:
 	       
		my_s = pop+ '.'+pop2
	        if (pop2> pop):    
		    vec=[]
		    for m in range(0, num_cohorts):
			sample_pair_root = str(num_samples)+'.'+str(m) +'.'+ pop + '.' + pop2
			t = './maxfiles/cross/'+str(num_samples) + '/' + sample_pair_root + '.max'
			
			print t
			try:
			    a = 3
			#    max_file = open(t, "r")
			    #   print "input file is", max_file.name
			except IOError:
			    print "Oh no, we couldn't open our file:", max_file.name
				
			for line0 in max_file:
			    tempnum = line0.rstrip()
			    num = int(tempnum)
			    if (num >0):
				vec.append(np.log(num))
			    else:
				pass	
				
			if len(vec)>0:
			    dist_dict.update({my_s: vec})	
			
			else:
			    pass
		
		    print 'asdf'
		    #plt.subplot(3, 1, i)
		    #plt.hist(dist_dict[my_s], bins = 100, alpha = .25, label =my_s, normed = True)
		    i+=1
		else:
		    pass
	
	
    		
	


#i = 1
	#pairs = len(my_list) * (len(my_list) - 1) / 2 + len(my_list)
	#for pop in my_list:
	#for pop2 in my_list:
		#if (pop2>= pop):    	    
		#my_s = pop+ '.'+pop2
		#plt.subplot( pairs, 1, i)
		#plt.legend(loc='upper left')
		#plt.ylabel('number of matches')
		#if i == pairs:
			#plt.xlabel('log(match length)')
		#else:
			#pass    
		#title = 'Cross matches between ' + pop + ' and ' + pop2
		#plt.title(title)
		
		#plt.hist(dist_dict[my_s], bins = 100, alpha = .25, label =my_s, normed = True)
		#i+=1
		#else:
		#pass
	
	
	
	#for pop in my_list:
	    #for pop2 in my_list:
		#if (pop2>= pop):    	    
		    #my_s = pop+ '.'+pop2
		    #plt.hist(dist_dict[my_s], bins = 100, alpha = .25, label =my_s, normed = True)
		
		#else:
		    #pass
    
    title = 'Internal matches, ' + str(int(2 * num_samples)) + ' individuals from each population'
    plt.title(title)
    plt.legend(loc='upper left')
    plt.ylabel('relative frequency')
    plt.xlabel('log(match length)')

    import pylab

    my_string =   str(num_samples) + '.internal'

    for pop in my_list:
        my_string = my_string  + '.' + pop 

    my_string = my_string +  '.png'
    print my_string
    pylab.savefig(my_string)
    plt.close("all")



#groupA = ['CHS', 'CHB']
#groupB = ['CHS', 'CHB']
##groupB = ['ESN', 'YRI']


#for pop in groupA:
    #print groupA
    #for pop2 in groupB:
	#my_s = pop+ '.'+pop2
	#if (pop2> pop):    
	    #vec=[]
	    #for m in range(0, num_cohorts):
		#sample_pair_root = str(num_samples)+'.'+str(m) +'.'+ pop + '.' + pop2
		#t = './maxfiles/pairs/'+str(num_samples) + '/' + sample_pair_root + '.max'
		
		##print t
		#try:
		    #max_file = open(t, "r")
		 ##   print "input file is", max_file.name
		#except IOError:
		    #print "Oh no, we couldn't open our file:", max_file.name
			
		#for line0 in max_file:
		    #tempnum = line0.rstrip()
		    #num = int(tempnum)
		    #if (num >0):
	        	#vec.append(np.log(num))
	    	    #else:
			#pass	
		    	
		#if len(vec)>0:
	    	    #dist_dict.update({my_s: vec})	
	    	    #print len(vec), 'hello'
		#else:
	    	    #pass
	
		
		
	#else:
	    #pass



#print dist_dict


#for pop in groupA:
    #print groupA
    #for pop2 in groupB:
        #if (pop2> pop):    	    
	    #my_s = pop+ '.'+pop2
	    #plt.hist(dist_dict[my_s], bins = 100, alpha = .25, label =my_s)
	#else:
	    #pass


#for pop in plotting_list:
#for pop in eastasia:
    ##plt.hist(dist_dict['CHS'], bins = 100, alpha = .25, label ='CHS')
    ##plt.hist(dist_dict['JPT'], bins = 100, alpha = .25)
    #plt.hist(dist_dict[pop], bins = 100, alpha = .1, label =pop)

#for pop in africa:
    ##plt.hist(dist_dict['CHS'], bins = 100, alpha = .25, label ='CHS')
    ##plt.hist(dist_dict['JPT'], bins = 100, alpha = .25)
    #plt.hist(dist_dict[pop], bins = 100, alpha = .1, label =pop)





#plt.show()


#plt.matshow(new_M)
#plt.summer()
##x = midpoint_dict.values()
##names = midpoint_dict.keys()
#plt.xticks(midpoint_list, pop_list)
#plt.yticks(midpoint_list, pop_list)
#title = ' L2 distance between deciles of length dist, '+str(2*num_samples)+' haps = '+str(num_samples)+' ind'
#plt.title(title)
##plt.subplot(121)
#plt.show()




