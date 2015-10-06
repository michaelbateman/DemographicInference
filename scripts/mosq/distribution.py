

#my_list = [ 'AK','ABM','ABS','AC', 'AJ','AN', 'AR', 'AS', 'AV'  ] 

my_list = [ 'ABM','ABS','AC', 'AJ','AN', 'AV'  ] 

#my_list = [ 'ABS','AC','AR', 'AK','AS'] 

for pop in my_list:
    f = pop + '.samples'
    num_lines = sum(1 for line in open(f))
    print pop, 'has', num_lines, 'samples'


#my_list = [    'AK', 'ABS' ] 



pop_dict = {'ABM':'Burkina Faso M','ABS':'Burkina Faso S','AC':'Uganda', 'AJ':'Guinea-Bissau','AK':'Kenya','AN':'Cameroon', 'AR':'Angola', 'AS':'Gabon', 'AV':'Guinea'}

#my_list = [ 'AS', 'ABS'	 ] 



print pop_dict['ABS']

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


dist_dict = {}

local_index_dict = {}
len_dict = {}
#start_dict = {}
midpoint_list = []

import sys
#num_samples = int(sys.argv[1])
num_cohorts = int(sys.argv[2])


num_samples_list = [int(sys.argv[1])]

num_samples_list = [.5, 1, 2, 4, 8, 16, 30]

num_samples_list = [1]

#starting_tile = int(sys.argv[3])

tile_range = [5,10,12, 15, 17,20,21]

#tile_range = []

#for j in range(0,9):
    #tile_range.append(j*1 + starting_tile)
#print tile_range
chrom = str(3) + 'L'

plotting_list = ['PEL']#, 'PEL']
print num_samples_list
print my_list
for num_samples in num_samples_list:
    plot_num = 1
    for pop in my_list:
	
	print pop
	my_s = pop + '.' + pop
	vec=[]	
	for i in range(0, num_cohorts):
	        pop_file = 'ag1000g.' + pop + '.phase1.AR2.' + chrom + '.PASS.vcf' 
		max_file = pop_file + '.' + str(num_samples) +'.'+ str(i) + '.max'
	
		t = './maxfiles/all/'+str(int(num_samples))+'/'+ max_file
		#print t
		try:
		    max_file = open(t, "r")
		    #print "input file is", max_file.name
		except IOError:
		    print "Oh no, we couldn't open our file:", max_file.name
		
		
		for line0 in max_file:
		    tempnum = line0.rstrip()
		    num = int(tempnum)
		    if (num >=0):
			#vec.append(np.sqrt(1 + np.log(1 +num) ) )
			vec.append(np.log(1 +num))
			#vec.append(num)
		    else:
			pass		
		if len(vec)>0:
		    dist_dict.update({pop: vec})	
		    #print len(vec)
		else:
		    pass
		
	print 'There were a total of ', len(vec), ' matches collected for', pop
	#plt.subplot(2, 1, plot_num)
	#(m, mbins, mpatches) = plt.hist(dist_dict[pop], bins = 100, alpha = .25, label ='adsf', cumulative = True, normed = True)
	#(n, bins, patches) = plt.hist(dist_dict[pop], bins = 100, alpha = .25, label =pop_dict[pop])#, normed = True)
	#print bins
	#print np.exp(bins)
	#print n
	
	#A = min(vec)
        #B = max(vec)
	##(n, BINS, patches) = plt.hist(vec, bins=range(A, B + 1, 1),  alpha = .25, label = pop)
	#bin = np.logspace(0, np.log10(B), 100)  
	#plt.clf()
	#Y = np.zeros(( B - A,1))
	#for i in range(A, B):
            #print A
            #print B
            #print np.shape(Y)
	    #print np.shape(n)
	    #print i
	    #Y[i] = (i+1) * n[i]
	#(M, MBINS, Mpatches) = plt.hist(Y, bin , alpha = .25, label = 'pop', normed = True)
	#plot_num +=1
	#plt.gca().set_xscale("log")
	
	(n, bins) = np.histogram(vec, bins = 100, normed = True)
        Y = np.zeros((100,))
        X = np.zeros((100,))
        delta = np.zeros((100,))
        mass = np.zeros((100,))

    
        total = 0
    
        for i in range(0, 100):
       	    Y[i] = np.exp(bins[i]) * n[i]
	    X[i] = bins[i]
	    delta[i] = bins[i+1] - bins[i]
            mass[i] = Y[i] * delta[i]
	    total +=mass[i]
	
    
        print 'Total = ', total
        plt.plot(X, np.divide(Y, total*1.0), label = pop_dict[pop])


    
    title = str(num_cohorts) + ' samples per country'
    plt.title(title)
    plt.ylabel('portion of genome spent in ROH of given length')
    plt.xlabel('log(length)')
    plt.legend(loc='upper right')

    import pylab

    my_string =   str(num_cohorts) + '.samples'

    for pop in my_list:
        my_string = my_string  + '.' + pop 

    my_string = my_string +  '.png'
    print my_string
    pylab.savefig(my_string)
    #plt.close("all")



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
plt.show()




