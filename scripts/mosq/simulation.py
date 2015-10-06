import subprocess
import numpy as np
from scipy import stats
import subprocess
import sys
import pylab
import matplotlib.pyplot as plt



num_haps = 2
chr_length = 1e8




#print mu, mu * pow(10,9)
#print rho, rho * pow(10,9)
#print t
#print r

het_dict = {}

scenario_list = ['flat', 'crash']

fig_dict = {}

for scenario in scenario_list:
    mu = 3.4e-9
    rho = 1.1e-8

    pop_size = 1e6
    t = 4 * mu * pop_size 
    r = 4 * rho * pop_size

    het_counter = 0
	

    if scenario == 'flat':
	print 'Pop size = ', pop_size
	print 't = ', t
	print 'r = ', r
	call_string = '../pbwtMay/pbwt/macs/macs ' + str(num_haps) +  ' '  + str(chr_length) + ' -t ' +str(t) + ' -r ' + str(r) + ' > simulation.macs 2>/dev/null' 
	print call_string
	subprocess.call(call_string, shell=True)

    elif scenario == 'crash':
	
	factor = .01
	new_pop_size = pop_size * factor

	t = 4 * mu * new_pop_size 
	r = 4 * rho * new_pop_size
	print 'new pop size = ', new_pop_size
	print 't = ', t
	print 'r = ', r
	call_string = '../pbwtMay/pbwt/macs/macs ' + str(num_haps) +  ' '  + str(chr_length) + ' -t ' +str(t) + ' -r ' + str(r) + ' -eN .5 ' + str(1/factor) + ' > simulation.macs 2>/dev/null' 
	print call_string
	subprocess.call(call_string, shell=True)




    call_string = '../pbwtCurrent/pbwt/pbwt -readMacs simulation.macs -write simulation.pbwt -writeSites simulation.sites -maxWithin >  simulation.max' 
    print call_string
    subprocess.call(call_string, shell=True)

    max_file = 'simulation.max'
	
    t = max_file
    print t
    try:
        max_file = open(t, "r")
        #print "input file is", max_file.name
    except IOError:
        print "Oh no, we couldn't open our file:", max_file.name

    vec = []	
    for line0 in max_file:
        tempnum = line0.rstrip()
        num = int(tempnum)
        if (num >=0):
            #vec.append(np.log(1+np.log(1+num)))
	    vec.append( np.log(1+num) )
	    #vec.append( num )
	    if num == 0:
		het_counter+=1
   	    else:
		pass
        else:
            a = 4
	    #het_counter +=1	
	    #vec.append( np.log(1+num) )
    #if len(vec)>0:
        #dist_dict.update({pop: vec})	
        ##print len(vec)
    #else:
        #pass
    #print len(vec)
    #print vec
    het_dict.update({scenario:het_counter})
    A = min(vec)
    B = max(vec)
    #(n, bins, patches) = plt.hist(vec, bins = 100, alpha = .25, label = scenario, normed = True)

    #(n, BINS, patches) = plt.hist(vec, bins=range(A, B + 1, 1),  alpha = .25, label = scenario)
    (a,) = np.shape(n)
    #(test, ) = bins
    #print bins
    #Y = np.zeros(( B - A,1))
    
    (n, bins) = np.histogram(vec, bins = 100, normed = True)
    Y = np.zeros((100,))
    X = np.zeros((100,))
    delta = np.zeros((100,))
    mass = np.zeros((100,))
    #plt.clf()
    #print n
    #print len(vec)
    #print np.sum(n)
    #print A
    #print B
    
    total = 0
    
    for i in range(0, 100):
        #print A
        #print B
        #print np.shape(Y)
	#print np.shape(n)
	#print i
	Y[i] = np.exp(bins[i]) * n[i]
	#print Y[i]
	X[i] = bins[i]
	delta[i] = bins[i+1] - bins[i]
        mass[i] = Y[i] * delta[i]
	total +=mass[i]
	
    
    print 'Total = ', total
    #print Y
    
    #bin = np.logspace(0, np.log10(B), 100)  
    #np.insert(bin, 0,0.5)
    #np.insert(bin, 0,0)
    #bin = [0,1000,1000000000]
    
    plt.plot(X, np.divide(Y, total*1.0), label = scenario)
    #pic_name = scenario + '_fig'
    #fig_dict.update({scenario:plt.figure()})
    #plt.plot(X,  Y, label = scenario)
    #new_bins = bins
    
    #print Y
    #print X
    #(M, MBINS, Mpatches) = plt.hist(Y, X, alpha = .25, label = 'new', normed = True)
    #plt.hist(Y, bins=100 , alpha = .25, label = 'new', normed = True)
    #print np.logspace(0.0, np.log10(B), 100)
    #print M
    
    #fig = plt.figure()
    #plt.gca().set_xscale("log")
    #plt.close('all')
    #np.delete(bins, -1)
    ##n.append(0)
    #print np.shape(bins)
    #print np.shape(n)
    #Y = np.zeros((100,1))
    #X = np.zeros((100,1))
    #for i in range(0, 100):
	#Y[i] = n[i]
        #Y[i] = np.multiply( n[i], np.exp(bins[i]) )
	##Y[i] = np.multiply( n[i], bins[i] )
	#X[i] = bins[i]
    #print np.shape(Y)
    #print np.shape(X)
    #print X
    #print Y
    #plt.hist(n, bins , alpha = .25, label = scenario)
    #plt.bar( X, Y, label = scenario)




#for scenario in fig_dict:
    #plt.show(fig_dict[scenario])


#call_string = '../pbwtMay/pbwt/pbwt -readMacs simulation.macs -write simulation.pbwt -writeSites simulation.sites -maxWithin >  simulation.max' 
#print call_string
#subprocess.call(call_string, shell=True)

#max_file = 'simulation.max'
	
#t = max_file
#print t
#try:
    #max_file = open(t, "r")
    ##print "input file is", max_file.name
#except IOError:
    #print "Oh no, we couldn't open our file:", max_file.name

#vec = []	
#for line0 in max_file:
    #tempnum = line0.rstrip()
    #num = int(tempnum)
    #if (num >=0):
	#vec.append(np.log(1+num))
    #else:
	#pass		
##if len(vec)>0:
    ##dist_dict.update({pop: vec})	
    ###print len(vec)
##else:
    ##pass
#print len(vec)

#plt.hist(vec, bins = 100, alpha = .25, label = 'simulation', normed = True)







for scenario in scenario_list:
    print 'There were', het_dict[scenario], 'heterozygous sites in the', scenario, 'scenario.'
    print 'HELLOOOOO   Since there are', chr_length, 'sites, this implies an average heterozygosity of', het_dict[scenario] / float(chr_length)
    		










title = 'Portion of genome spent in homozygous stretches of given length'
plt.title(title)
plt.legend(loc='upper right')
plt.ylabel('relative frequency')
plt.xlabel('log(length of ROH)')
pylab.savefig('simulation.png')
plt.show()