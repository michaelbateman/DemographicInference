s = 'ag1000g.phase1.AR2.samples.txt'
print s
try:
	input_file = open(s, "r")
	print "input file is", input_file.name
except IOError:
	print "Oh no, we couldn't open our file!"#
	

with open("ABM.samples", "w") as output:
   

    for line in input_file:
	vec = line.split('\t')
	sample = vec[0]
	country = vec[4]
	subspecies = vec[9]
	print sample, '\t', country,'\t',subspecies
	
	if country == 'Burkina Faso' and subspecies == 'M':
            output.write(sample + "\n")
	