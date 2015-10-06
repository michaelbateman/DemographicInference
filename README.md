# Demographic_Inference
This repo contains a program to carry out demographic inference from genotyped or haplotyped dna.  (More recent events can be inferred with haplotypes.)  

FIRE is the name of the program that does the inference; it is written in C.  Pre/post-processing scripts are written in python.

Some of the scripts are specific to either human or (anopheles) mosquito samples.  So far we have run inference on humans and mosquitoes, but in princpile the method should apply more broadly.  The human data we have used is from the 1000 Genomes Project (http://www.1000genomes.org/); the mosquito data was from the Anopheles Gambiae 1000 Genomes Project (https://www.malariagen.net/projects/vector/ag1000g). 
