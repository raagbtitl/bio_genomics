#!/usr/bin/python2.7

import sys,math
##from collections import defaultdict


mrna_seq=[]
nonCod_gtf=[]
hash_mrna= {}
hash_nonCod_gtf={}


mRNA=open(sys.argv[1],"r")
nonC=open(sys.argv[2],"r")

##open(sys.argv[4],"w")


for line in mRNA.readlines():
    mrna=line.split(\t)
    mrna_seq.append(mrna)

for elem in xrange(1,len(mrna_seq),1):
	mrna_seq[elem][-1]=mrna_seq[elem][-1].rstrip(\n)

for line in nonC.readlines():
    nonc=line.split(\t)
    nonCod_gtf.append(nonc)

for elem in xrange(1,len(nonCod_gtf),1):
	nonCod_gtf[elem][-1]=nonCod_gtf[elem][-1].rstrip(\n)

for i in xrange(1,len(mrna_seq),1):
	hash_mrna[mrna_seq[i][0]]=mrna_seq[i]

for i in xrange(1,len(nonCod_gtf),1):
	hash_nonCod_gtf[nonCod_gtf[i][8]]=nonCod_gtf[i]

##print hash_mrna

for k in hash_nonCod_gtf.keys():
	sp=hash_nonCod_gtf[k][8].split(;)
	sp1=sp[4].split("gene_name")
	sp2=sp1[1].split(")
	##i=sp2[1]
	try:
		print hash_mrna[sp2[1]]
	except:
		pass			
		

	
