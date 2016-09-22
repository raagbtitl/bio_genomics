#!/usr/bin/python2.7

import sys,re

gene=[]
ex=''
df_contr=[]
contr=0
##hippo_contr=[]
hash_contr={}
f=open("LPA_Hippo_targets_RA-1.1.csv",'r')

f_contr=open(sys.argv[1],'r')





for line in f.readlines():
    fs=line.split()
    gene.append(fs)



for line in f_contr.readlines():
    fs_contr=line.split('\t')
    df_contr.append(fs_contr)


for i in xrange(1,len(df_contr),1):
	hash_contr[df_contr[i][0]]=df_contr[i]

inside=[]
inside.append("gene")
inside.append(" ".join(map(str,df_contr[0])))
print " ".join(inside)	

for g in xrange(0,len(gene),1):	
	for k in hash_contr.keys():
		try:
			print " ".join(gene[g])+'\t'+" ".join(hash_contr[gene[g][0]])
		except:
			pass
		
