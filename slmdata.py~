#!/usr/bin/python2.7

import sys,os

def addhash(g, obj):
    ghash_slimmed[g]=ghash_slimmed.get(g,[])+[obj]

def writeDict(dict, filename, sep):
	with open(filename, "w") as f:
        	for i in dict.keys():            
            		f.write(i + "\t" + sep.join([str(x) for x in dict[i]]) + "\n")


file_to_slim=[]
ghash={}
ghash_slimmed={}

f=open(sys.argv[1],'rU')

ind=sys.argv[3]
for line in f.readlines():
    	fsp=line.split('\t')
	fsp[-1]=fsp[-1].rstrip('\n')
	file_to_slim.append(fsp)

for i in xrange(0,len(file_to_slim),1):
	ghash[file_to_slim[i][ind]]=file_to_slim[i]

for k in ghash.keys():
	for i in xrange(0,len(file_to_slim),1):
		if (k==file_to_slim[i][sys.argv[3]]):
			addhash(k,file_to_slim[i][sys.argv[4]])

writeDict(ghash_slimmed, sys.argv[2], "\t")




