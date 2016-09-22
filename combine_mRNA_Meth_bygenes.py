#!/usr/bin/python2.7

import sys,math
from collections import defaultdict

mrna_seq=[]
meth_seq=[]
hash_mrna= {}
hash_meth={}


mRNA=open(sys.argv[1],"r")
Meth=open(sys.argv[2],"r")

##open(sys.argv[4],"w")


for line in mRNA.readlines():
    mrna=line.split('\t')
    mrna_seq.append(mrna)

for elem in xrange(1,len(mrna_seq),1):
	mrna_seq[elem][-1]=mrna_seq[elem][-1].rstrip('\n')

for line in Meth.readlines():
    meth=line.split('\t')
    meth_seq.append(meth)

for elem in xrange(1,len(meth_seq),1):
	meth_seq[elem][-1]=meth_seq[elem][-1].rstrip('\n')

for i in xrange(1,len(mrna_seq),1):
	hash_mrna[mrna_seq[i][1]]=mrna_seq[i]

for i in xrange(1,len(meth_seq),1):
	hash_meth[meth_seq[i][0]]=meth_seq[i]





"""
log2_id_mRNA=range(13,29,1)
log2FoldChange_LY_2h_vs_control_2h 	 13	log2FoldChange_TGFb_2h_vs_control_2h    14	log2FoldChange_LY_TGFb_2h_vs_control_2h		15
log2FoldChange_LY_2h_vs_TGFb_2h		 16	log2FoldChange_LY_2h_vs_LY_TGFb_2h      17	log2FoldChange_TGFb_2h_vs_LY_TGFb_2h		18	
log2FoldChange_LY_48h_vs_control_48h 	 19	log2FoldChange_TGFb_48h_vs_control_48h	20	log2FoldChange_LY_TGFb_48h_vs_control_48h	21	
log2FoldChange_LY_48h_vs_TGFb_48h 	 22	log2FoldChange_LY_48h_vs_LY_TGFb_48h    23	log2FoldChange_TGFb_48h_vs_LY_TGFb_48h		24
log2FoldChange_LY_48h_vs_LY_2h		 25	log2FoldChange_TGFb_48h_vs_TGFb_2h	26 	log2FoldChange_LY_TGFb_48h_vs_LY_TGFb_2h	27
log2FoldChange_control_48h_vs_control_2h 28
"""
## 13,34 ; 14,40; 15,44 ; 16,50; 17,52; 18,38; 19,28; 20,56; 
##21,46; 22,36; 23,54;  24,32; 25,42 ; 26,48; 27, 30; 28,26 ##-1
"""
log2_id_meth=range(26,57,2)
"log2FoldChange_Control2h_Control48h" 	26		"log2FoldChange_LY48h_Control48h" 		28	"log2FoldChange_LY+TGFb2h_LY+TGFb48h"	30	
"log2FoldChange_LY+TGFb48h_TGFb48h"   	32		"log2FoldChange_LY2h_Control2h"	  		34	"log2FoldChange_LY48h_TGFb48h"	        36	
"log2FoldChange_LY+TGFb2h_TGFb2h"	38		"log2FoldChange_TGFb2h_Control2h_FC_PVAL"	40	"log2FoldChange_LY2h_LY48h" 		42		
"log2FoldChange_LY+TGFb2h_Control2h"	44		"log2FoldChange_LY+TGFb48h_Control48h"		46	
"log2FoldChange_TGFb2h_TGFb48h"		48		"log2FoldChange_LY2h_TGFb2h"			50	
"log2FoldChange_LY+TGFb2h_LY2h" 	52		"log2FoldChange_LY+TGFb48h_LY48h"		54	"log2FoldChange_TGFb48h_Control48h" 	56	

"""
temp=[]

for k in hash_mrna.keys():
	for l in hash_meth.keys():
		sp=hash_meth[l][0].split()		
		if k==sp[4]:
			ls=[k,sp[3],hash_mrna[k][24],hash_meth[l][41]]
			temp.append(ls)
			##print temp			
			##print "\t".join(ls) 
			

##[value for value in temp if not math.isnan(value)] 

d = defaultdict(int)
d1 = defaultdict(int)


for y in xrange(1, len(temp),1):
        
        try:
            geneid = str(temp[y][0])
	    mr1 = float(temp[y][2]) 
            me2 = float(temp[y][3])
        except ValueError:
            continue
        try:
		d[geneid] += mr1
		d1[geneid] += me2
		
	except:
		pass

for k in d.keys():
	for k1 in d1.keys():
		if (k==k1):
			print k+'\t'+'UTR5'+'\t'+str(d[k])+'\t'+str(d1[k])

