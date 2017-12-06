from chrm_state import *
import sys,math

query=sequence[int(inx[2][1]):int(inx[2][2])+1]

Vit= [ [ 0 for i in range(len(query)) ] for j in range(len(code)+1) ]

BTrk= [ [ 0 for i in range(len(query)) ] for j in range(len(code)+1) ]

state_predict=[0 for i in range(len(query))]

for s in xrange(0,len(Vit),1):
    Vit[s][0]=math.log(nI[s])+math.log(nE[s][Sym.index(query[0])],10)
    BTrk[s][0]=0



for col in xrange(1,len(query),1):
    for row in range(len(code)+1):
        poss=[Vit[0][col-1]+math.log(nS[0][row],10),Vit[1][col-1]+math.log(nS[1][row],10)]
        max_state=max(poss)
        Vit[row][col]=max_state+math.log(nE[row][Sym.index(query[col])],10)
        BTrk[row][col]=poss.index(max_state)
        
        
max_last=max([r[-1] for r in Vit])

state_predict[-1]=[r[-1] for r in Vit].index(max_last)


for b in xrange(len(query)-2,-1,-1):
    state_predict[b]=BTrk[state_predict[b+1]][b+1]
    
print state_predict   
    
    
    
