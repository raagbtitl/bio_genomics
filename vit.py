from chrm_state import *
import sys,math

query=sequence[int(inx[2][0]):int(inx[2][1])+1]

Vit= [ [ 0 for i in range(len(query)) ] for j in range(2) ]

BTrk= [ [ 0 for i in range(len(query)) ] for j in range(2) ]

state_predict=[0 for i in range(len(query))]

Vit[0][0]=math.log(nE[0][Sym.index(query[0])],10)

Vit[1][0]=math.log(nE[1][Sym.index(query[0])],10)

BTrk[0][0]=0
BTrk[1][0]=0


for col in xrange(1,len(query),1):
    for row in range(2):
        poss=[Vit[0][col-1]+math.log(nS[0][row],10),Vit[1][col-1]+math.log(nS[1][row],10)]
        max_state=max(poss)
        Vit[row][col]=max_state+math.log(nE[row][Sym.index(query[col])],10)
        BTrk[row][col]=poss.index(max_state)
        
        
max_last=max([r[-1] for r in Vit])

state_predict[-1]=[r[-1] for r in Vit].index(max_last)


for b in xrange(len(query)-2,-1,-1):
    state_predict[b]=BTrk[state_predict[b+1]][b+1]
    
print state_predict   
    
    
    
