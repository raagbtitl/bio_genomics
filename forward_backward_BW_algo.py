




def forward_algo(T, E, I, Seq):

    Sym=["A","C","G","T", "N"]
    forward = [ [ 0 for i in range(len(Seq)) ] for j in range(len(T)) ]
    
       
    for s in xrange(0,len(forward),1):
        forward[s][0] = math.log(I[s],10) +math.log(E[s][Sym.index(Seq[0])],10)

    
    for col in xrange(1,len(Seq),1):
        for row in range(len(T)):
            for row1 in range(len(T)):
                first_T=first_T+ (math.log(forward[row][col-1],10)+ math.log(T[row1][row],10))
            second_E= first_T+math.log(E[row][Sym.index(Seq[col])],10)
            forward[row][col]=second_E
            first_T=0
            

    return(forward)



def backward_algo(T, E, I, Seq):

    Sym=["A","C","G","T", "N"]

    backward = [ [ 0 for i in range(len(Seq)) ] for j in range(len(T)) ]


    for s in xrange(0,len(backward),1):
        backward[s][-1] = 1
    
    for col in xrange(len(Seq)-2,-1,-1):
        for row in range(len(T)):
            for row1 in range(len(T)):
                bfirst_T=bfirst_T+ (math.log(backward[row][col+1],10)+ math.log(T[row1][row],10))
            bsecond_E= bfirst_T+math.log(E[row][Sym.index(Seq[col])],10)
            backward[row][col]=second_E
            bfirst_T=0

    return(backward)

def prob_forward_backward(F,B,Seq_size):

    P_forw_back=[ [ 0 for i in range(Seq_size) ] for j in range(len(F)) ]
    
    for i in xrange(len(F)):
        for j in xrange(Seq_size):
            P_forw_back[i][j] = F[i][j])+ B[i][j]

    for i in xrange(len(F)):
        for j in xrange(Seq_size):
            P_forw_back[i][j] = P_forw_back[i][j]/sum(P_forw_back[i])

    return(P_forw_back)

def baum_welch(states, symbols, I,  Seq ):

    import sys,math,random,numpy

    # allocate

    Sym=["A","C","G","T", "N"]

    count_o=0

    P_forw_back=[ [ 0 for i in range(len(Seq)) ] for j in range(len(T)) ]
    
    Told = [ [ random.uniform(0,1) for i in range(states) ] for j in range(states) ]

    Told_norm = [ [ Told[j][i]/sum(Told[j]) for i in range(states) ] for j in range(states) ]
    
    Eold = [ [ random.uniform(0,1) for i in range(symbols) ] for j in range(states) ]

    Eold_norm = [ [ Eold[j][i]/sum(Eold[j]) for i in range(states) ] for j in range(states) ]
    
       
    theta= [ [ [0 for i in range(len(Seq))] for i in range(states) ] for j in range(states) ]
    
    while True:
        old_T = Told_norm
        old_E = Eold_norm
        Told_norm = [ [ 1 for i in range(states) ] for j in range(states) ]
        Eold_norm = [ [ 1 for i in range(symbols) ] for j in range(states) ]

        # called forward, backward and combined prob matrices

        F = forward_algo(old_T, old_E, I, Seq)
        B = backward_algo(old_T, old_E, I, Seq)
        P = prob_forward_backward(F,B,len(Seq))

        # form 3D theta matrix
        for a1 in xrange(states):
            for a2 in xrange(states):
                for e in xrange(len(Seq)):
                    theta[a1][a2][e] = F[a1][e] + B[a2][e+1] + math.log(old_A[a1][a2],10) +  math.log(old_O[a2][Sym.index(Seq[e])],10)

        # Told_norm and Eold_norm
        for a1 in xrange(states):
            for a2 in xrange(states):
                for e in xrange(len(Seq)):
                    Told_norm[a1][a2] = sum( theta[a1][a2][e] )/ sum(P[a1])

        for a1 in xrange(states):
            for a2 in xrange(states):
                Told_norm[a1][a2] = Told_norm[a1][a2] / sum(Told_norm[a1])

        for a1 in xrange(states):
            for o in xrange(symbols):
                for e in xrange(len(Seq)):
                    if (Sym.index(Seq[e])==o)==TRUE:
                        count_o=count_o+1
                Eold_norm[a1][o] = sum(P[a1])+ math.log(count_o,10)/ sum( P[a1])
                count_o=0

        for a1 in xrange(states):
            for o in xrange(symbols):
                Eold_norm[a1][o] =  Eold_norm[a1][o]/ sum(Eold_norm[a1])

        # test condition
        if numpy.linalg.norm(old_T-Told_norm) < .00001 and numpy.linalg.norm(old_E-Eold_norm) < .00001:
            break
    # satisfied and exit
    return Told_norm, Eold_norm

