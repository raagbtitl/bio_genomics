def generate_unique_kmer_OfChoice(Ifile,kmin,kmax,Ofile):

    seq='';seq1=''

    init=0;final=0

    kmer_limit_wind=''

    kmer_list=[]    
    
    
    f=open(Ifile,'r')

    data=f.readline()
    
    for line in f:
        if re.search('^>',line):
            pass  
        else:
            seq=seq+line.rstrip('\n')
    seq="".join(seq)
    print len(seq),seq
    c=0
    for k in xrange(kmin,kmax+1,1):
        print c; c=0
        for wind in xrange((len(seq)+1)-k):
            init=wind+1
            final=wind+k
            seq1=seq[wind:wind+k]
            kmer_limit_wind=','.join([seq1,str(init),str(final)])
            kmer_list.append(kmer_limit_wind)
            init=0;final=0;seq1=''
            c=c+1
    fi=open(Ofile,'w')        
    for item in kmer_list:
        fi.write(item)
        fi.write('\n')
    fi.close()
    f.close()

    
## total_kmers = len(dna) - k + 1
##    regex = re.compile(defined_pattern, re.IGNORECASE)   
##    
##    
##    for key in kmer_hash.keys():
##        kmer_values=kmer_hash[key].split(',')
##        for m in regex.finditer(kmer_values[0]):
##            size=(m.end()+1)-(m.start()+1)+1
##              coord1=int(kmer_values[1])+m.start()
##              coord2=int(kmer_values[2])+m.end()
##            print '%02d %02d %02d %s' % (size,coord1,coord2, m.group(0))

import sys,re,os,time
start_time = time.time()
generate_unique_kmer_OfChoice("Human_FMO3.fa",20,30,"patls.txt")
print time.time() - start_time,"seconds"
