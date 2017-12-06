def obtain_flanking_snp_seq(genome_seq,SNPPos,Len,OFile):
   
    snp=[]
    seq=''
    flank_seq=[]
    sq=[]
   
   
   
    




    f=open(genome_seq,'r')
    data=f.readline()

    for line in f:
        if re.search('^>',line):
            pass  
        else:
            seq=seq+line.rstrip('\n')

    seq="".join(seq)

    f.close()
    
    f1=open(SNPPos,'r')
    frd1=f1.read()
    fsp1=frd1.split('\n')

    for line in fsp1:
        fs1=line.split()
        snp.append(fs1)
    snp.pop()
    
    f1.close()


    for r in range(0,len(snp),1):
        n=int(snp[r][0])-1
        
        x=seq[n-Len:n]
        
        y=seq[n+1:n+(Len+1)]
        
        z=''.join(['>',snp[r][0]])
        
        z1=''.join(['[',snp[r][1],'/',snp[r][2],']'])
       
        sq=''.join([z,'\n',x,z1,y])
        
        flank_seq.append(sq)

    fi=open(OFile,'w+')

    for item in flank_seq:

        fi.write(''.join(item))

        fi.write('\n')
    
    
    fi.close()


def main():
    GenomeSeqFile=''
    
    SNP_Coord=[]
    
    Out_File=''
    
    Flank_Len=0
    
    t=sys.argv[0:]
    
    try:
        opts, args = getopt.getopt(t[1:], "hi:o:p:f:v", ["--help", "input=", "output=", "--SNP_Position", "--FlankLen"])
    except getopt.GetoptError, err:
        sys.stderr.write(str(err)) # will print error msg
        usage()
        sys.exit(2)
    output = None
    verbose = False
    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
          usage()
          sys.exit()
        elif o in ("-i", "--input"):
         
          GenomeSeqFile = a
      
        elif o in ("-o", "--output"):
     
          Out_File = a
      
        elif o in ("-p", "--SNP_Position"):
	
            SNP_Coord = a
        elif o in ("-f", "--FlankLen"):
            Flank_Len = int(a)

        else:
            sys.stderr.write("unhandled option\n")
            assert False, "unhandled option"
    obtain_flanking_snp_seq(GenomeSeqFile,SNP_Coord,Flank_Len,Out_File)

def usage():
    print "\n"+ "This script is used for extracting SNP flanking sequence from particular chromosome of a genomic sequence"
    print "\n"+ " -h for help, -i for input , -o for output, -p for SNP Positions, -f for length of flanking sequence"
    print "\n"+ " Format of SNP Position file : SNP_Coordinate Reference_Allele Alternate_Allele e.g. 12186390 A C"

       
      

if __name__ == "__main__":
    import getopt,sys,re,os
    main()  
 

