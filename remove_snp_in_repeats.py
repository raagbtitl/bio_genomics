

def remove_SNPs_in_repeats(vcffile,repeats_coordinate,OFile):
    snpp=[]
    rep=[]
    k=0
    a=[]


    f=open(vcffile,'rU')
    frd=f.read()
    fsp=frd.split('\n')
    fsp.pop()

    for line in fsp:
        if line.startswith("#"):
            continue
        else:
            fs=line.split()
            snpp.append(fs)

    d_snp={i[1]:i for i in snpp}


    f1=open(repeats_coordinate,'rU')
    frd1=f1.read()
    fsp1=frd1.split('\n')
    fsp1.pop()
        
    for line in fsp1:
        fs1=line.split()
        rep.append(fs1)

    for i in xrange(len(rep)-1,-1,-1):
        rep_range=xrange(int(rep[i][3]),int(rep[i][4])+1)
        a=a+rep_range

    for t in a:
        d_snp[str(t)]=0

    fi=open(OFile,'w+')    

    for k in d_snp.keys():
        if not d_snp[k]==0:
            
             fi.write('\t'.join(d_snp[k]))
             fi.write('\n')


    fi.close()
    f1.close()
    f.close()

def main():
    Vcf_file=''
    Repeat_file=''
    Output_file=''
    
    t=sys.argv[0:]
   
    try:
        opts, args = getopt.getopt(t[1:], "hv:r:o:", ["--help", "input_vcf_file=","input_repeat_file=", "output_file="])
    except getopt.GetoptError, err:
        sys.stderr.write(str(err)) 
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
        elif o in ("-v", "--input_vcf_file"):          
          Vcf_file = a      
        elif o in ("-r", "--input_repeat_file"):     
          Repeat_file = a
        elif o in ("-o","--output_file"):
            Output_file=a            
      
        else:
            sys.stderr.write("unhandled option\n")
            assert False, "unhandled option"
    remove_SNPs_in_repeats(Vcf_file,Repeat_file,Output_file)

def usage():
    print "\n"+ "This script is used for removing SNPs located in repetative regions"
    print "\n"+ " -h for help, -v for input vcf file , -r for input repeat file , -o for output file"
		

if __name__ == "__main__":
    import getopt,sys,re,os
    main()        
        
