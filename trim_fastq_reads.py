##! /mnt/python         # change the python path if needed and remove one # to make it active

def filter_reads(fqFile,fqoFile,adap,minlen):
    fq=open(fqFile,'r')
	
    
    fqls=[]
    i=0
    d={}
    for line in fq.xreadlines():
        fqls.append(line.rstrip())

    for x in xrange(i, len(fqls),4):			
        d[fqls[x+1]]=fqls[x:x+4]				# discard duplicated reads
        i=i+4

    

    for k in d.keys():

        try:
            sp1=re.split(adap,d[k][1])			# removal of adaptor
            l=len(sp1[0])
            d[k][1]=sp1[0]
            d[k][3]=d[k][3][:l]
        except:
            pass

        try:
            sp2=re.split('(#+)$',d[k][3])			# removal of low quality bases from end of reads
            L=len(sp2[0])
            d[k][3]=sp2[0]
            d[k][1]=d[k][1][:L]
        except:
            pass
        
        if len(d[k][1])< minlen:					# discard reads whose length below than minLen
            del d[k]

    print len(d)

    fqw=open(fqoFile,"w")

    for k in d.keys():
        fqw.write('\n'.join(d[k]))
        fqw.write('\n')

def main():
    fqFile=''
    fqoFile=''
    adap=''
    minlen=0
    t=sys.argv[0:]
    ##t=t.split()
    try:
        opts, args = getopt.getopt(t[1:], "hi:o:a:l:v", ["--help", "input=", "output=", "--adaptor", "--minLen"])
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
          #print "s="+a
          fqFile = a
      
        elif o in ("-o", "--output"):
     
          fqoFile = a
      
        elif o in ("-a", "--adaptor"):
	
            adap = a
        elif o in ("-l", "--minLen"):
            minlen = int(a)

        else:
            sys.stderr.write("unhandled option\n")
            assert False, "unhandled option"
    filter_reads(fqFile,fqoFile,adap,minlen)

def usage():
    print "\n"+ "This script is used for removal of a) duplicate Illumina reads, b)adaptor from reads,  c)low quality bases from each reads and d) short reads below threshold length"
    print "\n"+ " -h for help, -i for input , -o for output, -a for adaptor, -l for cuttofflength of the read"
	print "\n"+ " Caution: Should be no blank line space between read sequences"

       
      

if __name__ == "__main__":
    import getopt,sys,re,os
    main()        
        
