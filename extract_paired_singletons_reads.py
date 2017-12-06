def filter_common_pe_and_singletons(fqFile,fqFile1):
    fq=open(fqFile,'r')
    fq1=open(fqFile1,'r')

    Output_Comm=[fqFile.split(".")[0]+"_C1."+fqFile.split(".")[1],fqFile1.split(".")[0]+"_C2."+fqFile1.split(".")[1]]
    Output_Uniq=[fqFile.split(".")[0]+"_U1."+fqFile.split(".")[1],fqFile1.split(".")[0]+"_U2."+fqFile1.split(".")[1]]

    fqls=[];fqls1=[]
	
    i=0; i1=0
	
    d={};d1={};flag={}
    c1=[];c2=[];s1=[];s2=[]

    pair=0; U1=0 ; U2=0
    
    for line in fq.xreadlines():
        fqls.append(line.rstrip())
	
    for line1 in fq1.xreadlines():
        fqls1.append(line1.rstrip())

    for x in xrange(i, len(fqls),4):
        d[fqls[x][:-1]]=fqls[x:x+4]
        i=i+4



    for x1 in xrange(i1, len(fqls1),4):
        d1[fqls1[x1][:-1]]=fqls1[x1:x1+4]
        i1=i1+4
    
    flag=d1.copy()

    for k in d:
        if(k in d1):
            c1.append(d[k])
            c2.append(d1[k])
            pair=pair+1
            del flag[k] 
        else:
            s1.append(d[k])
            U2=U2+1
   

    for v in flag.keys():
        try:
            s2.append(flag[v])
            U1=U1+1        
        except:
            pass

    fqc1=open(Output_Comm[0],"w+")
    fqc2=open(Output_Comm[1],"w+")
    fqs1=open(Output_Uniq[0],"w+")
    fqs2=open(Output_Uniq[1],"w+")

    for item in c1:
        fqc1.write('\n'.join(item))
        fqc1.write('\n')
    
    for item in c2:
        fqc2.write('\n'.join(item))
        fqc2.write('\n')

    for item in s1:
        fqs1.write('\n'.join(item))
        fqs1.write('\n')
    
    for item in s2:
        fqs2.write('\n'.join(item))
        fqs2.write('\n')

    print "No. of reads in " +fqFile + " are:" +str(len(d))+ "\n"+ "No. of reads in " +fqFile1 + " are:" +str(len(d1))
    print "No. of common paired end reads are : " + str(pair) + "\n"+ "No. of singleton in " +fqFile + " are :"+   str(U1) +"\n"+ "No. of singleton in " +fqFile1 + " are:"+ str(U2)

def main():
    fqFile=''
    fqFile1=''
     
    t=sys.argv[0:]
    ##t=t.split()
    try:
        opts, args = getopt.getopt(t[1:], "hf:r:", ["--help", "fread1=", "rread2="])
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
        elif o in ("-f", "--input first paired end file"):
          #print "s="+a
          fqFile = a      
        elif o in ("-r", "--input second paired end file"):     
          fqFile1 = a
        else:
            sys.stderr.write("illegal option\n")
            assert False, "illegal option"
    filter_common_pe_and_singletons(fqFile,fqFile1)

def usage():
    print "\n"+ "This script is used for separate common paired reads from singletons"
    print "\n"+ " -h for help, -f for first paired end reads , -r for second paired end reads"
    print "\n"+ " Caution: Should be no blank line space between read sequences"      
      

if __name__ == "__main__":
    import getopt,sys,re,os
    main()        
        
