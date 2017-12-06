

def split_and_rename_fasta_file(sequence_file,pattern):
	header=''
	seq=''
	sequence='';sequence1=''
	j=0
	name = ""

	f=open(sequence_file,'r')
	while True:
		line = name or f.readline()
		if not line:
			break
		seq = []
		header = line[1:].rstrip()
		
		
		while True:
			name = f.readline()
			if not name or name.startswith(">"):
				break
			else:
				seq.append(name.rstrip())
            
		sequence = "".join(seq)
		sequence1=sequence.split('pattern')
		
		j=0
		for p in range(0,len(sequence1),1):
			j=j+1
			print '>'+ header+'-'+str(j)
			print sequence1[p]
        
      
def main():
    File=''
    pattern=''
    t=sys.argv[0:]
   
    try:
        opts, args = getopt.getopt(t[1:], "hi:p:", ["--help", "input=", "pattern="])
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
          File = a      
        elif o in ("-p", "--pattern"):     
          pattern = a
      
        else:
            sys.stderr.write("unhandled option\n")
            assert False, "unhandled option"
    split_and_rename_fasta_file(File,pattern)

def usage():
    print "\n"+ "This script is used for splitting the fasta file into multiple file by providing the desired pattern"
    print "\n"+ " -h for help, -i for input , -p for pattern"
		

if __name__ == "__main__":
    import getopt,sys,re,os
    main()        
        
    
        
    

