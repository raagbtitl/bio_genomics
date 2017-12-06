
def rename_read_file(sequence_file,list_header):
	header=[]; seq=''
	sequence='';p=0
	header=list_header

	f=open(fasta_file,'r')

	name = ""
	while True:
		line = name or f.readline()
		if not line:
			break
			seq = []
			while True:
				name = f.readline()
				if not name or name.startswith(">|@"):
					break
				else:
					seq.append(name.rstrip())
            p=p+1
			sequence = "".join(seq)
			print '>'+ header[p]
			print sequence

def main():
    File=''
    header_name_list=[]
    t=sys.argv[0:]
   
    try:
        opts, args = getopt.getopt(t[1:], "hi:r:", ["--help", "input=", "header_name_list="])
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
        elif o in ("-r", "--header_name_list"):     
          header_name_list = a
      
        else:
            sys.stderr.write("unhandled option\n")
            assert False, "unhandled option"
    rename_fasta_file(File,header_name_list)

def usage():
    print "\n"+ "This script is used for renaming the fasta/fastq file"
    print "\n"+ " -h for help, -i for input , -r for fasta/fastq file"
	

if __name__ == "__main__":
    import getopt,sys,re,os
    main()        
        
