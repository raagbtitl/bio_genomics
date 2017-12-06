#!/usr/bin/env python

import re, sys
  
bed_fh = open(sys.argv[1], 'rU')
    
for line in bed_fh: 
        line = line.strip( '\n\r' ).split( '\t' )
        if re.match('#', line[0]):continue
        if len(line) != 12: # considering BED lines with 12 fields
            line = '\t'.join(line)
            sys.stdout.write('Invalid BED line found') 
            continue
        if len(line[-1].split(',')) != len(line[-2].split(',')):continue # checking the consistency b/w relative start of exon and its length
        rstart = line[-1].split(',')
        if rstart[-1] == '': rstart.pop()
        exon_len = line[-2].split(',')
        if exon_len[-1] == '': exon_len.pop()
        if len(rstart) != int(line[-3]): continue # checking the number of exons and block count are same
        if line[5] != '' and line[5] != '-':line[5] = '.' # replace the unknown starnd with '.' 
        # write feature lines to the result file 
        print line[0] + '\tbed2gff\tgene\t' + str(int(line[1]) + 1)+ '\t'  +line[2] + '\t' + line[4] + '\t' + line[5] + '\t.\t'+  'gene_id ' + '"'+ line[3]+'"' + ';gene_name ' + '"'+ line[3]+'"' 
        print line[0] + '\tbed2gff\ttranscript\t' + str(int(line[1]) + 1) + '\t'+  line[2] + '\t' + line[4] + '\t' + line[5]+  '\t.\t' + 'gene_id '+  '"'+ line[3]+'"'+  ';transcript_id '+ '"'+ line[3]+'"'+';gene_name '+  '"'+ line[3]+'"'
        st = int(line[1])
        for ex_cnt in range(int(line[-3])):
            start = st + int(rstart[ex_cnt]) + 1
	    stop = start + int(exon_len[ex_cnt]) - 1
            print line[0] + '\tbed2gff\texon\t'+  str(start) + '\t' + str(stop) + '\t'+  line[4]+  '\t'+  line[5] + '\t.\t'+  'gene_id '+  '"'+ line[3]+'"'+ ';'
bed_fh.close()


