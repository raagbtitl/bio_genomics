##directories <- list.files(main_dir, full.names=T);filenames <- list.files(main_dir, recursive = TRUE, full.names=T)

main_dir <- '/media/SingleWorkingDisk2/rahul/miRNA_ITLiver/miRNA_fastq/temp_nico/HLF*.csv' ; filenames <- Sys.glob(main_dir)

datalist = lapply(filenames, function(x){read.delim(file=x,sep="\t",check.names=FALSE,header=T)})

m1 <- Reduce(function(a, b) { merge(a, b,by=1,all=T) },datalist);write.table(m1,"merged.csv",sep="\t",row.names=F)
