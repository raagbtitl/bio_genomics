library(NOISeq)

count <- as.matrix(read.table("count2g.txt",sep = "\t", header = T, row.names = 1))

group <- data.frame(group = factor(c("A", "A", "A", "B", "B", "B")))

ndata <- readData(count, factors = group); nres <- noiseqbio(ndata, norm = "tmm", factor = "group")

table <- nres@results[[1]] ; write.table(table, file = "result.txt", col.names = T, row.names = T, sep = "\t")

##res <- noiseq(ndata, norm = "tmm", replicates = "no", factor = "group")