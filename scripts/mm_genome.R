setwd("genome")
library("Biostrings")
mm9.masked <- readDNAStringSet(filepath = "Mus_musculus.NCBIM37.67.dna_rm.toplevel.fa.gz")
mm9.masked <- mm9.masked[191:211]
names(mm9.masked) <- c("chrY", "chr19", "chr18", "chr17", "chr16", "chr15", 
	"chr13", "chr12", "chr11", "chr9", "chr14", "chr10", "chr8", "chr6", "chr7", 
	"chr5", "chr4", "chr3", "chrX", "chr2", "chr1")
setwd("../RData")
save(mm9.masked, file = "mm9_masked.RData")
