setwd("genome")
library("Biostrings")
sc.genome <- readDNAStringSet("S288C_reference_sequence_R64-1-1_20110203.fsa")	# bug fixed
sc.genome <- sc.genome[1:16]
names(sc.genome) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
	"chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16")
setwd("../RData")
sc.genome.2011 <- sc.genome	# sacCer3 (Apr.2011)
save(sc.genome.2011, file = "sc_genome_2011.RData")
