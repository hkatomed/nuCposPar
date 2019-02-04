setwd("genome")
library("Biostrings")
sp.genome <- readDNAStringSet(filepath = "Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz")
sp.genome <- sp.genome[1:3]
names(sp.genome) <- c("chr1", "chr2", "chr3")
setwd("../RData")
save(sp.genome, file = "sp_genome.RData")
