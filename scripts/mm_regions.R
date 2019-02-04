loaddir <- paste(getwd(), "Mm_map", sep = "/")
savedir <- paste(getwd(), "RData", sep = "/")
source("scripts/R_functions.R")
load("RData/mm9_masked.RData")
chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
"chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
"chr16", "chr17", "chr18", "chr19", "chrY", "chrX") 
for(i in 1:length(chromosomes)){
setwd(loaddir)
obj.name1 <- paste("chem.mm9.", chromosomes[i], sep = "")
assign(obj.name1, add.nuc.seq.mm9.chr(genome = mm9.masked, chr = chromosomes[i], nuc.size = 147))
obj.name2 <- paste("chem.mm9.NucDNA.", chromosomes[i], sep = "")
obj.name3 <- paste("chem.mm9.LinkerDNA.", chromosomes[i], sep = "")
assign(obj.name2, get(obj.name1)[[1]])
assign(obj.name3, get(obj.name1)[[2]])
setwd(savedir)
file.name2 <- paste("chem_mm9_NucDNA_", chromosomes[i], ".RData", sep = "")
file.name3 <- paste("chem_mm9_LinkerDNA_", chromosomes[i], ".RData", sep = "")
save(list = obj.name2, file = file.name2)
save(list = obj.name3, file = file.name3)
rm(list = c(obj.name1, obj.name2, obj.name3))
}
