setwd("RData")
library("Biostrings")
load(file = "nature11142_s2_147.RData")
load(file = "sd01_147.RData")
load(file = "chem_mm9_NucDNA.RData")
nature11142_s2.147.nDSS <- append(DNAStringSet(nature11142_s2.147$seq), 
			reverseComplement(DNAStringSet(nature11142_s2.147$seq)))
sd01.147.nDSS <- append(DNAStringSet(sd01.147$seq), 
			reverseComplement(DNAStringSet(sd01.147$seq)))
chem.mm9.nDSS <- c(chem.mm9.NucDNA, reverseComplement(chem.mm9.NucDNA))
START <- 1
END <- START + 3
nature11142_s2.147.freqN4SA <- round(matrix(oligonucleotideFrequency(
				subseq(nature11142_s2.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
sd01.147.freqN4SA <- round(matrix(oligonucleotideFrequency(
				subseq(sd01.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
chem.mm9.freqN4SA <- round(matrix(oligonucleotideFrequency(
				subseq(chem.mm9.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
START <- 12
END <- START + 3
nature11142_s2.147.freqN4SB <- round(matrix(oligonucleotideFrequency(
				subseq(nature11142_s2.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
sd01.147.freqN4SB <- round(matrix(oligonucleotideFrequency(
				subseq(sd01.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
chem.mm9.freqN4SB <- round(matrix(oligonucleotideFrequency(
				subseq(chem.mm9.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
START <- 22
END <- START + 3
nature11142_s2.147.freqN4SC <- round(matrix(oligonucleotideFrequency(
				subseq(nature11142_s2.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
sd01.147.freqN4SC <- round(matrix(oligonucleotideFrequency(
				subseq(sd01.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
chem.mm9.freqN4SC <- round(matrix(oligonucleotideFrequency(
				subseq(chem.mm9.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
START <- 33
END <- START + 3
nature11142_s2.147.freqN4SD <- round(matrix(oligonucleotideFrequency(
				subseq(nature11142_s2.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
sd01.147.freqN4SD <- round(matrix(oligonucleotideFrequency(
				subseq(sd01.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
chem.mm9.freqN4SD <- round(matrix(oligonucleotideFrequency(
				subseq(chem.mm9.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
START <- 43
END <- START + 3
nature11142_s2.147.freqN4SE <- round(matrix(oligonucleotideFrequency(
				subseq(nature11142_s2.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
sd01.147.freqN4SE <- round(matrix(oligonucleotideFrequency(
				subseq(sd01.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
chem.mm9.freqN4SE <- round(matrix(oligonucleotideFrequency(
				subseq(chem.mm9.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
START <- 54
END <- START + 3
nature11142_s2.147.freqN4SF <- round(matrix(oligonucleotideFrequency(
				subseq(nature11142_s2.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
sd01.147.freqN4SF <- round(matrix(oligonucleotideFrequency(
				subseq(sd01.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
chem.mm9.freqN4SF <- round(matrix(oligonucleotideFrequency(
				subseq(chem.mm9.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
START <- 64
END <- START + 3
nature11142_s2.147.freqN4SG <- round(matrix(oligonucleotideFrequency(
				subseq(nature11142_s2.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
sd01.147.freqN4SG <- round(matrix(oligonucleotideFrequency(
				subseq(sd01.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
chem.mm9.freqN4SG <- round(matrix(oligonucleotideFrequency(
				subseq(chem.mm9.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
START <- 75
END <- START + 3
nature11142_s2.147.freqN4SH <- round(matrix(oligonucleotideFrequency(
				subseq(nature11142_s2.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
sd01.147.freqN4SH <- round(matrix(oligonucleotideFrequency(
				subseq(sd01.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
chem.mm9.freqN4SH <- round(matrix(oligonucleotideFrequency(
				subseq(chem.mm9.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
START <- 85
END <- START + 3
nature11142_s2.147.freqN4SI <- round(matrix(oligonucleotideFrequency(
				subseq(nature11142_s2.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
sd01.147.freqN4SI <- round(matrix(oligonucleotideFrequency(
				subseq(sd01.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
chem.mm9.freqN4SI <- round(matrix(oligonucleotideFrequency(
				subseq(chem.mm9.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
START <- 96
END <- START + 3
nature11142_s2.147.freqN4SJ <- round(matrix(oligonucleotideFrequency(
				subseq(nature11142_s2.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
sd01.147.freqN4SJ <- round(matrix(oligonucleotideFrequency(
				subseq(sd01.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
chem.mm9.freqN4SJ <- round(matrix(oligonucleotideFrequency(
				subseq(chem.mm9.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
START <- 106
END <- START + 3
nature11142_s2.147.freqN4SK <- round(matrix(oligonucleotideFrequency(
				subseq(nature11142_s2.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
sd01.147.freqN4SK <- round(matrix(oligonucleotideFrequency(
				subseq(sd01.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
chem.mm9.freqN4SK <- round(matrix(oligonucleotideFrequency(
				subseq(chem.mm9.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
START <- 117
END <- START + 3
nature11142_s2.147.freqN4SL <- round(matrix(oligonucleotideFrequency(
				subseq(nature11142_s2.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
sd01.147.freqN4SL <- round(matrix(oligonucleotideFrequency(
				subseq(sd01.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
chem.mm9.freqN4SL <- round(matrix(oligonucleotideFrequency(
				subseq(chem.mm9.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
START <- 127
END <- START + 3
nature11142_s2.147.freqN4SM <- round(matrix(oligonucleotideFrequency(
				subseq(nature11142_s2.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
sd01.147.freqN4SM <- round(matrix(oligonucleotideFrequency(
				subseq(sd01.147.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
chem.mm9.freqN4SM <- round(matrix(oligonucleotideFrequency(
				subseq(chem.mm9.nDSS, start = START, end = END), 
				width = 4, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE), 
				byrow = TRUE, ncol = 4), digits = 5)
save(nature11142_s2.147.freqN4SA, file = "nature11142_s2_147.freqN4SA.RData")
save(nature11142_s2.147.freqN4SB, file = "nature11142_s2_147.freqN4SB.RData")
save(nature11142_s2.147.freqN4SC, file = "nature11142_s2_147.freqN4SC.RData")
save(nature11142_s2.147.freqN4SD, file = "nature11142_s2_147.freqN4SD.RData")
save(nature11142_s2.147.freqN4SE, file = "nature11142_s2_147.freqN4SE.RData")
save(nature11142_s2.147.freqN4SF, file = "nature11142_s2_147.freqN4SF.RData")
save(nature11142_s2.147.freqN4SG, file = "nature11142_s2_147.freqN4SG.RData")
save(sd01.147.freqN4SA, file = "sd01_147.freqN4SA.RData")
save(sd01.147.freqN4SB, file = "sd01_147.freqN4SB.RData")
save(sd01.147.freqN4SC, file = "sd01_147.freqN4SC.RData")
save(sd01.147.freqN4SD, file = "sd01_147.freqN4SD.RData")
save(sd01.147.freqN4SE, file = "sd01_147.freqN4SE.RData")
save(sd01.147.freqN4SF, file = "sd01_147.freqN4SF.RData")
save(sd01.147.freqN4SG, file = "sd01_147.freqN4SG.RData")
save(chem.mm9.freqN4SA, file = "chem.mm9.freqN4SA.RData")
save(chem.mm9.freqN4SB, file = "chem.mm9.freqN4SB.RData")
save(chem.mm9.freqN4SC, file = "chem.mm9.freqN4SC.RData")
save(chem.mm9.freqN4SD, file = "chem.mm9.freqN4SD.RData")
save(chem.mm9.freqN4SE, file = "chem.mm9.freqN4SE.RData")
save(chem.mm9.freqN4SF, file = "chem.mm9.freqN4SF.RData")
save(chem.mm9.freqN4SG, file = "chem.mm9.freqN4SG.RData")
save(nature11142_s2.147.freqN4SH, file = "nature11142_s2_147.freqN4SH.RData")
save(nature11142_s2.147.freqN4SI, file = "nature11142_s2_147.freqN4SI.RData")
save(nature11142_s2.147.freqN4SJ, file = "nature11142_s2_147.freqN4SJ.RData")
save(nature11142_s2.147.freqN4SK, file = "nature11142_s2_147.freqN4SK.RData")
save(nature11142_s2.147.freqN4SL, file = "nature11142_s2_147.freqN4SL.RData")
save(nature11142_s2.147.freqN4SM, file = "nature11142_s2_147.freqN4SM.RData")
save(sd01.147.freqN4SH, file = "sd01_147.freqN4SH.RData")
save(sd01.147.freqN4SI, file = "sd01_147.freqN4SI.RData")
save(sd01.147.freqN4SJ, file = "sd01_147.freqN4SJ.RData")
save(sd01.147.freqN4SK, file = "sd01_147.freqN4SK.RData")
save(sd01.147.freqN4SL, file = "sd01_147.freqN4SL.RData")
save(sd01.147.freqN4SM, file = "sd01_147.freqN4SM.RData")
save(chem.mm9.freqN4SH, file = "chem.mm9.freqN4SH.RData")
save(chem.mm9.freqN4SI, file = "chem.mm9.freqN4SI.RData")
save(chem.mm9.freqN4SJ, file = "chem.mm9.freqN4SJ.RData")
save(chem.mm9.freqN4SK, file = "chem.mm9.freqN4SK.RData")
save(chem.mm9.freqN4SL, file = "chem.mm9.freqN4SL.RData")
save(chem.mm9.freqN4SM, file = "chem.mm9.freqN4SM.RData")

rm(list = ls())
load(file = "nature11142_s2_147.freqN4SA.RData")
load(file = "nature11142_s2_147.freqN4SB.RData")
load(file = "nature11142_s2_147.freqN4SC.RData")
load(file = "nature11142_s2_147.freqN4SD.RData")
load(file = "nature11142_s2_147.freqN4SE.RData")
load(file = "nature11142_s2_147.freqN4SF.RData")
load(file = "nature11142_s2_147.freqN4SG.RData")
load(file = "sd01_147.freqN4SA.RData")
load(file = "sd01_147.freqN4SB.RData")
load(file = "sd01_147.freqN4SC.RData")
load(file = "sd01_147.freqN4SD.RData")
load(file = "sd01_147.freqN4SE.RData")
load(file = "sd01_147.freqN4SF.RData")
load(file = "sd01_147.freqN4SG.RData")
load(file = "chem.mm9.freqN4SA.RData")
load(file = "chem.mm9.freqN4SB.RData")
load(file = "chem.mm9.freqN4SC.RData")
load(file = "chem.mm9.freqN4SD.RData")
load(file = "chem.mm9.freqN4SE.RData")
load(file = "chem.mm9.freqN4SF.RData")
load(file = "chem.mm9.freqN4SG.RData")
load(file = "nature11142_s2_147.freqN4SH.RData")
load(file = "nature11142_s2_147.freqN4SI.RData")
load(file = "nature11142_s2_147.freqN4SJ.RData")
load(file = "nature11142_s2_147.freqN4SK.RData")
load(file = "nature11142_s2_147.freqN4SL.RData")
load(file = "nature11142_s2_147.freqN4SM.RData")
load(file = "sd01_147.freqN4SH.RData")
load(file = "sd01_147.freqN4SI.RData")
load(file = "sd01_147.freqN4SJ.RData")
load(file = "sd01_147.freqN4SK.RData")
load(file = "sd01_147.freqN4SL.RData")
load(file = "sd01_147.freqN4SM.RData")
load(file = "chem.mm9.freqN4SH.RData")
load(file = "chem.mm9.freqN4SI.RData")
load(file = "chem.mm9.freqN4SJ.RData")
load(file = "chem.mm9.freqN4SK.RData")
load(file = "chem.mm9.freqN4SL.RData")
load(file = "chem.mm9.freqN4SM.RData")
setwd("../nuCpos_parameters")
save(list = ls(), file = "sysdata_local.rda")
