source("scripts/R_functions.R")
setwd("RData")
load(file = "nature11142_s2_sacCer3.RData")
load(file = "sd01.RData")
load(file = "sc_genome_2011.RData")
load(file = "sp_genome.RData")
library("parallel")
library("Biostrings")
library("TTR")
nature11142_s2.147 <- add.nuc.seq(dyad.table = nature11142_s2.sacCer3, genome = sc.genome.2011, 
				species = "sc", nuc.size = 147)
sd01.147 <- add.nuc.seq(dyad.table = sd01, genome = sp.genome, 
				species = "sp", nuc.size = 147)
nature11142_s2.linker.147 <- get.linker.table(dyad.table = nature11142_s2.sacCer3, 
					genome = sc.genome.2011, nuc.size = 147)
sd01.linker.147 <- get.linker.table(dyad.table = sd01, genome = sp.genome, nuc.size = 147)
nature11142_s2.linker.147.seq <- get.linker.seq(linker.table = nature11142_s2.linker.147, 
							genome = sc.genome.2011)
sd01.linker.147.seq <- get.linker.seq(linker.table = sd01.linker.147, genome = sp.genome)
nature11142_s2.linker.147.seq.1_500 <- subset(nature11142_s2.linker.147.seq, 
				right - left + 1 > 0 & right - left + 1 < 501)
sd01.linker.147.seq.1_500 <- subset(sd01.linker.147.seq, 
				right - left + 1 > 0 & right - left + 1 < 501)
nature11142_s2.linker.147.seq.7_500 <- subset(nature11142_s2.linker.147.seq, 
				right - left + 1 > 6 & right - left + 1 < 501)

sd01.linker.147.seq.7_500 <- subset(sd01.linker.147.seq, 
				right - left + 1 > 6 & right - left + 1 < 501)
nature11142_s2.linker.147.length <- nature11142_s2.linker.147.seq.1_500$right - 
					nature11142_s2.linker.147.seq.1_500$left + 1
sd01.linker.147.length <- sd01.linker.147.seq.1_500$right - sd01.linker.147.seq.1_500$left + 1
nature11142_s2.linker.147.count <- get.linker.count(nature11142_s2.linker.147.length)
sd01.linker.147.count <- get.linker.count(sd01.linker.147.length)
nature11142_s2.linker.147.prob <- nature11142_s2.linker.147.count/
					sum(nature11142_s2.linker.147.count)
sd01.linker.147.prob <- sd01.linker.147.count/sum(sd01.linker.147.count)
nature11142_s2.linker.147.prob_SMA <- SMA(nature11142_s2.linker.147.prob, n = 3)
sd01.linker.147.prob_SMA <- SMA(sd01.linker.147.prob, n = 3)
nature11142_s2.linker.147.prob_SMA[1:(length(nature11142_s2.linker.147.prob_SMA) - 1)] <- 
	nature11142_s2.linker.147.prob_SMA[2:length(nature11142_s2.linker.147.prob_SMA)]
sd01.linker.147.prob_SMA[1:(length(sd01.linker.147.prob_SMA) - 1)] <- 
	sd01.linker.147.prob_SMA[2:length(sd01.linker.147.prob_SMA)]
nature11142_s2.linker.147.prob_SMA[1] <- nature11142_s2.linker.147.prob[1]
nature11142_s2.linker.147.prob_SMA[length(nature11142_s2.linker.147.prob)] <- 
		nature11142_s2.linker.147.prob[length(nature11142_s2.linker.147.prob)]
sd01.linker.147.prob_SMA[1] <- sd01.linker.147.prob[1]
sd01.linker.147.prob_SMA[length(sd01.linker.147.prob)] <- 
		sd01.linker.147.prob[length(sd01.linker.147.prob)]
nature11142_s2.linker.147.prob_SMA <- correct.prob_SMA(nature11142_s2.linker.147.prob_SMA)
sd01.linker.147.prob_SMA <- correct.prob_SMA(sd01.linker.147.prob_SMA)
uniform.linker.prob <- integer(500)
uniform.linker.prob[1:500] <- 1
uniform.linker.prob[1:500] <- uniform.linker.prob[1:500] / sum(uniform.linker.prob[1:500])

nature11142_s2.147.lDSS <- append(DNAStringSet(nature11142_s2.linker.147.seq.7_500$seq), 
			reverseComplement(DNAStringSet(nature11142_s2.linker.147.seq.7_500$seq)))
sd01.147.lDSS <- append(DNAStringSet(sd01.linker.147.seq.7_500$seq), 
			reverseComplement(DNAStringSet(sd01.linker.147.seq.7_500$seq)))
nature11142_s2.147.freqL <- oligonucleotideFrequency(nature11142_s2.147.lDSS, width = 1, as.prob = TRUE, 
							simplify.as = "collapsed", as.array = FALSE)
sd01.147.freqL <- oligonucleotideFrequency(sd01.147.lDSS, width = 1, as.prob = TRUE, 
							simplify.as = "collapsed", as.array = FALSE)
nature11142_s2.147.tranL <- oligonucleotideTransitions(nature11142_s2.147.lDSS, left = 1, as.prob = TRUE)
sd01.147.tranL <- oligonucleotideTransitions(sd01.147.lDSS, left = 1, as.prob = TRUE)
nature11142_s2.147.tranL2 <- oligonucleotideTransitions(nature11142_s2.147.lDSS, left = 2, as.prob = TRUE)
sd01.147.tranL2 <- oligonucleotideTransitions(sd01.147.lDSS, left = 2, as.prob = TRUE)
nature11142_s2.147.tranL3 <- oligonucleotideTransitions(nature11142_s2.147.lDSS, left = 3, as.prob = TRUE)
sd01.147.tranL3 <- oligonucleotideTransitions(sd01.147.lDSS, left = 3, as.prob = TRUE)
nature11142_s2.147.tranL4 <- oligonucleotideTransitions(nature11142_s2.147.lDSS, left = 4, as.prob = TRUE)
sd01.147.tranL4 <- oligonucleotideTransitions(sd01.147.lDSS, left = 4, as.prob = TRUE)
attr(nature11142_s2.147.freqL, "names") <- NULL
attr(sd01.147.freqL, "names") <- NULL
attr(nature11142_s2.147.tranL, "dimnames") <- NULL
attr(nature11142_s2.147.tranL2, "dimnames") <- NULL
attr(nature11142_s2.147.tranL3, "dimnames") <- NULL
attr(nature11142_s2.147.tranL4, "dimnames") <- NULL
attr(sd01.147.tranL, "dimnames") <- NULL
attr(sd01.147.tranL2, "dimnames") <- NULL
attr(sd01.147.tranL3, "dimnames") <- NULL
attr(sd01.147.tranL4, "dimnames") <- NULL
nature11142_s2.147.freqL <- round(nature11142_s2.147.freqL, digits = 3)
nature11142_s2.147.tranL <- round(nature11142_s2.147.tranL, digits = 3)
nature11142_s2.147.tranL2 <- round(nature11142_s2.147.tranL2, digits = 3)
nature11142_s2.147.tranL3 <- round(nature11142_s2.147.tranL3, digits = 3)
nature11142_s2.147.tranL4 <- round(nature11142_s2.147.tranL4, digits = 3)
sd01.147.freqL <- round(sd01.147.freqL, digits = 3)
sd01.147.tranL <- round(sd01.147.tranL, digits = 3)
sd01.147.tranL2 <- round(sd01.147.tranL2, digits = 3)
sd01.147.tranL3 <- round(sd01.147.tranL3, digits = 3)
sd01.147.tranL4 <- round(sd01.147.tranL4, digits = 3)
nature11142_s2.147.nDSS <- append(DNAStringSet(nature11142_s2.147$seq), 
			reverseComplement(DNAStringSet(nature11142_s2.147$seq)))
sd01.147.nDSS <- append(DNAStringSet(sd01.147$seq), 
			reverseComplement(DNAStringSet(sd01.147$seq)))
nature11142_s2.147.freqN <- oligonucleotideFrequency(
				subseq(nature11142_s2.147.nDSS, start = 1, end = 1), 
				width = 1, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE)
sd01.147.freqN <- oligonucleotideFrequency(
				subseq(sd01.147.nDSS, start = 1, end = 1), 
				width = 1, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE)
nature11142_s2.147.tranN1 <- oligonucleotideTransitions(
				subseq(nature11142_s2.147.nDSS, start = 1, end = 2), 
				left = 1, as.prob = TRUE)
sd01.147.tranN1 <- oligonucleotideTransitions(
				subseq(sd01.147.nDSS, start = 1, end = 2), 
				left = 1, as.prob = TRUE)
nature11142_s2.147.tranN2 <- oligonucleotideTransitions(
				subseq(nature11142_s2.147.nDSS, start = 1, end = 3), 
				left = 2, as.prob = TRUE)
sd01.147.tranN2 <- oligonucleotideTransitions(
				subseq(sd01.147.nDSS, start = 1, end = 3), 
				left = 2, as.prob = TRUE)
nature11142_s2.147.tranN3 <- oligonucleotideTransitions(
				subseq(nature11142_s2.147.nDSS, start = 1, end = 4), 
				left = 3, as.prob = TRUE)
sd01.147.tranN3 <- oligonucleotideTransitions(
				subseq(sd01.147.nDSS, start = 1, end = 4), 
				left = 3, as.prob = TRUE)
nature11142_s2.147.freqN2 <- getFreqN2(nature11142_s2.147.freqN, nature11142_s2.147.tranN1)
sd01.147.freqN2 <- getFreqN2(sd01.147.freqN, sd01.147.tranN1)
nature11142_s2.147.freqN3 <- getFreqN3(nature11142_s2.147.freqN2, nature11142_s2.147.tranN2)
sd01.147.freqN3 <- getFreqN3(sd01.147.freqN2, sd01.147.tranN2)
nature11142_s2.147.freqN4 <- getFreqN4(nature11142_s2.147.freqN3, nature11142_s2.147.tranN3)
sd01.147.freqN4 <- getFreqN4(sd01.147.freqN3, sd01.147.tranN3)
nature11142_s2.147.freqN4 <- round(nature11142_s2.147.freqN4, digits = 5)
sd01.147.freqN4 <- round(sd01.147.freqN4, digits = 5)
nature11142_s2.147.tranN4 <- getTranN4(nature11142_s2.147.nDSS, nuc.size = 147)
sd01.147.tranN4 <- getTranN4(sd01.147.nDSS, nuc.size = 147)
nature11142_s2.147.tranN4 <- round(nature11142_s2.147.tranN4, digits = 3)
sd01.147.tranN4 <- round(sd01.147.tranN4, digits = 3)
nature11142_s2.147.tranN4_SMA <- getTranN4_SMA(nature11142_s2.147.nDSS, nuc.size = 147)
sd01.147.tranN4_SMA <- getTranN4_SMA(sd01.147.nDSS, nuc.size = 147)
nature11142_s2.147.tranN4_SMA <- round(nature11142_s2.147.tranN4_SMA, digits = 3)
sd01.147.tranN4_SMA <- round(sd01.147.tranN4_SMA, digits = 3)

save(nature11142_s2.147, file = "nature11142_s2_147.RData")
save(nature11142_s2.linker.147, file = "nature11142_s2_linker_147.RData")
save(nature11142_s2.linker.147.seq, file = "nature11142_s2_linker_147_seq.RData")
save(sd01.147, file = "sd01_147.RData")
save(sd01.linker.147, file = "sd01_linker_147.RData")
save(sd01.linker.147.seq, file = "sd01_linker_147_seq.RData")
save(uniform.linker.prob, file = "uniform_linker_prob.RData")
save(nature11142_s2.linker.147.seq.1_500, file = "nature11142_s2_linker_147_seq_1_500.RData")
save(nature11142_s2.linker.147.seq.7_500, file = "nature11142_s2_linker_147_seq_7_500.RData")
save(nature11142_s2.linker.147.length, file = "nature11142_s2_linker_147_length.RData")
save(nature11142_s2.linker.147.count, file = "nature11142_s2_linker_147_count.RData")
save(nature11142_s2.linker.147.prob, file = "nature11142_s2_linker_147_prob.RData")
save(nature11142_s2.linker.147.prob_SMA, file = "nature11142_s2_linker_147_prob_SMA.RData")
save(sd01.linker.147.seq.1_500, file = "sd01_linker_147_seq_1_500.RData")
save(sd01.linker.147.seq.7_500, file = "sd01_linker_147_seq_7_500.RData")
save(sd01.linker.147.length, file = "sd01_linker_147_length.RData")
save(sd01.linker.147.count, file = "sd01_linker_147_count.RData")
save(sd01.linker.147.prob, file = "sd01_linker_147_prob.RData")
save(sd01.linker.147.prob_SMA, file = "sd01_linker_147_prob_SMA.RData")
save(nature11142_s2.147.freqL, file = "nature11142_s2_147_freqL.RData")
save(nature11142_s2.147.tranL, file = "nature11142_s2_147_tranL.RData")
save(nature11142_s2.147.tranL2, file = "nature11142_s2_147_tranL2.RData")
save(nature11142_s2.147.tranL3, file = "nature11142_s2_147_tranL3.RData")
save(nature11142_s2.147.tranL4, file = "nature11142_s2_147_tranL4.RData")
save(sd01.147.freqL, file = "sd01_147_freqL.RData")
save(sd01.147.tranL, file = "sd01_147_tranL.RData")
save(sd01.147.tranL2, file = "sd01_147_tranL2.RData")
save(sd01.147.tranL3, file = "sd01_147_tranL3.RData")
save(sd01.147.tranL4, file = "sd01_147_tranL4.RData")
save(nature11142_s2.147.freqN, file = "nature11142_s2_147_freqN.RData")
save(nature11142_s2.147.freqN4, file = "nature11142_s2_147_freqN4.RData")
save(nature11142_s2.147.tranN1, file = "nature11142_s2_147_tranN1.RData")
save(nature11142_s2.147.tranN2, file = "nature11142_s2_147_tranN2.RData")
save(nature11142_s2.147.tranN3, file = "nature11142_s2_147_tranN3.RData")
save(nature11142_s2.147.tranN4, file = "nature11142_s2_147_tranN4.RData")
save(nature11142_s2.147.tranN4_SMA, file = "nature11142_s2_147_tranN4_SMA.RData")
save(sd01.147.freqN, file = "sd01_147_freqN.RData")
save(sd01.147.freqN4, file = "sd01_147_freqN4.RData")
save(sd01.147.tranN1, file = "sd01_147_tranN1.RData")
save(sd01.147.tranN2, file = "sd01_147_tranN2.RData")
save(sd01.147.tranN3, file = "sd01_147_tranN3.RData")
save(sd01.147.tranN4, file = "sd01_147_tranN4.RData")
save(sd01.147.tranN4_SMA, file = "sd01_147_tranN4_SMA.RData")

rm(list = ls())
load(file = "nature11142_s2_linker_147_prob_SMA.RData")
load(file = "sd01_linker_147_prob_SMA.RData")
load(file = "nature11142_s2_147_freqL.RData")
load(file = "nature11142_s2_147_tranL.RData")
load(file = "nature11142_s2_147_tranL2.RData")
load(file = "nature11142_s2_147_tranL3.RData")
load(file = "nature11142_s2_147_tranL4.RData")
load(file = "sd01_147_freqL.RData")
load(file = "sd01_147_tranL.RData")
load(file = "sd01_147_tranL2.RData")
load(file = "sd01_147_tranL3.RData")
load(file = "sd01_147_tranL4.RData")
load(file = "nature11142_s2_147_freqN4.RData")
load(file = "nature11142_s2_147_tranN4_SMA.RData")
load(file = "sd01_147_freqN4.RData")
load(file = "sd01_147_tranN4_SMA.RData")
nature11142_s2.147.freqN4SA <- nature11142_s2.147.freqN4
sd01.147.freqN4SA <- sd01.147.freqN4
rm(nature11142_s2.147.freqN4)
rm(sd01.147.freqN4)

setwd("../nuCpos_parameters")
save(list = ls(), file = "sysdata_ScSp.rda")

