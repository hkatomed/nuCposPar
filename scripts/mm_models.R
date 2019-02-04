source("scripts/R_functions.R")
setwd("RData")
library("TTR")
library("Biostrings")
NucDNA.files <- list.files(pattern = "chem_mm9_NucDNA")
for(i in 1:21) load(file = NucDNA.files[i])
NucDNA.obj.names <- ""
for(i in 1:21)	NucDNA.obj.names[i] <- strsplit(NucDNA.files[i], split = "\\.")[[1]][1]
for(i in 1:21)	NucDNA.obj.names[i] <- gsub(pattern = "_", replacement = "\\.", x = NucDNA.obj.names[i])
chem.mm9.NucDNA <- get(NucDNA.obj.names[1])
for(i in 2:21) chem.mm9.NucDNA <- c(chem.mm9.NucDNA, get(NucDNA.obj.names[i]))
rm(list = NucDNA.obj.names)
chem.mm9.NucDNA <- c(chem.mm9.NucDNA, reverseComplement(chem.mm9.NucDNA))
LinkerDNA.files <- list.files(pattern = "chem_mm9_LinkerDNA")
for(i in 1:21) load(file = LinkerDNA.files[i])
LinkerDNA.obj.names <- ""
for(i in 1:21)	LinkerDNA.obj.names[i] <- strsplit(LinkerDNA.files[i], split = "\\.")[[1]][1]
for(i in 1:21)	LinkerDNA.obj.names[i] <- gsub(pattern = "_", replacement = "\\.", x = LinkerDNA.obj.names[i])
chem.mm9.LinkerDNA <- get(LinkerDNA.obj.names[1])
for(i in 2:21) chem.mm9.LinkerDNA <- c(chem.mm9.LinkerDNA, get(LinkerDNA.obj.names[i]))
rm(list = LinkerDNA.obj.names)
chem.mm9.LinkerDNA.1_500 <- chem.mm9.LinkerDNA[width(chem.mm9.LinkerDNA) >= 1 & 
						width(chem.mm9.LinkerDNA) <= 500]
chem.mm9.LinkerDNA.7_500 <- chem.mm9.LinkerDNA[width(chem.mm9.LinkerDNA) >= 7 & 
						width(chem.mm9.LinkerDNA) <= 500]
chem.mm9.LinkerDNA.linker.count <- integer(500)
for(i in 1:length(width(chem.mm9.LinkerDNA.1_500))){
	chem.mm9.LinkerDNA.linker.count[width(chem.mm9.LinkerDNA.1_500)[i]] <- 
		chem.mm9.LinkerDNA.linker.count[width(chem.mm9.LinkerDNA.1_500)[i]] + 1
}
chem.mm9.LinkerDNA.prob <- chem.mm9.LinkerDNA.linker.count/
					sum(chem.mm9.LinkerDNA.linker.count)
chem.mm9.LinkerDNA.prob_SMA <- SMA(chem.mm9.LinkerDNA.prob, n = 3)
chem.mm9.LinkerDNA.prob_SMA[1:(length(chem.mm9.LinkerDNA.prob_SMA) - 1)] <- 
	chem.mm9.LinkerDNA.prob_SMA[2:length(chem.mm9.LinkerDNA.prob_SMA)]
chem.mm9.LinkerDNA.prob_SMA[1] <- chem.mm9.LinkerDNA.prob[1]
chem.mm9.LinkerDNA.prob_SMA[length(chem.mm9.LinkerDNA.prob)] <- 
		chem.mm9.LinkerDNA.prob[length(chem.mm9.LinkerDNA.prob)]
correct.prob_SMA <- function(Pd){
	Pd.sort <- Pd[order(Pd)]
	Pd[which(Pd == 0)] <- Pd.sort[Pd.sort > 0][1]
	Pd <- Pd/sum(Pd)
}
chem.mm9.LinkerDNA.prob_SMA <- correct.prob_SMA(chem.mm9.LinkerDNA.prob_SMA)
save(chem.mm9.NucDNA, file = "chem_mm9_NucDNA.RData")
save(chem.mm9.LinkerDNA, file = "chem_mm9_LinkerDNA.RData")
save(chem.mm9.LinkerDNA.1_500, file = "chem_mm9_LinkerDNA_1_500.RData")
save(chem.mm9.LinkerDNA.7_500, file = "chem_mm9_LinkerDNA_7_500.RData")
save(chem.mm9.LinkerDNA.linker.count, file = "chem_mm9_LinkerDNA_linker_count.RData")
save(chem.mm9.LinkerDNA.prob, file = "chem_mm9_LinkerDNA_prob.RData")
save(chem.mm9.LinkerDNA.prob_SMA, file = "chem_mm9_LinkerDNA_prob_SMA.RData")
chem.mm9.lDSS <- c(chem.mm9.LinkerDNA.7_500, 
			reverseComplement(chem.mm9.LinkerDNA.7_500))
chem.mm9.freqL <- oligonucleotideFrequency(chem.mm9.lDSS, width = 1, as.prob = TRUE, 
							simplify.as = "collapsed", as.array = FALSE)
chem.mm9.tranL <- oligonucleotideTransitions(chem.mm9.lDSS, left = 1, as.prob = TRUE)
chem.mm9.tranL2 <- oligonucleotideTransitions(chem.mm9.lDSS, left = 2, as.prob = TRUE)
chem.mm9.tranL3 <- oligonucleotideTransitions(chem.mm9.lDSS, left = 3, as.prob = TRUE)
chem.mm9.tranL4 <- oligonucleotideTransitions(chem.mm9.lDSS, left = 4, as.prob = TRUE)
attr(chem.mm9.freqL, "names") <- NULL
attr(chem.mm9.tranL, "dimnames") <- NULL
attr(chem.mm9.tranL2, "dimnames") <- NULL
attr(chem.mm9.tranL3, "dimnames") <- NULL
attr(chem.mm9.tranL4, "dimnames") <- NULL
chem.mm9.freqL <- round(chem.mm9.freqL, digits = 3)
chem.mm9.tranL <- round(chem.mm9.tranL, digits = 3)
chem.mm9.tranL2 <- round(chem.mm9.tranL2, digits = 3)
chem.mm9.tranL3 <- round(chem.mm9.tranL3, digits = 3)
chem.mm9.tranL4 <- round(chem.mm9.tranL4, digits = 3)
chem.mm9.nDSS <- chem.mm9.NucDNA
chem.mm9.freqN <- oligonucleotideFrequency(
				subseq(chem.mm9.nDSS, start = 1, end = 1), 
				width = 1, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE)
chem.mm9.tranN1 <- oligonucleotideTransitions(
				subseq(chem.mm9.nDSS, start = 1, end = 2), 
				left = 1, as.prob = TRUE)
chem.mm9.tranN2 <- oligonucleotideTransitions(
				subseq(chem.mm9.nDSS, start = 1, end = 3), 
				left = 2, as.prob = TRUE)
chem.mm9.tranN3 <- oligonucleotideTransitions(
				subseq(chem.mm9.nDSS, start = 1, end = 4), 
				left = 3, as.prob = TRUE)
getFreqN2 <- function(freqN, tranN1){
	for(i in 1:4){
		tranN1[i,] <- tranN1[i,]*freqN[i]
	}
	return(tranN1)
}
chem.mm9.freqN2 <- getFreqN2(chem.mm9.freqN, chem.mm9.tranN1)
getFreqN3 <- function(freqN2, tranN2){
	for(i in 1:nrow(tranN2)){
		tranN2[i,] <- tranN2[i,]*t(freqN2)[i]
	}
	return(tranN2)
}
chem.mm9.freqN3 <- getFreqN3(chem.mm9.freqN2, chem.mm9.tranN2)
getFreqN4 <- function(freqN3, tranN3){
	for(i in 1:nrow(tranN3)){
		tranN3[i,] <- tranN3[i,]*t(freqN3)[i]
	}
	return(tranN3)
}
chem.mm9.freqN4 <- getFreqN4(chem.mm9.freqN3, chem.mm9.tranN3)
chem.mm9.freqN4 <- round(chem.mm9.freqN4, digits = 5)
getTranN4 <- function(nDSS, nuc.size = 147){
	require(Biostrings)
	TranN4 <- matrix(numeric(0), ncol = 4)
	for(i in 1:(nuc.size-4)){
		iTranN4 <- oligonucleotideTransitions(subseq(nDSS, start = i, end = i + 4),
				 left = 4, as.prob = TRUE)
		TranN4 <- rbind(TranN4, iTranN4)
	}
	return(TranN4)
}
chem.mm9.tranN4 <- getTranN4(chem.mm9.nDSS, nuc.size = 147)
chem.mm9.tranN4 <- round(chem.mm9.tranN4, digits = 3)
getTranN4_SMA <- function(nDSS, nuc.size = 147){
	require(Biostrings)
	require(TTR)
	TranN4 <- matrix(numeric(0), nrow = 0, ncol = 4)
	TranN4.table <- matrix(numeric((nuc.size-4)*256), nrow = (nuc.size-4), ncol = 1024)
	TranN4.table.SMA <- matrix(numeric(0), nrow = (nuc.size-4), ncol = 1024)
	for(i in 1:(nuc.size-4)){
		iTranN4 <- oligonucleotideTransitions(subseq(nDSS, start = i, end = i + 4),
				 left = 4, as.prob = TRUE)
		TranN4.table[i,] <- as.numeric(t(iTranN4))
	}
	for(i in 1:1024){
		TranN4.table.SMA[,i] <- SMA(TranN4.table[,i], n = 3)
	}
	
	TranN4.table.SMA[2:(nuc.size-5),] <- TranN4.table.SMA[3:(nuc.size-4),]
	TranN4.table.SMA[1,] <- TranN4.table[1,]
	TranN4.table.SMA[(nuc.size-4),] <- TranN4.table[(nuc.size-4),]

	for(i in 1:(nuc.size-4)){
		TranN4 <- rbind(TranN4, matrix(TranN4.table.SMA[i,], byrow = TRUE, ncol = 4))
	}
	
	# SMA後の補正
	for(i in 1:nrow(TranN4)){
		TranN4[i,] <- TranN4[i,]/sum(TranN4[i,])
	}

	return(TranN4)
}
chem.mm9.tranN4_SMA <- getTranN4_SMA(chem.mm9.nDSS, nuc.size = 147)
chem.mm9.tranN4_SMA <- round(chem.mm9.tranN4_SMA, digits = 3)

save(chem.mm9.freqL, file = "chem_mm9_freqL.RData")
save(chem.mm9.tranL, file = "chem_mm9_tranL.RData")
save(chem.mm9.tranL2, file = "chem_mm9_tranL2.RData")
save(chem.mm9.tranL3, file = "chem_mm9_tranL3.RData")
save(chem.mm9.tranL4, file = "chem_mm9_tranL4.RData")
save(chem.mm9.freqN, file = "chem_mm9_freqN.RData")
save(chem.mm9.freqN4, file = "chem_mm9_freqN4.RData")
save(chem.mm9.tranN1, file = "chem_mm9_tranN1.RData")
save(chem.mm9.tranN2, file = "chem_mm9_tranN2.RData")
save(chem.mm9.tranN3, file = "chem_mm9_tranN3.RData")
save(chem.mm9.tranN4, file = "chem_mm9_tranN4.RData")
save(chem.mm9.tranN4_SMA, file = "chem_mm9_tranN4_SMA.RData")

rm(list = ls())
load(file = "chem_mm9_LinkerDNA_prob_SMA.RData")
load(file = "chem_mm9_freqL.RData")
load(file = "chem_mm9_tranL.RData")
load(file = "chem_mm9_tranL2.RData")
load(file = "chem_mm9_tranL3.RData")
load(file = "chem_mm9_tranL4.RData")
load(file = "chem_mm9_freqN4.RData")
load(file = "chem_mm9_tranN4_SMA.RData")
chem.mm9.freqN4SA <- chem.mm9.freqN4
rm(chem.mm9.freqN4)
setwd("../nuCpos_parameters")
save(list = ls(), file = "sysdata_Mm.rda")
