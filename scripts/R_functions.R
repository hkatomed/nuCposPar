add.nuc.seq.mm9.chr <- function(genome, chr = "chrY", nuc.size = 147){
	require("Biostrings")
	require("parallel")
	file.name <- paste(chr, "_unique.map_95pc.txt", sep = "")

	cat("chr: ", chr, "\n")
	LOG <- paste("chr: ", chr)
	cat("File loading, file.name: ", file.name, "\n")
	LOG[2] <- paste("File loading, file.name: ", file.name)
	dyad.table <- read.table(file.name, stringsAsFactors = FALSE)
	names(dyad.table) <- c("chr", "pos")

	cat("Before trimming, nrow(dyad.table): ", nrow(dyad.table), "\n")
	LOG[3] <- paste("Before trimming, nrow(dyad.table): ", nrow(dyad.table))

	LEN <- (nuc.size - 1)/2
	dyad.table$seq <- as.character(NA)
	if(chr == "chr1")	dyad.table <- subset(dyad.table, chr == "chr1" & pos > LEN & pos < (197195432-LEN))
	if(chr == "chr2")	dyad.table <- subset(dyad.table, chr == "chr2" & pos > LEN & pos < (181748087-LEN))
	if(chr == "chr3")	dyad.table <- subset(dyad.table, chr == "chr3" & pos > LEN & pos < (159599783-LEN))
	if(chr == "chr4")	dyad.table <- subset(dyad.table, chr == "chr4" & pos > LEN & pos < (155630120-LEN))
	if(chr == "chr5")	dyad.table <- subset(dyad.table, chr == "chr5" & pos > LEN & pos < (152537259-LEN))
	if(chr == "chr6")	dyad.table <- subset(dyad.table, chr == "chr6" & pos > LEN & pos < (149517037-LEN))
	if(chr == "chr7")	dyad.table <- subset(dyad.table, chr == "chr7" & pos > LEN & pos < (152524553-LEN))
	if(chr == "chr8")	dyad.table <- subset(dyad.table, chr == "chr8" & pos > LEN & pos < (131738871-LEN))
	if(chr == "chr9")	dyad.table <- subset(dyad.table, chr == "chr9" & pos > LEN & pos < (124076172-LEN))
	if(chr == "chr10")	dyad.table <- subset(dyad.table, chr == "chr10" & pos > LEN & pos < (129993255-LEN))
	if(chr == "chr11")	dyad.table <- subset(dyad.table, chr == "chr11" & pos > LEN & pos < (121843856-LEN))
	if(chr == "chr12")	dyad.table <- subset(dyad.table, chr == "chr12" & pos > LEN & pos < (121257530-LEN))
	if(chr == "chr13")	dyad.table <- subset(dyad.table, chr == "chr13" & pos > LEN & pos < (120284312-LEN))
	if(chr == "chr14")	dyad.table <- subset(dyad.table, chr == "chr14" & pos > LEN & pos < (125194864-LEN))
	if(chr == "chr15")	dyad.table <- subset(dyad.table, chr == "chr15" & pos > LEN & pos < (103494974-LEN))
	if(chr == "chr16")	dyad.table <- subset(dyad.table, chr == "chr16" & pos > LEN & pos < (98319150-LEN))
	if(chr == "chr17")	dyad.table <- subset(dyad.table, chr == "chr17" & pos > LEN & pos < (95272651-LEN))
	if(chr == "chr18")	dyad.table <- subset(dyad.table, chr == "chr18" & pos > LEN & pos < (90772031-LEN))
	if(chr == "chr19")	dyad.table <- subset(dyad.table, chr == "chr19" & pos > LEN & pos < (61342430-LEN))
	if(chr == "chrX")	dyad.table <- subset(dyad.table, chr == "chrX" & pos > LEN & pos < (166650296-LEN))
	if(chr == "chrY")	dyad.table <- subset(dyad.table, chr == "chrY" & pos > LEN & pos < (15902555-LEN))

	cat("After trimming, nrow(dyad.table): ", nrow(dyad.table), "\n")
	LOG[4] <- paste("After trimming, nrow(dyad.table): ", nrow(dyad.table))

		get.seq <- function(chr, pos){
			left <- pos - LEN
			right <- pos + LEN
			out.seq <- as.character(genome[[chr]][left:right])
			return(out.seq)
		}

	dyad.table$seq <- unlist(mcMap(get.seq, chr = chr, pos = dyad.table$pos, 
				mc.cores = 20))

	check.N <- function(seq){
		N.freq <- letterFrequency(DNAString(seq), letters = "N", as.prob = TRUE)
		if(N.freq == 0)	return(TRUE)
		if(N.freq > 0)	return(FALSE)
	}

	dyad.table$noN <- unlist(mcMap(f = check.N, seq = dyad.table$seq, mc.cores = 20))

	dyad.table$seq[dyad.table$noN == FALSE] <- as.character(NA)

	cat("nrow(subset(dyad.table, noN == TRUE)): ", nrow(subset(dyad.table, noN == TRUE)), "\n")
	LOG[5] <- paste("nrow(subset(dyad.table, noN == TRUE)): ", nrow(subset(dyad.table, noN == TRUE)))
	cat("nrow(dyad.table): ", nrow(dyad.table), "\n")
	LOG[6] <- paste("nrow(dyad.table): ", nrow(dyad.table))

	get.linker.length <- function(pos1, pos2, noN1, noN2, LEN){
		if(noN1 == TRUE & noN2 == TRUE){
			leftEndOfSecondNuc <- pos2 - LEN
			rightEndOfFirstNuc <- pos1 + LEN
			linker.length <- as.integer((leftEndOfSecondNuc - 1) - (rightEndOfFirstNuc + 1) + 1)
		}else{
			linker.length <- as.integer(NA)
		}
		return(linker.length)
	}

	dyad.table$linker.length <- as.integer(NA)
	dyad.table$linker.length[2:nrow(dyad.table)] <- unlist(mcMap(f = get.linker.length, 
		pos1 = dyad.table$pos[1:(nrow(dyad.table)-1)], 
		pos2 = dyad.table$pos[2:nrow(dyad.table)], 
		noN1 = dyad.table$noN[1:(nrow(dyad.table)-1)], 
		noN2 = dyad.table$noN[2:nrow(dyad.table)], 
		LEN = LEN, 
		mc.cores = 20))

	get.linker.seq <- function(pos1, pos2, noN1, noN2, LEN, chr){
		if(noN1 == TRUE & noN2 == TRUE){
			leftEndOfSecondNuc <- pos2 - LEN
			rightEndOfFirstNuc <- pos1 + LEN
			if((rightEndOfFirstNuc + 1) <= (leftEndOfSecondNuc - 1)){
				linker.seq <- as.character(genome[[chr]][(rightEndOfFirstNuc + 1):(leftEndOfSecondNuc - 1)])
				}else{
				linker.seq <- as.character(NA)
			}
		}else{
			linker.seq <- as.character(NA)
		}
		return(linker.seq)
	}

	dyad.table$linker.seq <- as.character(NA)
	dyad.table$linker.seq[2:nrow(dyad.table)] <- unlist(mcMap(f = get.linker.seq, 
		pos1 = dyad.table$pos[1:(nrow(dyad.table)-1)], 
		pos2 = dyad.table$pos[2:nrow(dyad.table)], 
		noN1 = dyad.table$noN[1:(nrow(dyad.table)-1)], 
		noN2 = dyad.table$noN[2:nrow(dyad.table)], 
		LEN = LEN, chr = chr, 
		mc.cores = 20))

	cat("nrow(subset(dyad.table, is.na(linker.length) == FALSE)): ", 
		nrow(subset(dyad.table, is.na(linker.length) == FALSE)), "\n")
	LOG[7] <- paste("nrow(subset(dyad.table, is.na(linker.length) == FALSE)): ", 
		nrow(subset(dyad.table, is.na(linker.length) == FALSE)))

	check.N.linker <- function(seq){
		if(is.na(seq) == TRUE){
			return(FALSE)
		}else{
		N.freq <- letterFrequency(DNAString(seq), letters = "N", as.prob = TRUE)
		if(N.freq == 0)	return(TRUE)
		if(N.freq > 0)	return(FALSE)
		}
	}

	dyad.table$linker.noN <- unlist(mcMap(f = check.N.linker, seq = dyad.table$linker.seq, mc.cores = 20))

	cat("nrow(subset(dyad.table, linker.noN == TRUE)): ", 
		nrow(subset(dyad.table, linker.noN == TRUE)), "\n")
	LOG[8] <- paste("nrow(subset(dyad.table, linker.noN == TRUE)): ", 
		nrow(subset(dyad.table, linker.noN == TRUE)))

	log.file.name <- paste(chr, "_add_nuc_seq_mm9_chr_log.txt", sep = "")
	write(LOG, file = log.file.name)
	out.seq <- DNAStringSet(subset(dyad.table, noN == TRUE)$seq)
	linker.seq <- DNAStringSet(subset(dyad.table, linker.noN == TRUE)$linker.seq)
	out.list <- list(out.seq, linker.seq)
	names(out.list) <- c("Nuc.DNA", "Linker.DNA")
	return(out.list)
}

add.nuc.seq <- function(dyad.table, genome, species = "sc", nuc.size = 147){
	require("Biostrings")
	require("parallel")
	LEN <- (nuc.size - 1)/2
	dyad.table$seq <- as.character(NA)
	if(species == "sc"){
		dyad.table.chr1 <- subset(dyad.table, chr == "chr1" & pos > LEN & pos < (230208-LEN))
		dyad.table.chr2 <- subset(dyad.table, chr == "chr2" & pos > LEN & pos < (813178-LEN))
		dyad.table.chr3 <- subset(dyad.table, chr == "chr3" & pos > LEN & pos < (316617-LEN))
		dyad.table.chr4 <- subset(dyad.table, chr == "chr4" & pos > LEN & pos < (1531919-LEN))
		dyad.table.chr5 <- subset(dyad.table, chr == "chr5" & pos > LEN & pos < (576869-LEN))
		dyad.table.chr6 <- subset(dyad.table, chr == "chr6" & pos > LEN & pos < (270148-LEN))
		dyad.table.chr7 <- subset(dyad.table, chr == "chr7" & pos > LEN & pos < (1090947-LEN))
		dyad.table.chr8 <- subset(dyad.table, chr == "chr8" & pos > LEN & pos < (562643-LEN))
		dyad.table.chr9 <- subset(dyad.table, chr == "chr9" & pos > LEN & pos < (439885-LEN))
		dyad.table.chr10 <- subset(dyad.table, chr == "chr10" & pos > LEN & pos < (745742-LEN))
		dyad.table.chr11 <- subset(dyad.table, chr == "chr11" & pos > LEN & pos < (666454-LEN))
		dyad.table.chr12 <- subset(dyad.table, chr == "chr12" & pos > LEN & pos < (1078175-LEN))
		dyad.table.chr13 <- subset(dyad.table, chr == "chr13" & pos > LEN & pos < (924429-LEN))
		dyad.table.chr14 <- subset(dyad.table, chr == "chr14" & pos > LEN & pos < (784333-LEN))
		dyad.table.chr15 <- subset(dyad.table, chr == "chr15" & pos > LEN & pos < (1091289-LEN))
		dyad.table.chr16 <- subset(dyad.table, chr == "chr16" & pos > LEN & pos < (948062-LEN))
		dyad.table <- rbind(dyad.table.chr1, dyad.table.chr2, dyad.table.chr3, 
			dyad.table.chr4, dyad.table.chr5, dyad.table.chr6, dyad.table.chr7, 
			dyad.table.chr8, dyad.table.chr9, dyad.table.chr10, dyad.table.chr11, 
			dyad.table.chr12, dyad.table.chr13, dyad.table.chr14, dyad.table.chr15, 
			dyad.table.chr16)
	}else if(species == "sp"){
		dyad.table.chr1 <- subset(dyad.table, chr == "chr1" & pos > LEN & pos < (5579133-LEN))
		dyad.table.chr2 <- subset(dyad.table, chr == "chr2" & pos > LEN & pos < (4539804-LEN))
		dyad.table.chr3 <- subset(dyad.table, chr == "chr3" & pos > LEN & pos < (2452883-LEN))
		dyad.table <- rbind(dyad.table.chr1, dyad.table.chr2, dyad.table.chr3)
	}
	
		get.seq <- function(chromosome, pos, strand){
			left <- pos - LEN
			right <- pos + LEN
			if(strand == "+")	out.seq <- 
				as.character(genome[[chromosome]][left:right])
			if(strand == "-")	out.seq <- 
				as.character(reverseComplement(genome[[chromosome]][left:right]))
			return(out.seq)
		}

	dyad.table$seq <- unlist(mcMap(get.seq, chromosome = dyad.table$chr, pos = dyad.table$pos, 
				strand = dyad.table$strand, mc.cores = 20))
	return(dyad.table)
}

correct.prob_SMA <- function(Pd){
	Pd.sort <- Pd[order(Pd)]
	Pd[which(Pd == 0)] <- Pd.sort[Pd.sort > 0][1]
	Pd <- Pd/sum(Pd)
}

getFreqN2 <- function(freqN, tranN1){
	for(i in 1:4){
		tranN1[i,] <- tranN1[i,]*freqN[i]
	}
	return(tranN1)
}

getFreqN3 <- function(freqN2, tranN2){
	for(i in 1:nrow(tranN2)){
		tranN2[i,] <- tranN2[i,]*t(freqN2)[i]
	}
	return(tranN2)
}

getFreqN4 <- function(freqN3, tranN3){
	for(i in 1:nrow(tranN3)){
		tranN3[i,] <- tranN3[i,]*t(freqN3)[i]
	}
	return(tranN3)
}

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

get.linker.table <- function(dyad.table, genome, nuc.size = 147){
	require("Biostrings")
	require("parallel")
	LEN <- (nuc.size - 1)/2
	chr.names <- names(genome)
	chr.lengths <- width(genome)
	linker.table <- data.frame(chr = character(0), left = integer(0), right = integer(0))
	for(i in 1:length(chr.names)){
		cat("chr.name: ", chr.names[i], "\n")
		NucLinker.index <- vector(length = chr.lengths[i])
		dyad.table.chr <- subset(dyad.table, chr == chr.names[i] & 
					pos > LEN & pos < (chr.lengths[i] - LEN))
		if(nrow(dyad.table.chr) == 0) next
		dyad.table.chr <- dyad.table.chr[order(dyad.table.chr$pos),]
		for(j in 1:nrow(dyad.table.chr)){
			left <- dyad.table.chr$pos[j] - LEN
			right <- dyad.table.chr$pos[j] + LEN
			NucLinker.index[left:right] <- TRUE
		}
		linker.start <- as.numeric(1)
		linker.end <- numeric(0)
		for(k in 2:(chr.lengths[i]-1)){
			if(NucLinker.index[k] == FALSE & NucLinker.index[k+1] == TRUE){
				linker.end <- c(linker.end, k)
			}
			if(NucLinker.index[k] == TRUE & NucLinker.index[k+1] == FALSE){
				linker.start <- c(linker.start, k + 1)
			}
		}
		cat("Before...\n")
		cat(" linker.start: ", length(linker.start), "\n")
		cat(" linker.end: ", length(linker.end), "\n")

		# 最後に TRUE then FALSE になるケース（chr11）に対応
		if(NucLinker.index[chr.lengths[i]-1] == TRUE & 
			NucLinker.index[chr.lengths[i]] == FALSE){
			linker.end <- c(linker.end, chr.lengths[i])
		}
		# 最初に FALSE then TRUE になるケース（chr12）に対応
		if(NucLinker.index[1] == FALSE & 
			NucLinker.index[2] == TRUE){
			linker.end <- c(1, linker.end)
		}
		# 最初が TRUE のケース（chr13）に対応
		if(NucLinker.index[1] == TRUE){
			linker.start <- linker.start[2:length(linker.start)]
		}
		# 最後が TRUE のケースに対応
		if(NucLinker.index[chr.lengths[i]] == TRUE){
			linker.start <- linker.start[1:(length(linker.start)-1)]
		}

		cat("After...\n")
		cat(" linker.start: ", length(linker.start), "\n")
		cat(" linker.end: ", length(linker.end), "\n")

		if(length(linker.start) == length(linker.end) + 1){
			linker.end <- c(linker.end, chr.lengths[i])
		}

		cat("Final...\n")
		cat(" linker.start: ", length(linker.start), "\n")
		cat(" linker.end: ", length(linker.end), "\n")

		# if(i == 13) return(list(dyad.table.chr, linker.start, linker.end, 
		# 			NucLinker.index))

		linker.table.chr <- data.frame(chr = chr.names[i], 
					left = linker.start, 
					right = linker.end)
		linker.table <- rbind(linker.table, linker.table.chr)
		linker.table$chr <- as.character(linker.table$chr)
	} # i-loop
	
	return(linker.table)
}

get.linker.seq <- function(linker.table, genome){
	linker.table$seq <- as.character(NA)
	# for(i in 1:nrow(linker.table)){
	# 	linker.table$seq[i] <- 
	# 		as.character(genome[[linker.table$chr[i]]][linker.table$left[i]:linker.table$right[i]])
	# }

	get.seq <- function(chromosome, left, right){
		out.seq <- as.character(genome[[chromosome]][left:right])
		return(out.seq)
		}

	linker.table$seq <- unlist(mcMap(get.seq, chromosome = linker.table$chr, 
				left = linker.table$left, right = linker.table$right, 
				mc.cores = 20))

	return(linker.table)
}

get.linker.count <- function(linker.length){
	linker.count <- numeric(500)
	for(i in 1:length(linker.length)){
		linker.count[linker.length[i]] <- linker.count[linker.length[i]] + 1
	}
	return(linker.count)
}



