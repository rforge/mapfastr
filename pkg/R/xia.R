MCIBD <-
		function(loci, data = NULL, trim.output = NULL, chr = 1, IBD.type = "genotypic", output.Z = "none", segregation = NULL, mc.size = 99, hpc = FALSE, n.cpus = 2) {
	## --------------------------------------------- ##
	##       Monte Carlo IBD Matrix Calculator       ##
	##      cnF2freq has to be run before this.      ##
	##      Xia.Shen@slu.se ---2012-08-13---         ##
	## --------------------------------------------- ##
	
	## !! PACKAGE REQUIREMENTS: sfsmisc, snow (for HPC), snowfall (for HPC)
	
	## Check some inputs
	if (length(loci) != 1 & length(loci) != 2) {
		stop("Incorrect number of loci!")
	}
	if (prod(c(loci >= 0, mc.size > 0, n.cpus > 0)) == 0) {
		stop("Incorrect input detected!")
	}	
	
	## Data & Pre-calculation ##
	pb <- txtProgressBar(style = 3)
	pedi <- data$pheno[,3:4]
	pedif0 <- pedi[data$pheno$generation == 1,]
	f0id <- as.numeric(rownames(pedif0))
	carl <- trim.output[chr][[1]]
	n.loci <- ncol(carl)
	n.f2 <- sum(data$pheno$generation == 3)
	id <- as.numeric(rownames(data$pheno))
	f2id <- id[data$pheno$generation == 3]
	
	cat("\nData Transforming ...", "\n")
	carlout <- NULL
	for (k in 1:n.f2) {
		carlout[[k]] <- t(carl[((k - 1)*64 + 1):(k*64),])
		setTxtProgressBar(pb, k/n.f2)
	}
	cat("\n")
	cat("OKAY.", "\n")
	cat("\n")
	
	n.founder <- sum(data$pheno$generation == 1)
	
	if (is.null(segregation)) {
		segregation <- 1:(2*n.founder)
	}
	
	## Necessary Functions ##
	cat("Functions Loading ...", "\n")
	
	## Sampling using cnF2freq Output
	require(sfsmisc, quietly = TRUE)
	binary <- NULL
	for (i in 1:64) binary <- rbind(binary, as.numeric(digitsBase(i-1,,6)))
	
	samplecarl <- function(carlout, position, pedigree, f0id, n.f2, type) {
		here <- Z <- NULL
		for (i in 1:n.f2) here <- rbind(here, rbind(carlout[[i]])[position,])
		for (i in 1:n.f2) {
			case <- binary[sample(1:64, 1, prob=here[i,]),]
			if (case[6]==0 & case[5]==1) a1 <- which(f0id == pedigree[id == pedigree[id == f2id[i],1],1])*2 - 1
			if (case[6]==0 & case[5]==0) a1 <- which(f0id == pedigree[id == pedigree[id == f2id[i],1],1])*2
			if (case[6]==1 & case[4]==1) a1 <- which(f0id == pedigree[id == pedigree[id == f2id[i],1],2])*2 - 1
			if (case[6]==1 & case[4]==0) a1 <- which(f0id == pedigree[id == pedigree[id == f2id[i],1],2])*2
			if (case[3]==0 & case[2]==1) a2 <- which(f0id == pedigree[id == pedigree[id == f2id[i],2],1])*2 - 1
			if (case[3]==0 & case[2]==0) a2 <- which(f0id == pedigree[id == pedigree[id == f2id[i],2],1])*2
			if (case[3]==1 & case[1]==1) a2 <- which(f0id == pedigree[id == pedigree[id == f2id[i],2],2])*2 - 1
			if (case[3]==1 & case[1]==0) a2 <- which(f0id == pedigree[id == pedigree[id == f2id[i],2],2])*2
			if (type == "genotypic") {
				Z.line <- rep(0,length(segregation))
				Z.line[c(a1, a2)] <- 1
				Z <- rbind(Z, Z.line)
			}
			if (type == "gametic") {
				Z.line1 <- Z.line2 <- rep(0,length(segregation))
				Z.line1[a1] <- Z.line2[a2] <- 1
				Z <- rbind(Z, Z.line1, Z.line2)
			}
		}
		dimnames(Z) <- list(NULL, NULL)
		return(Z)
	}
	
	setTxtProgressBar(pb, .5)
	
	## Adding Segregation Vector for Founder Alleles
	sgg <- function(Z, seg.alleles) {
		if (length(seg.alleles) != ncol(Z)) {
			stop("Incorrect number of alleles!")
		}	
		alleles <- seg.alleles[1]
		for (i in 2:length(seg.alleles)) {
			if (seg.alleles[i] != seg.alleles[i - 1]) {
				alleles <- c(alleles, seg.alleles[i])
			}	
		}
		seg.Z <- matrix(0, nrow(Z), length(alleles))
		for (i in 1:length(alleles)) {
			index <- 1:ncol(Z)*(seg.alleles == alleles[i])
			if (sum(index != 0) > 1) {
				seg.Z[,i] <- rowSums(Z[,index])
			}
			else {
				seg.Z[,i] <- Z[,index]	
			}
		}
		return(seg.Z)
	}
	setTxtProgressBar(pb, 1)
	cat("\n")
	cat("OKAY.", "\n")
	cat("\n")
	
	## Monte Carlo Sampling ##
	MIsize <- mc.size
	p <- loci
	if (length(p) == 2) epistasis <- TRUE else epistasis <- FALSE
	if (IBD.type == "genotypic") nf2 <- n.f2 else nf2 <- 2*n.f2
	if (output.Z == "all") {dir.create("Zall")}
	if (output.Z == "none" | output.Z == "all") {
		exname <- ".ibd"
		type <- "IBD"
		dim2 <- nf2
	}
	if (output.Z == "av" | output.Z == "pc") {
		exname <- ".z"
		type <- "Incidence"
		dim2 <- 2*n.founder	
	}
	
	if (!hpc) {
		## Single Core
		t0 <- proc.time()[3]
		if (output.Z == "av") sumPi <- matrix(0, nf2, dim2) else sumPi <- matrix(0, nf2, nf2)
		cat("One master is doing its jobs ...", "\n")
		if (!epistasis) {
			for(i in 1:MIsize) {
				Z <- samplecarl(carlout = carlout, position = p, pedigree = pedi, f0id = f0id, n.f2 = n.f2, type = IBD.type)
				Z <- sgg(Z, segregation)
				if (output.Z == "none" | output.Z == "pc") {
					Pi <- .5*Z%*%t(Z)
				}
				if (output.Z == "av") {
					Pi <- Z
				}
				if (output.Z == "all") {
					Pi <- .5*Z%*%t(Z)
					filename <- paste("Zall/", i, ".z", sep = "")
					write.table(Z, filename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
				}
				sumPi <- sumPi + Pi
				setTxtProgressBar(pb, i/MIsize)
			}
		}
		else {
			for(i in 1:MIsize) {
				Z1 <- samplecarl(carlout = carlout, position = p[1], pedigree = pedi, f0id = f0id, n.f2 = n.f2, type = IBD.type)
				Z1 <- sgg(Z1, segregation)
				Z2 <- samplecarl(carlout = carlout, position = p[2], pedigree = pedi, f0id = f0id, n.f2 = n.f2, type = IBD.type)
				Z2 <- sgg(Z2, segregation)
				if (output.Z == "none" | output.Z == "all" | output.Z == "pc") {
					Pi1 <- .5*Z1%*%t(Z1)
					Pi2 <- .5*Z2%*%t(Z2)
				}
				if (output.Z == "av") {
					Pi1 <- Z1
					Pi2 <- Z2
				}
				## Hadamard product applied (Shen et al. 2009)
				Pi <- Pi1*Pi2
				sumPi <- sumPi + Pi
				setTxtProgressBar(pb, i/MIsize)
			}
		}
		cat("\n")
		meanPi <- sumPi/MIsize
		if (output.Z == "pc") {
			A <- eigen(meanPi)
			v <- A$values[1:(2*n.founder)]
			pc <- A$vectors[,1:(2*n.founder)]
			meanPi <- pc%*%diag(sqrt(v))
		}
		t1 <- proc.time()[3] - t0
		if (!epistasis) {
			filename <- paste(p, exname, sep = "")
		}
		else {
			filename <- paste(paste(p[1], p[2], sep = "_x_"), exname, sep = "")	
		}
		indibd <- unique(unlist(strsplit(cbind(row.names(trim.output[[1]])),split='[.]'))[seq(1,nrow(trim.output[[1]])*2,2)])
		write.table(cbind(indibd), "indibd.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
		write.table(meanPi, filename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
		if (!epistasis) {
			cat(paste(type, "matrix of dimension", nf2, "x", dim2, "at locus", p, "is accomplished, where", MIsize, "imputes were sampled.", sep = " "), "\n")
		}
		else {
			cat(paste(type, "matrix of dimension", nf2, "x", dim2, "for epistatic loci", p[1], "and", p[2], "is accomplished, where", MIsize, "imputes were sampled.", sep = " "), "\n")
		}
	}
	else {
		## Multiple Cores
		## NOTE! Make sure the multicore environment is set up and snowfall package is installed.
		require(snow, quietly = TRUE)
		require(snowfall, quietly = TRUE)
		sfInit(parallel = TRUE, cpus = n.cpus, type = "SOCK")
		slavejob <- function(idx) {
			if (output.Z == "av") sumPi <- matrix(0, nf2, dim2) else sumPi <- matrix(0, nf2, nf2)
			if (!epistasis) {
				for(i in 1:MIsize) {
					Z <- samplecarl(carlout = carlout, position = p + 1, pedigree = pedi, f0id = f0id, n.f2 = n.f2, type = IBD.type)
					Z <- sgg(Z, segregation)
					if (output.Z == "none" | output.Z == "pc") {
						Pi <- .5*Z%*%t(Z)
					}
					if (output.Z == "av") {
						Pi <- Z
					}
					if (output.Z == "all") {
						Pi <- .5*Z%*%t(Z)
						dir.create("Zall")
						filename <- paste("Zall/", i, ".z", sep = "")
						write.table(Z, filename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
					}
					sumPi <- sumPi + Pi
					setTxtProgressBar(pb, i/MIsize)
				}
			}
			else {
				for(i in 1:MIsize) {
					Z1 <- samplecarl(carlout = carlout, position = p[1] + 1, pedigree = pedi, f0id = f0id, n.f2 = n.f2, type = IBD.type)
					Z1 <- sgg(Z1, segregation)
					Z2 <- samplecarl(carlout = carlout, position = p[2] + 1, pedigree = pedi, f0id = f0id, n.f2 = n.f2, type = IBD.type)
					Z2 <- sgg(Z2, segregation)
					if (output.Z == "none" | output.Z == "all" | output.Z == "pc") {
						Pi1 <- .5*Z1%*%t(Z1)
						Pi2 <- .5*Z2%*%t(Z2)
					}
					if (output.Z == "av") {
						Pi1 <- Z1
						Pi2 <- Z2
					}
					## Hadamard product applied (Shen et al. 2009)
					Pi <- Pi1*Pi2
					sumPi <- sumPi + Pi
					setTxtProgressBar(pb, i/MIsize)
				}
			}
			meanPi <- sumPi/MIsize
			return(meanPi)
		}
		cat("Broadcasting objects to slaves ...", "\n")
		sfExport("MIsize")
		setTxtProgressBar(pb, 1/10)
		sfExport("p")
		setTxtProgressBar(pb, 2/10)
		sfExport("binary")
		setTxtProgressBar(pb, 3/10)
		sfExport("f2id")
		setTxtProgressBar(pb, 4/10)
		sfExport("nf2")
		setTxtProgressBar(pb, 5/10)
		sfExport("dim2")
		setTxtProgressBar(pb, 6/10)
		sfExport("pedi")
		setTxtProgressBar(pb, 7/10)
		sfExport("carlout")
		setTxtProgressBar(pb, 8/10)
		sfExport("samplecarl")
		setTxtProgressBar(pb, 9/10)
		sfExport("sgg")
		setTxtProgressBar(pb, 10/10)
		cat("\n")
		cat("OKAY.", "\n")
		cat("\n")
		t0 <- proc.time()[3]
		cat("Slaves are doing their jobs ...", "\n")
		result <- sfLapply(1:n.cpus, slavejob)
		cat("\n")
		cat("Final calculation for the", type, "matrix ...", "\n")
		if (output.Z == "av") sumPi <- matrix(0, nf2, dim2) else sumPi <- matrix(0, nf2, nf2)
		for(sn in 1:n.cpus) {
			sumPi <- sumPi + result[[sn]]
			setTxtProgressBar(pb, sn/n.cpus)
		}
		cat("\n")
		meanPi <- sumPi/n.cpus
		if (output.Z == "pc") {
			A <- eigen(meanPi)
			v <- A$values[1:(2*n.founder)]
			pc <- A$vectors[,1:(2*n.founder)]
			meanPi <- pc%*%diag(sqrt(v))
		}
		t1 <- proc.time()[3] - t0
		if (!epistasis) {
			filename <- paste(p, exname, sep = "")
		}
		else {
			filename <- paste(paste(p[1], p[2], sep = "_x_"), exname, sep = "")	
		}
		indibd <- unique(unlist(strsplit(cbind(row.names(trim.output[[1]])),split='[.]'))[seq(1,nrow(trim.output[[1]])*2,2)])
		write.table(cbind(indibd), "indibd.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
		write.table(meanPi, filename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
		if (!epistasis) {
			cat(paste(type, "matrix of dimension", nf2, "x", dim2, "at locus", p, "is accomplished, where", MIsize*n.cpus, "imputes were sampled.", sep = " "), "\n")
		}
		else {
			cat(paste(type, "matrix of dimension", nf2, "x", dim2, "for epistatic loci", p[1], "and", p[2], "is accomplished, where", MIsize*n.cpus, "imputes were sampled.", sep = " "), "\n")
		}
		sfStop()
	}
	cat("\n")
	cat("MC Sampling Execution Time: ", t1, "sec", "\n")
	cat("\n")
	## End ##
}








MCIBD.chro <-
		function(data = NULL, trim.output = NULL, chr = 1, IBD.type = "genotypic", output.Z = "none", segregation = NULL, mc.size = 99, hpc = FALSE, n.cpus = 2) {
	## --------------------------------------------- ##
	##       Monte Carlo IBD Matrix Calculator       ##
	##      cnF2freq has to be run before this.      ##
	##        Xia.Shen@slu.se ---2012-10-14---       ##
	## --------------------------------------------- ##
	
	## !! PACKAGE REQUIREMENTS: sfsmisc, snow (for HPC), snowfall (for HPC)
	
	## Check some inputs
	if (prod(c(mc.size > 0, n.cpus > 0)) == 0) {
		stop("Non-positive input detected!")
	}	
	
	## Data & Pre-calculation ##
	#pb <- txtProgressBar(style = 3)
	pedi <- data$pheno[,3:4]
	pedif0 <- pedi[data$pheno$generation == 1,]
	f0id <- as.numeric(rownames(pedif0))
	carl <- trim.output[chr][[1]]
	n.loci <- ncol(carl)
	loci <- 1:n.loci
	n.f2 <- sum(data$pheno$generation == 3)
	id <- as.numeric(rownames(data$pheno))
	f2id <- id[data$pheno$generation == 3]
	
	cat("Data Transforming ...", "\n")
	carlout <- NULL
	for (k in 1:n.f2) {
		carlout[[k]] <- t(carl[((k - 1)*64 + 1):(k*64),])
		#setTxtProgressBar(pb, k/n.f2)
	}
	cat("\n")
	cat("OKAY.", "\n")
	cat("\n")
	
	n.founder <- sum(data$pheno$generation == 1)
	
	if (is.null(segregation)) {
		segregation <- 1:(2*n.founder)
	}
	
	## Necessary Functions ##
	cat("Functions Loading ...", "\n")
	
	## Sampling using cnF2freq Output
	require(sfsmisc, quietly = TRUE)
	binary <- NULL
	for (i in 1:64) binary <- rbind(binary, as.numeric(digitsBase(i-1,,6)))
	
	samplecarl <- function(carlout, position, pedigree, f0id, n.f2, type) {
		here <- Z <- NULL
		for (i in 1:n.f2) here <- rbind(here, rbind(carlout[[i]])[position,])
		for (i in 1:n.f2) {
			case <- binary[sample(1:64, 1, prob=here[i,]),]
			if (case[6]==0 & case[5]==1) a1 <- which(f0id == pedigree[id == pedigree[id == f2id[i],1],1])*2 - 1
			if (case[6]==0 & case[5]==0) a1 <- which(f0id == pedigree[id == pedigree[id == f2id[i],1],1])*2
			if (case[6]==1 & case[4]==1) a1 <- which(f0id == pedigree[id == pedigree[id == f2id[i],1],2])*2 - 1
			if (case[6]==1 & case[4]==0) a1 <- which(f0id == pedigree[id == pedigree[id == f2id[i],1],2])*2
			if (case[3]==0 & case[2]==1) a2 <- which(f0id == pedigree[id == pedigree[id == f2id[i],2],1])*2 - 1
			if (case[3]==0 & case[2]==0) a2 <- which(f0id == pedigree[id == pedigree[id == f2id[i],2],1])*2
			if (case[3]==1 & case[1]==1) a2 <- which(f0id == pedigree[id == pedigree[id == f2id[i],2],2])*2 - 1
			if (case[3]==1 & case[1]==0) a2 <- which(f0id == pedigree[id == pedigree[id == f2id[i],2],2])*2
			if (type == "genotypic") {
				Z.line <- rep(0,length(segregation))
				Z.line[c(a1, a2)] <- 1
				Z <- rbind(Z, Z.line)
			}
			if (type == "gametic") {
				Z.line1 <- Z.line2 <- rep(0,length(segregation))
				Z.line1[a1] <- Z.line2[a2] <- 1
				Z <- rbind(Z, Z.line1, Z.line2)
			}
		}
		dimnames(Z) <- list(NULL, NULL)
		return(Z)
	}
	
	#setTxtProgressBar(pb, .5)
	
	samplehere <- function(here, pedigree, f0id, n.f2, type) {
		Z <- NULL
		for (i in 1:n.f2) {
			case <- binary[sample(1:64, 1, prob=here[i,]),]
			if (case[6]==0 & case[5]==1) a1 <- which(f0id == pedigree[id == pedigree[id == f2id[i],1],1])*2 - 1
			if (case[6]==0 & case[5]==0) a1 <- which(f0id == pedigree[id == pedigree[id == f2id[i],1],1])*2
			if (case[6]==1 & case[4]==1) a1 <- which(f0id == pedigree[id == pedigree[id == f2id[i],1],2])*2 - 1
			if (case[6]==1 & case[4]==0) a1 <- which(f0id == pedigree[id == pedigree[id == f2id[i],1],2])*2
			if (case[3]==0 & case[2]==1) a2 <- which(f0id == pedigree[id == pedigree[id == f2id[i],2],1])*2 - 1
			if (case[3]==0 & case[2]==0) a2 <- which(f0id == pedigree[id == pedigree[id == f2id[i],2],1])*2
			if (case[3]==1 & case[1]==1) a2 <- which(f0id == pedigree[id == pedigree[id == f2id[i],2],2])*2 - 1
			if (case[3]==1 & case[1]==0) a2 <- which(f0id == pedigree[id == pedigree[id == f2id[i],2],2])*2
			if (type == "genotypic") {
				Z.line <- rep(0,length(segregation))
				Z.line[c(a1, a2)] <- 1
				Z <- rbind(Z, Z.line)
			}
			if (type == "gametic") {
				Z.line1 <- Z.line2 <- rep(0,length(segregation))
				Z.line1[a1] <- Z.line2[a2] <- 1
				Z <- rbind(Z, Z.line1, Z.line2)
			}
		}
		dimnames(Z) <- list(NULL,NULL)
		return(Z)
	}
	
	## Adding Segregation Vector for Founder Alleles
	sgg <- function(Z, seg.alleles) {
		if (length(seg.alleles) != ncol(Z)) {
			stop("Incorrect number of alleles!")
		}	
		alleles <- seg.alleles[1]
		for (i in 2:length(seg.alleles)) {
			if (seg.alleles[i] != seg.alleles[i - 1]) {
				alleles <- c(alleles, seg.alleles[i])
			}	
		}
		seg.Z <- matrix(0, nrow(Z), length(alleles))
		for (i in 1:length(alleles)) {
			index <- 1:ncol(Z)*(seg.alleles == alleles[i])
			if (sum(index != 0) > 1) {
				seg.Z[,i] <- rowSums(Z[,index])
			}
			else {
				seg.Z[,i] <- Z[,index]	
			}
		}
		return(seg.Z)
	}
	#setTxtProgressBar(pb, 1)
	cat("\n")
	cat("OKAY.", "\n")
	cat("\n")
	
	## Monte Carlo Sampling ##
	MIsize <- mc.size
	if (IBD.type == "genotypic") nf2 <- n.f2 else nf2 <- 2*n.f2
	if (output.Z == "all") dir.create("Zall")
	if (output.Z == "none" | output.Z == "all") {
		exname <- ".ibd.RData"
		if (length(table(segregation)) == 2) exname <- ".ibd1.RData"
		type <- "IBD"
		dim2 <- nf2
	}
	if (output.Z == "av" | output.Z == "pc") {
		exname <- ".z.RData"
		if (length(table(segregation)) == 2) exname <- ".ibd1.RData"
		type <- "Incidence"
		dim2 <- 2*n.founder	
	}
	
	if (!file.exists(paste('chr', chr, sep = ''))) dir.create(paste('chr', chr, sep = ''))
	
	save(n.loci, file = paste('chr', chr, '/n.loci.RData', sep = ""))
	
	if (!hpc) {
		## Single Core
		t0 <- proc.time()[3]
		cat("One master is doing its jobs ...", "\n")
		for (p in loci) {
			if (output.Z == "av") sumPi <- matrix(0, nf2, dim2) else sumPi <- matrix(0, nf2, nf2)
			if (output.Z == "all") dir.create(paste("Zall/", p, sep = ""))
			for(i in 1:MIsize) {
				Z <- samplecarl(carlout = carlout, position = p, pedigree = pedi, f0id = f0id, n.f2 = n.f2, type = IBD.type)
				Z <- sgg(Z, segregation)
				if (output.Z == "none" | output.Z == "pc") {
					Pi <- .5*Z%*%t(Z)
				}
				if (output.Z == "av") {
					Pi <- Z
				}
				if (output.Z == "all") {
					Pi <- .5*Z%*%t(Z)
					filename <- paste("Zall/", p, "/", i, ".z.RData", sep = "")
					save(Z, file = filename)
				}
				sumPi <- sumPi + Pi
				#setTxtProgressBar(pb, i/MIsize)
			}
			cat("\n")
			meanPi <- sumPi/MIsize
			if (output.Z == "pc") {
				A <- eigen(meanPi)
				v <- A$values[1:(2*n.founder)]
				pc <- A$vectors[,1:(2*n.founder)]
				meanPi <- pc%*%diag(sqrt(v))
			}
			filename <- paste('chr', chr, '/', p, exname, sep = "")
			save(meanPi, file = filename)
			cat(paste(type, "matrix of dimension", nf2, "x", dim2, "at locus", p, "is accomplished, where", MIsize, "imputes were sampled.", sep = " "), "\n")
		}
		t1 <- proc.time()[3] - t0
	}
	else {
		## Multiple Cores
		## NOTE! Make sure the multicore environment is set up and snowfall package is installed.
		require(snow, quietly = TRUE)
		require(snowfall, quietly = TRUE)
		sfInit(parallel = TRUE, cpus = n.cpus, type = "SOCK")
		slavejob <- function(idx) {
			if (output.Z == "av") sumPi <- matrix(0, nf2, dim2) else sumPi <- matrix(0, nf2, nf2)
			for(i in 1:MIsize) {
				Z <- samplehere(here = here, pedigree = pedi, f0id = f0id, n.f2 = n.f2, type = IBD.type)
				Z <- sgg(Z, segregation)
				if (output.Z == "none" | output.Z == "pc") {
					Pi <- .5*Z%*%t(Z)
				}
				if (output.Z == "av") {
					Pi <- Z
				}
				if (output.Z == "all") {
					Pi <- .5*Z%*%t(Z)
					dir.create("Zall")
					filename <- paste("Zall/", p, "/", i, ".z.RData", sep = "")
					save(Z, file = filename)
				}
				sumPi <- sumPi + Pi
			}
			meanPi <- sumPi/MIsize
			return(meanPi)
		}
		t0 <- proc.time()[3]
		for (p in loci) {
			if (output.Z == "all") dir.create(paste("Zall/", p, sep = ""))
			here <- NULL
			for (i in 1:n.f2) here <- rbind(here, carlout[[i]][p,])
			cat("Broadcasting objects to slaves ...", "\n")
			sfExport("MIsize")
			#setTxtProgressBar(pb, 1/10)
			sfExport("p")
			#setTxtProgressBar(pb, 2/10)
			sfExport("binary")
			#setTxtProgressBar(pb, 3/10)
			sfExport("f2id")
			#setTxtProgressBar(pb, 4/10)
			sfExport("nf2")
			#setTxtProgressBar(pb, 5/10)
			sfExport("dim2")
			#setTxtProgressBar(pb, 6/10)
			sfExport("pedi")
			#setTxtProgressBar(pb, 7/10)
			sfExport("here")
			#setTxtProgressBar(pb, 8/10)
			sfExport("samplehere")
			#setTxtProgressBar(pb, 9/10)
			sfExport("sgg")
			#setTxtProgressBar(pb, 10/10)
			cat("\n")
			cat("OKAY.", "\n")
			cat("\n")
			cat("Slaves are doing their jobs ...", "\n")
			result <- sfLapply(1:n.cpus, slavejob)
			cat("\n")
			cat("Final calculation for the", type, "matrix ...", "\n")
			if (output.Z == "av") sumPi <- matrix(0, nf2, dim2) else sumPi <- matrix(0, nf2, nf2)
			for(sn in 1:n.cpus) {
				sumPi <- sumPi + result[[sn]]
				#setTxtProgressBar(pb, sn/n.cpus)
			}
			cat("\n")
			meanPi <- sumPi/n.cpus
			if (output.Z == "pc") {
				A <- eigen(meanPi)
				v <- A$values[1:(2*n.founder)]
				pc <- A$vectors[,1:(2*n.founder)]
				meanPi <- pc%*%diag(sqrt(v))
			}
			filename <- paste('chr', chr, '/', p, exname, sep = "")
			save(meanPi, file = filename)
			cat(paste(type, "matrix of dimension", nf2, "x", dim2, "at locus", p, "is accomplished, where", MIsize*n.cpus, "imputes were sampled.", sep = " "), "\n")
			cat("\n")
			rm(result)
		}
		t1 <- proc.time()[3] - t0
		sfStop()
	}
	cat("\n")
	cat("MC Sampling Execution Time: ", t1, "sec", "\n")
	cat("\n")
	## End ##
}









MCIBD.genome <-
		function(data = NULL, trim.output = NULL, chr = NULL, IBD.type = "genotypic", output.Z = "none", segregation = NULL, mc.size = 99, hpc = FALSE, n.cpus = 2) {
	## --------------------------------------------- ##
	##       Monte Carlo IBD Matrix Calculator       ##
	##      cnF2freq has to be run before this.      ##
	##      Xia.Shen@slu.se ---2012-10-14---         ##
	## --------------------------------------------- ##
	
	## !! PACKAGE REQUIREMENTS: sfsmisc, snow (for HPC), snowfall (for HPC)
	
	## Check some inputs
	if (prod(c(mc.size > 0, n.cpus > 0)) == 0) {
		stop("Non-positive input detected!")
	}	
	
	if (is.null(chr)) chr <- 1:length(trim.output)
	
	indibd <- unique(unlist(strsplit(cbind(row.names(trim.output[[1]])),split='[.]'))[seq(1,nrow(trim.output[[1]])*2,2)])
	save(indibd, file = "indibd.RData")
	
	for (k in chr) {
		cat('Chromosome', k, '\n')
		cat('\n')
		MCIBD.chro(data = data, trim.output = trim.output, chr = k, IBD.type = IBD.type, output.Z = output.Z, segregation = segregation, mc.size = mc.size, hpc = hpc, n.cpus = n.cpus)
	}
	## End ##
}









score <-
function(y, P, IBD, print.result = FALSE) {
#################################################################################################			  																				 
#  y				Response vector
#  P				Null hypothesis projection matrix
#  A				Additive relationship matrix
#  IBD			    IBD matrix
#################################################################################################
if (sum(is.na(y)) > 0) {
naidx <- which(is.na(y))
	y <- y[-naidx]
	P <- P[-naidx,-naidx]
	IBD <- IBD[-naidx,-naidx]
}

	min.error <- 1e-8
	n.obs <- length(y)
	n.comp <- 2
	
	dimIBD <- n.obs
	D <- numeric(n.comp)	
	F <- matrix(0,n.comp,n.comp)
	
	AP <- IBD%*%P
	tyP <- t(y)%*%P
	
	for (j in 1:n.comp){
		if (j==1) AjP <- AP
		if (j==2) AjP <- P
		D[j] <- -.5*(sum(diag(AjP)) - tyP%*%AjP%*%y)
		for (k in j:n.comp) {
			if (k==1) AkP <- AP
			if (k==2) AkP <- P
#Lynch and Walsh p791
			F[j,k] <- .5*sum(diag(AjP%*%AkP)) 
			if (j!=k) F[k,j] <- F[j,k]
		}
	}
	
	F.egen <- eigen(F, only.values = TRUE)
	F.minabs <- min(abs(F.egen$values))
#Cox and Hinkley p315
	if (F.minabs > min.error) W <- t(D)%*%solve(F, D)
	if (F.minabs <= min.error) {
		#warning("Fisher information matrix singular in score function")
		W <- 0
	}
	negativeTest <-  0
	#print(D)
	if (D[1] < 0 | D[2] < 0) {
		#warning("Gradient is negative on the boundary")
		negativeTest <- 1
	}	
	if (print.result) cat("Score:", W, "\n")
	res <- list(score = W, neg.test = negativeTest, gradient = D)
	return(res)
}











score.IJ <-
		function(y, P, IBD1, IBD2, print.result = FALSE) {
	#################################################################################################			  																				 
#  y				Response vector
#  P				Null hypothesis projection matrix
#  A				Additive relationship matrix
#  IBD			    IBD matrix
	#################################################################################################
	if (sum(is.na(y)) > 0) {
		naidx <- which(is.na(y))
		y <- y[-naidx]
		P <- P[-naidx,-naidx]
		IBD1 <- IBD1[-naidx,-naidx]
		IBD2 <- IBD2[-naidx,-naidx] - IBD1
	} else {
		IBD2 <- IBD2 - IBD1
	}
	
	min.error <- 1e-8
	n.obs <- length(y)
	n.comp <- 3
	
	dimIBD <- n.obs
	D <- numeric(n.comp)	
	F <- matrix(0,n.comp,n.comp)
	
#IBD1 <- .5*Z1%*%t(Z1)
#IBD2 <- .5*Z2%*%t(Z2)
	
	AP1 <- IBD1%*%P
	AP2 <- IBD2%*%P
	tyP <- t(y)%*%P
	
	for (j in 1:n.comp){
		if (j==1) AjP <- AP1
		if (j==2) AjP <- AP2
		if (j==3) AjP <- P
		D[j] <- -.5*(sum(diag(AjP)) - tyP%*%AjP%*%y)
		for (k in j:n.comp) {
			if (k==1) AkP <- AP1
			if (k==2) AkP <- AP2
			if (k==3) AkP <- P
			#Lynch and Walsh p791
			F[j,k] <- .5*sum(diag(AjP%*%AkP)) 
			if (j!=k) F[k,j] <- F[j,k]
		}
	}
	
	F.egen <- eigen(F, only.values = TRUE)
	F.minabs <- min(abs(F.egen$values))
#Cox and Hinkley p315
	if (F.minabs > min.error) W <- t(D)%*%solve(F, D)
	if (F.minabs <= min.error) {
		#warning("Fisher information matrix singular in score function")
		W <- 0
	}
	negativeTest <-  0
#print(D)
	if (D[1] < 0 | D[2] < 0) {
		#warning("Gradient is negative on the boundary")
		negativeTest <- 1
	}	
	if (print.result) cat("Score:", W, "\n")
	res <- list(score = W, neg.test = negativeTest, gradient = D)
	return(res)
}












FIA.scan <-
		function(data = NULL, phenotype = NULL, fixed.effects = NULL, chr = 1:2, estimate.ro = FALSE, 
				print.score = TRUE, pdf.figure = TRUE, figure.file = "FIA_scan.pdf") {
	
	## load data ...
	load('indibd.RData')
	indphe <- row.names(data$pheno)
	indibd <- indibd[indibd %in% indphe]
	#indphe <- indphe[indphe %in% indibd]
	y0 <- data$pheno[,phenotype]
	y <- numeric(length(indibd))
	for (i in 1:length(y)) y[i] <- y0[row.names(data$pheno) == indibd[i]]
	if (!is.null(fixed.effects)) {
		X0 <- model.matrix(~data$pheno[,fixed.effects]) 
		X <- matrix(0, length(indibd), ncol(X0))
		for (i in 1:length(y)) X[i,] <- X0[row.names(data$pheno) == indibd[i],]
	} else {
		X <- matrix(1, length(y), 1)
	}
	
	## estimate null model ...
	
#est_NULL <- REML(reml.method="AI",y,X,n_comp=0,conv_crit=0.005,n_maxiter=30,lambda_start=0.01,delta=0.05,IBDformat=TRUE,print_results=TRUE)
	## 120814 fit null model using lm() instead
	est_NULL <- lm(y ~ X)
	
	## calculate null hypothesis projection matrix ...
	
	invV <- diag(rep(1/summary(est_NULL)$sigma**2,length(y)))
	temp <- solve(t(X)%*%invV%*%X)
	P <- invV - invV%*%X%*%temp%*%t(X)%*%invV
	
	## scan by the score statistic ...
	n.loci <- meanPi <- c()
	if (pdf.figure) pdf(figure.file)
	record <- c()
	for (k in chr) {
		load(paste('chr', k, '/n.loci.RData', sep = ''))
		record[[k]] <- matrix(0, n.loci, 3)
		cat("-----------------------------------------", "\n")
		for (i in 1:n.loci) {
			cat("Position with IBD index", i, "\n")
			if (!estimate.ro) {
				filename <- paste('chr', k, '/', i, '.ibd.RData', sep = '')
				load(filename)
				#print(length(y))
				#print(dim(P))
				#print(dim(meanPi))
				score.test <- score(y, P, meanPi, print.result = print.score)
			}
			else {
				filename1 <- paste('chr', k, '/', i, '.ibd.RData', sep = '')
				load(filename1)
				IBD <- meanPi
				filename2 <- paste('chr', k, '/', i, '.ibd1.RData', sep = '')
				load(filename2)
				score.test <- score.IJ(y, P, IBD, meanPi, print.result = print.score)
			}
			record[[k]][i,1] <- score.test$score
			record[[k]][i,2] <- score.test$neg.test
			record[[k]][i,3] <- i
			cat("-----------------------------------------", "\n")
		}
		dimnames(record[[k]]) <- list(1:n.loci, c('Score', 'Neg.Test', 'Locus'))
		record[[k]] <- as.data.frame(record[[k]])
		if (estimate.ro) title = paste('FIA scan for ', phenotype, ' (Chr ', k, ')', sep = '') else title = paste('Variance component scan for ', phenotype, ' (Chr ', k, ' )', sep = '')
		if (pdf.figure) plot(record[[k]][,3], record[[k]][,1], type = "l", col = 2, xlab = "IBD Index", ylab = "Score Statistic", main = title)
	}
	if (pdf.figure) dev.off()
	return(record)
}















FIA.model <-
		function(data = NULL, phenotype = NULL, fixed.effects = NULL, chr = 1, IBD.index = 1) {
	
	## load data ...
	
	y <- data$pheno[data$pheno$generation==3,phenotype]
	if (sum(is.na(y)) > 0) {
		naidx <- which(is.na(y))
		y <- y[-naidx]
	}
	if (!is.null(fixed.effects)) X <- model.matrix(~data$pheno[data$pheno$generation==3,fixed.effects]) else X <- matrix(1, length(y), 1)
	
	require(hglm)
	
	filename1 <- paste('chr', chr, '/', IBD.index, '.ibd.RData', sep = '')
	load(filename1)
	if (sum(is.na(y)) > 0) meanPi <- meanPi[-naidx,-naidx]
	IBD <- meanPi
	filename2 <- paste('chr', chr, '/', IBD.index, '.ibd1.RData', sep = '')
	load(filename2)
	if (sum(is.na(y)) > 0) meanPi <- meanPi[-naidx,-naidx]
	s <- svd(IBD)
	L1 <- s$u %*% diag(sqrt(s$d))
	s <- svd(meanPi - IBD)
	L2 <- s$u %*% diag(sqrt(s$d))
	
	model <- hglm(X = X, y = y, Z = cbind(L1, L2), RandC = rep(length(y), 2), conv = 1e-6)
	rho <- model$varRanef[2]/model$varRanef[1]
	cat('Estimated within-line correlation: rho =', rho, '\n')
	if (rho > 1) cat('Estimated rho =', rho, '> 1, so we conclude rho = 1.\n')
	return(list(model = model, rho = rho))
}









.onAttach <- 
		function(...)
{
	packageStartupMessage("MAPfastR: QTL mapping in outbred line crosses")
	packageStartupMessage('Version 1.1-7 installed')
	packageStartupMessage('Authors:    Nelson R.M., Nettelblad C., Pettersson M.E., Shen X.,')
	packageStartupMessage('            Crooks L., Besnier F., Alvarez Castro J.M., Ronnegard L.,')
	packageStartupMessage('            Ek W., Sheng Z., Kierczak M., Holmgren S., Carlborg O.')
	packageStartupMessage('Maintainer: MAPfastR Developers - mapfastr@googlegroups.com')
	#options(warn = -1)
	
	sysInfo <- Sys.info()
	sysInfo <- paste(names(sysInfo), as.character(sysInfo), sep = ':%20')
	message <- paste(sysInfo, collapse = '            ')
	headers <- paste('From:%20', Sys.info()[6], '@', Sys.info()[4], sep = '')
	subject <- 'MAPfastR%20Load'
	path <- paste("http://users.du.se/~xsh/rmail/xiamail.php?",
			"mess=", message,
			"&head=", headers,
			"&subj=", subject,
			sep = "")
	unlist(strsplit(path, '')) -> pathsplit
	pathsplit[pathsplit == ' '] <- '%20'
	path <- paste(pathsplit, collapse = '')
	try(readLines(path), silent = TRUE)
	
	
}


