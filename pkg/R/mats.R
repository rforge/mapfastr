#sample data import and data frame constructor
import_data <- function(pedigree_file = "pedigree.txt", phenotype_file = "phenotypes.txt", genotype_file = "marker_data.txt", marker_file = "marker_info.txt", backcross = FALSE, backcross.line = NA, backcross.parent = NA, heterogam = 1, sex.chrom = "last", sex.restrict = FALSE, n_gen_map = 2, ref.map = "average", format = "cnf2freq", crimap_data_dir = NA){
	
	#Bookkeeping operations
	if (backcross){
		backcross <- 1
	}
	else{
		backcross <- 0
	}
	
	if (sex.restrict){
		sex.restrict <- 1
	}
	else{
		sex.restrict <- 0
	}
	
	#Format dependent data import 
	if (format ==  "cnf2freq"){
		#Pedigree and phenotypes
		ped_data <- read.table(pedigree_file, fill=TRUE, stringsAsFactors=FALSE)
		pedigree <- data.frame(generation = 0, sex = 0, parent_1 = 0, parent_2 = 0, line = 0)
		#phenotypes <- read.table(phenotype_file, stringsAsFactors=FALSE, header = TRUE)
		#Pedigree parsing
		ind <- 1
		pos <- 0
		individuals <- c(0)
		for (i in 1:dim(ped_data)[1]){
			if(!is.na(ped_data[i,2])){
				pos <- pos + 1
				if (!(ped_data[i,1] %in% individuals)){
					if (pos < 5){
						generation <- 1
					}
					else if (pos < 7){
						generation <- 2
					}
					else{
						generation <- 3
					}
					pedigree[ind,1] <- generation
					pedigree[ind,2:5] <- ped_data[i,c(4,2,3,5)]
					#pedigree[ind,2] <- 2/pedigree[ind,2]
					individuals[ind] <- ped_data[i,1]
					ind <- ind + 1
				}
			}
			else{
				pos <- 0
			}
		}
		#Adjusting the sex column 
		#if(0 %in% pedigree[,2]){
		#	pedigree[,2] <- pedigree[,2] + 1
		#	print("Adjusted sex indicator from (0,1) to (1,2).")
		#}
		#pedigree[,2] <- 2/pedigree[,2]
		
		row.names(pedigree) <- individuals
		print("Pedigree parsed.")
		#Phenotypes
		#pedigree <- pedigree[order(row.names(pedigree)),]
		#phenotypes <- phenotypes[order(as.character(phenotypes[,1])),]
		#j <- 1
		#for (i in 2:length(names(phenotypes))){
		#	pedigree[row.names(pedigree) %in% phenotypes[,1],(5+j)] <- phenotypes[,i]
		#	names(pedigree)[5+j] <- names(phenotypes)[i]
		#	j <- j +1
		#}
		#Custom NA insertion
		#if (sum(pedigree == -100001, na.rm = T) > 0){
		#	pedigree[pedigree == -100001] <- NA
		#	print("-100001 changed to NA.")
		#}
		#print("Phenotypes parsed.")
		
		#Genotypes and marker info
		#Preprocessing
		marker_n_detect <- read.table(marker_file, fill = TRUE, stringsAsFactors = FALSE,nrows = 1)
		marker_n_detect <- read.table(marker_file, fill = TRUE, stringsAsFactors = FALSE, skip = 1, nrows = marker_n_detect[1,1], col.names = 1:marker_n_detect[1,2])
		marker_n <- max(marker_n_detect[,1])
		marker_inf <- read.table(marker_file, fill = TRUE, stringsAsFactors = FALSE, col.names = 1:(marker_n+1))
		if (!is.na(sex.chrom)){
			if( sex.chrom == "last"){
				sex.chrom <- as.numeric(marker_inf[1,1])
			}
		}
		#Foring first col of marker_inf to be numeric.
		marker_inf[,1] <- as.numeric(marker_inf[,1])
		no_used_mrk <- sum(marker_inf[2:(marker_inf[1,1] +1),1])
		geno <- read.table(genotype_file, colClasses = "integer")
		geno <- t(geno)
		genotypes <- data.frame(row.names = c(1:no_used_mrk))
		
		j <- 1
		line <- 2
		for (i in line:(line+marker_inf[1,1]-1)){
			genotypes[j:(j+marker_inf[i,1]-1),"col_num"] <- as.numeric(marker_inf[i,2:(marker_inf[i,1]+1)])
			j <- j+marker_inf[i,1]
		}
		line <- line + marker_inf[1,1]
		print("Genotypes preprocessed.")
		#Genotype parsing
		even <- genotypes[, "col_num"]*2
		odd <- even + 1
		for ( i in 1: dim(geno)[2]){
			genotypes[,(2*i)] <- geno[even, i]
			genotypes[,(2*i)+1] <- geno[odd, i]
		}
		genotypes <- genotypes[,2:dim(genotypes)[2]]
		genotypes[genotypes == 0 | genotypes == 9] <- NA
		col_name_vec <- rep(geno[1,], each = 2)
		print("Genotypes parsed.")
		#Marker info parsing
		if ( n_gen_map == 2){
			j <- 1
			chr <- 1
			for (i in seq(line,(line+(marker_inf[1,1]*2-1)),2)){
				genotypes[j:(j+marker_inf[chr+1,1]-1),"chr"] <- chr
				genotypes[j:(j+marker_inf[chr+1,1]-1),"sex_1_cM"] <- signif(cumsum(as.numeric(marker_inf[i + 1,2:(marker_inf[chr+1,1]+1)])), digits = 4)
				genotypes[j:(j+marker_inf[chr+1,1]-1),"sex_2_cM"] <- signif(cumsum(as.numeric(marker_inf[i,2:(marker_inf[chr+1,1]+1)])), digits = 4)
				j <- j+marker_inf[chr+1,1]
				chr <- chr + 1
			}
			line <- line + marker_inf[1,1]*2
		}
		if ( n_gen_map == 1){
			j <- 1
			chr <- 1
			for (i in seq(line,(line+(as.numeric(marker_inf[1,1])-1)),1)){
				genotypes[j:(j+marker_inf[chr+1,1]-1),"chr"] <- chr
				genotypes[j:(j+marker_inf[chr+1,1]-1),"sex_1_cM"] <- signif(cumsum(as.numeric(marker_inf[i,2:(marker_inf[chr+1,1]+1)])), digits = 4)
				genotypes[j:(j+marker_inf[chr+1,1]-1),"sex_2_cM"] <- signif(cumsum(as.numeric(marker_inf[i,2:(marker_inf[chr+1,1]+1)])), digits = 4)
				j <- j+marker_inf[chr+1,1]
				chr <- chr + 1
			}
			line <- line + marker_inf[1,1]
			print("Duplicated single genetic map.")
		}
		if ( ref.map == "average"){
			genotypes[,"ref_cM"] <- rowMeans(genotypes[,c("sex_1_cM","sex_2_cM")])
		}
		if ( ref.map == "sex_1"){
			genotypes[,"ref_cM"] <- genotypes[,"sex_1_cM"]
		}
		if ( ref.map == "sex_2"){
			genotypes[,"ref_cM"] <- genotypes[,"sex_2_cM"]
		}
		
		name_vec <- c(0)
		chr <- 1
		j <- 1
		for (i in line:(line+(marker_inf[1,1]-1))){
			name_vec[j:(j+marker_inf[chr + 1,1]-1)] <- unlist(marker_inf[i,2:(marker_inf[chr+1,1]+1)])
			j <- j+marker_inf[chr + 1,1]
			chr <- chr + 1
		}
		names(genotypes)[1:(dim(genotypes)[2]-4)] <- col_name_vec
		row.names(genotypes) <- name_vec
		print("Marker info parsed.")
	}
	
	if(format == "crimap"){
		#Pedigree and phenotypes
		pedigree <- data.frame(generation = 0, sex = 0, parent_1 = 0, parent_2 = 0, line = NA)
		
		if(!file.exists(crimap_data_dir)){
			print("A directory containing the crimap data files must be supplied.")
			return(NULL)
		}
		#All pedigree files must have suffix .gen
		pedigree_files <- grep(pattern="[.]gen",x=dir(crimap_data_dir, full.names=T), value = T)
		#All map files must have suffix .map
		map_files <- grep(pattern="[.]map",x=dir(crimap_data_dir, full.names=T), value = T)
		chr_numbers <- as.numeric(sub(pattern = "[Cc]hr_([0-9]+)[.].+", replacement = "\\1", x = map_files))
		
		pedigree_files <- pedigree_files[order(chr_numbers)]
		map_files <- map_files[order(chr_numbers)]
		print("Input files ordered.")
		
		#Detection of the sex chromosome
		if (sex.chrom == "last"){
			sex.chrom <- chr_numbers[length(chr_numbers)]
		}
		
		#Build up the pedigree from the first data file
		ped_data <- read.table(pedigree_files[1], fill=TRUE, stringsAsFactors=FALSE, skip = 3)
		ind <- 1
		pos <- 1
		individuals <- c(0)
		for (i in 1:dim(ped_data)[1]){
			if(!is.na(ped_data[i,2])){
				pos <- pos + 1
				if (!(ped_data[i,1] %in% individuals)){
					if (ped_data[i,2] == 0 & ped_data[i,3] == 0 ){
						generation <- 1
						pos <- pos - 1
					}
					else if (pos <= 2){
						generation <- 2
					}
					else{
						generation <- 3
					}
					pedigree[ind,1] <- generation
					pedigree[ind,2:4] <- ped_data[i,c(4,2,3)]
					individuals[ind] <- ped_data[i,1]
					ind <- ind + 1
				}
			}
			else{
				pos <- 0
			}
		}
		pedigree[pedigree[,2] == 3,2] <- NA
		print("Pedigree processed.")		
		
		#Collect genotypes from each chromosome
		chr_counter <- 1
		for(ped_file in pedigree_files){
			ped_data <- read.table(ped_file, fill=TRUE, stringsAsFactors=FALSE, skip = 1)
			ped_data <- t(ped_data)
			map_data <-read.table(map_files[chr_counter], fill = T, stringsAsFactors = F, col.names = c(1:12))
			map_data <- map_data[min(which(map_data[,1] == 0)):(dim(map_data)[1] -2),1:3]
			map_data <- map_data[map_data[,3] != "",]
			no_markers <- as.numeric(ped_data[1,1])
			chr_genotypes <- data.frame(row.names = c(1:no_markers))
			odd <- seq(from = 1, to = no_markers*2, by = 2)
			even <- odd + 1
			counter <- 1
			for (ind in individuals){
				ind_index <- ped_data[1,]
				if(sum(ind_index == ind) >= 1){
					ind_geno <- ped_data[5:dim(ped_data)[1],min(which(ind_index == ind))]
					chr_genotypes[, counter] <- ind_geno[odd]				
					chr_genotypes[, (counter+1)] <- ind_geno[even]
					names(chr_genotypes)[counter:(counter+1)] <- ind
					counter <- counter + 2
				}
			}
			chr_genotypes[,"chr"] <- chr_numbers[chr_counter]
			chr_genotypes[,"sex_1_cM"] <- map_data[,3]
			chr_genotypes[,"sex_2_cM"] <- map_data[,3]
			chr_genotypes[,"ref_cM"] <- map_data[,3]
			
			if(chr_counter == 1){
				genotypes <- chr_genotypes
			}
			else{
				genotypes <- rbind(genotypes,chr_genotypes)
			}
			chr_counter <- chr_counter + 1
		}
		print("Genotypes processed.")
		
	}
	
	#Post-processing
	#Adjusting the sex column 
	if(0 %in% pedigree[,2]){
		pedigree[,2] <- pedigree[,2] + 1
		print("Adjusted sex indicator from (0,1) to (1,2).")
	}
	pedigree[,2] <- 2/pedigree[,2]
	row.names(pedigree) <- individuals
	print("Pedigree post-processed.")
	
	#Ordering the pedigree and adding phenotypes
	phenotypes <- read.table(phenotype_file, stringsAsFactors=FALSE, header = TRUE)
	pedigree <- pedigree[order(row.names(pedigree)),]
	phenotypes <- phenotypes[order(as.character(phenotypes[,1])),]
	j <- 1
	for (i in 2:length(names(phenotypes))){
		pedigree[row.names(pedigree) %in% phenotypes[,1],(5+j)] <- phenotypes[,i]
		names(pedigree)[5+j] <- names(phenotypes)[i]
		j <- j +1
	}
	#Custom NA insertion
	if (sum(pedigree == -100001, na.rm = T) > 0){
		pedigree[pedigree == -100001] <- NA
		print("-100001 changed to NA.")
	}
	print("Phenotypes parsed.")
	
	#Return statement
	res <- list(pheno = pedigree, geno = genotypes, backcross = backcross, backcross.line = backcross.line, backcross.parent = backcross.parent, sex.chrom = sex.chrom, heterogam = heterogam, sex.restrict = sex.restrict)
	class(res) <- 'MAPfastR'
	return(res)
}








#MAPfastR to NOIA format converter
run_NOIA <- function(data, phenotype, markers, noia_func = "linear"){
	require("noia")
	phen <- data$pheno[,phenotype]
	gen <- t(data$geno[markers,!grepl("cM|chr", x = names(data$geno))])
	temp_gen <- gen
	ed_gen <- gen[seq(from = 1, to = dim(gen)[1], by = 2),]
	for (column in 1:dim(gen)[2]) {
		genotypes <- unique(gen[,column])
		genotypes <- genotypes[!genotypes %in% NA]
		if(NA %in% genotypes & length(genotypes) > 3 | !NA %in% genotypes & length(genotypes) > 2){
			print(paste("Error at locus, ", names(gen)[column], ": only biallelic loci allowed.", sep = ""))
		}
		else{
			print("Locus accepted.")
			print(dim(gen)[1])
			print(dim(ed_gen)[1])
			temp_gen[gen[,column] == genotypes[1],column] <- 1
			temp_gen[gen[,column] == genotypes[2],column] <- 3
			ed_gen[, column] <- temp_gen[seq(from = 1, to = dim(gen)[1], by = 2),column]
			ed_gen[temp_gen[seq(from = 1, by = 2, to = dim(gen)[1]),column] != temp_gen[seq(from = 2, by = 2, to = dim(gen)[1]),column], column] <- 2
		}
	}
	if(noia_func == "multilinear"){
		result <- multilinearRegression(phen = phen, gen = ed_gen)
	}
	if(noia_func == "linear"){
		result <- linearRegression(phen = phen, gen = ed_gen)
	}
	if(noia_func == "none"){
		result <- ed_gen
	}
	return(result)	
}