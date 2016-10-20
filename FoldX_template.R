###########################
### 	Environment 	###
###########################

library(foreach)
library(data.table)
library(ggplot2)

args=commandArgs(trailingOnly=TRUE)
nRuns <- as.numeric(args[2])

###########################
### 	FoldX 			###
###########################

###########################
# to execute type: 
# Rscript --vanilla FoldX_template.R 
#		[filename.pdb] [nRuns] [file of aa's to omit]
# *.pdb must be in same directory as foldx executable
# working dir must contain foldx executable and rotbase.txt
###########################
## Insert downloaded pdb
pdbName <- args[1]

repDir <- paste("./", gsub(".pdb", "", pdbName), "_repair", sep="")
system(paste0("mkdir -p ", repDir))

## Repair PDBs
system(paste0("./foldx --pdb=", pdbName, " --command=RepairPDB --water=CRYSTAL --output-dir=", repDir, sep=""))
repFile <- list.files(repDir, pattern='Repair.pdb')
file.copy(paste0(repDir, "/", repFile), repFile, overwrite=TRUE)

###########################
### 	Mutations 		###
###########################

G <- c("Y", "I", "M", "W", "L")
G3 <- c("TYR", "ILE", "MET", "TRP", "LEU")
GG <- structure(G3, .Names=G)

fname <- repFile
system(paste0("sed -i 's/ 1[.]00/ 1.00 /g' ",fname))
nLines <- grep("TER", readLines(fname))[1] -1
if (length(which(nLines > 0)) != 0) {
system(paste0('sed -i "s/A1/A /g" ',fname))
system(paste0('head -n',nLines, ' ', fname, ' > temp.pdb'))
system(paste0('cat temp.pdb > ', fname))
system(paste0('rm temp.pdb'))
	temp <- fread(fname)
	# x <- gsub("A1", "A ", temp$V5)
	# x <- data.table(do.call('rbind', strsplit(x, " ")))
	# names(x) <- c("V5", "V6")
	# temp <- cbind(temp[,c(1:4), with=FALSE], x, temp[,c(6:10), with=FALSE])
	# names(temp) <- c("V1", "V2", "V3", "V4", "V5", 
	# 		"V6", "V7", "V8", "V9", "V10", "V11")
	# temp$V6 <- as.numeric(temp$V6)
	# write(temp, fname)
} else { 
	temp <- fread(fname)
}


res <- foreach(I=GG, .combine=rbind) %do% {
	idx <- which(temp$V4==I)
	if (length(idx) > 0) {
	mut <- paste(names(which(GG==I)), "A", temp$V6[idx], "F", sep="")
	mut <- unique(mut)
	data.table(mut)
	}
}

if (length(which(res$mut=="MA1F")) != 0) {
	idx <- which(res$mut=="MA1F")
	res <- res[-idx] 
}

if (length(args) >= 3) {
	omits <- as.numeric(fread(args[3]))
	pos <- as.numeric(gsub("[[:alpha:]]", "", res$mut))
	idx <- which(pos %in% omits)
	res <- res[-idx]
}

res <- res$mut

getSummary <- function(data, ID) {
	data.table(data[1], Run=ID)
}

getFile <- function(data, pdbfile, dirDir) {
	temp <- fread(paste0(dirDir, "/Raw_",pdbfile,".fxout", sep=""))
	idx <- grep(data[1]$Pdb, temp$Pdb)[1:nRuns]
	t <- temp[idx]
	names(t) <- gsub(" ", "_", names(t))
	t[order(total_energy),][1]$Pdb
}
## Single mutation

resDir <- paste("./", gsub(".pdb", "", pdbName), "_result", sep="")
PDB <- gsub(".pdb", "", pdbName)

write(paste(res, ";", sep=""), "individual_list.txt")
system(paste0("mkdir -p ",resDir))
system(paste0("./foldx --pdb=", fname, " --command=BuildModel --mutant-file=individual_list.txt --water=CRYSTAL --numberOfRuns=",nRuns, " --output-dir=", resDir, " --out-pdb=1", sep=""))
# Read result and select top 2
dx <- fread(paste0(resDir, "/Average_", PDB,"_Repair.fxout", sep=""))
dx <- dx[,1:3, with=FALSE]
names(dx) <- gsub(" ", "_", names(dx))
dx$Mut <- res
dx <- dx[order(total_energy),]
topMut <- dx[1:2]
fname <- getFile(topMut[1], paste0(PDB, "_Repair"), resDir)
system(paste0("cp ", resDir, "/",fname," ", PDB,"_One.pdb", sep=""))

## First mutation path test
nres <- res

path_One <- foreach(Run=seq(1:length(res)), .combine=rbind) %do% {
	fname <- paste0(PDB,"_One.pdb")
	system(paste0("sed -i 's/ 1[.]00/ 1.00 /g' ",fname))
	if (which(nres==topMut$Mut[1])>=0) {
		nres <- nres[-which(nres==topMut$Mut[1])]
	} 

	if (length(nres) != 0) {
	write(paste(nres, ";", sep=""), "individual_list.txt")
	system(paste0("mkdir -p ",PDB, "_One"))
	outDir <- paste0("./",PDB,"_One")
	system(paste0("./foldx --pdb=", fname, " --command=BuildModel --mutant-file=individual_list.txt --water=CRYSTAL --numberOfRuns=",nRuns, " --output-dir=", outDir, " --out-pdb=1", sep=""))
	# Read nresult and select top 2
	dx <- fread(paste0(outDir, "/Average_",PDB,"_One.fxout", sep=""))
	dx <- dx[,1:3, with=FALSE]
	names(dx) <- gsub(" ", "_", names(dx))
	dx$Mut <- nres
	dx <- dx[order(total_energy),]
	topMut <- dx[1]
	filename <- getFile(topMut[1], paste0(PDB,"_One"), outDir)
	system(paste0("cp ", outDir, "/",filename," ",PDB, "_One.pdb", sep=""))
	system(paste0("rm ", outDir, "/*", sep=""))
	getSummary(topMut, Run)
	}
}

write.csv(path_One, paste0(PDB,"_", Run, "_mut.csv"), row.names=FALSE)

###########################
### 	Plot energy 	###
###########################

## Load data
WT <- fread(paste0(resDir,"/Raw_",PDB, "_Repair.fxout"))
mut50 <- fread(paste0(PDB,"_", Run, "_mut.csv"))

## Get WT energy
WT <- WT[grep("WT", WT$Pdb)]
names(WT) <- gsub(" ", "_", names(WT))
energy <- mean(WT$total_energy)

## Get ene diff
cene <- energy

mut50$WT_mut <- foreach(I=1:nrow(mut50), .combine=rbind) %do% {
	cene <- cene + mut50$total_energy[I]
	cene
}

write.csv(mut50, paste0(PDB,"_", Run, "_mut.csv"), row.names=FALSE)

## Plot
png(paste0(PDB,"_energy.png"), height=700, width=700, pointsize=16)
# plot(mut50$WT_mut)
# abline(h=energy, col='red', lwd=2)
# abline(v=which(mut50$WT_mut > energy)[1], col='blue', lwd=2)
# text(mut50$WT_mut, labels=1:nrow(mut50), cex=0.7, pos=3)
ggplot() +
	geom_point(data=mut50, aes(x=Run, y=WT_mut)) +
	geom_hline(aes(yintercept=energy, color="#FF0000")) +
	geom_point(data=mut50, aes(x=which(mut50$WT_mut > energy)[1], 
		y=mut50$WT_mut[which(mut50$WT_mut > energy)[1]],
		color="#FF0000",
		shape="circle"), 
		size=2,
		stat="identity") +
	scale_shape(solid=TRUE, guide=FALSE) +
	scale_color_manual(values=c("#FF0000"), guide=FALSE) +
	geom_text(aes(x=which(mut50$WT_mut > energy)[1],
		y=mut50$WT_mut[which(mut50$WT_mut > energy)[1]],
		label=which(mut50$WT_mut > energy)[1]), 
		vjust=-2) +
	theme_bw(base_size = 12) +
	theme(axis.line = element_line(size=1, color = "black")) +
	theme(text=element_text(size=12, family="Raleway Light")) +
	xlab("Number of mutations") +
	ylab("dG")
dev.off()

x <- mut50$total_energy/mut50$Run

a <- foreach(I=2:(length(x)-1), .combine=rbind) %do% {
	if (x[I-1] > x[I] & x[I] < x[I+1]) {
		I
	}
} 

png(paste0(PDB,"_pits.png"), height=700, width=700, pointsize=16)
# plot(mut50$total_energy/mut50$Run)
# abline(v=a, lwd=2)
# text(mut50$total_energy/mut50$Run, labels=1:nrow(mut50), cex=0.7, pos=3)
ggplot() +
	geom_point(data=mut50, aes(x=Run, 
		y=mut50$total_energy/mut50$Run)) +
	geom_point(aes(x=a, y=mut50$total_energy[a]/mut50$Run[a],
		shape="circle",
		color="#FF0000"),
	size=2) +
	scale_shape(solid=TRUE, guide=FALSE) +
	scale_color_manual(values=c("#FF0000"), guide=FALSE) +
	geom_text(aes(x=a, y=mut50$total_energy[a]/mut50$Run[a],
		label=mut50$Run[a]),
		vjust=2) +
	theme_bw(base_size = 12) +
	theme(axis.line = element_line(size=1, color = "black")) +
	theme(text=element_text(size=12, family="Raleway Light")) +
	xlab("Number of mutations") +
	ylab("ddG / Number of mutations")
dev.off()

system(paste0("mkdir -p out_",PDB))
system(paste0("mv ",PDB ,"* out_",PDB))


###########################
### 	Get .fasta	 	###
###########################
### For manual use ###
# repFile <- repaired PDB (WT)
# PDB <- "example" (not "example.pdb")
# mut50 <- fread("example_n_mut.csv")
# create arguments a and x from previous section
infile <- fread(repFile)

sequence <- unique(infile[,c(4,6), with=FALSE])[,V4]
sequence <- paste(sequence, collapse="")

convert <- function(l) {

  map <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
           "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

  names(map) <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN",
                  "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                  "PRO", "SER", "THR", "TRP", "TYR", "VAL")

  sapply(strsplit(l, "(?<=[A-Z]{3})", perl = TRUE),
         function(x) paste(map[x], collapse = ""))
}

newseq <- convert(sequence)

fastav <- foreach(I=1:nchar(newseq), .combine=rbind) %do% {
	t <- substr(newseq, I, I)
	data.table(WT=t)
}

rownames(fastav) <- unique(infile$V6)

a <- as.numeric(a)
Mut <- as.numeric(gsub("[[:alpha:]]", "", mut50$Mut))

muts <- foreach(I=a, .combine=cbind) %do% {
	idx <- Mut[1:I]
	sidx <- match(idx, as.numeric(rownames(fastav)))
	t <- c(fastav$WT)
	t[sidx] <- "F"
	data.table(t)
}

names(muts) <- paste("m_",a,sep="")
fastav <- cbind(fastav, muts)
fastav <- t(fastav)

foreach(I=1:nrow(fastav), .combine=cbind) %do% {
	fr <- rownames(fastav)[I]
	A <- paste(">",PDB,"-",fr, sep="")
	# B <- paste(">1CZD-",fr,":B|PDBID|CHAIN|SEQUENCE", sep="")
	# C <- paste(">1CZD-",fr,":C|PDBID|CHAIN|SEQUENCE", sep="")
	t1 <- paste(fastav[I,1:ncol(fastav)], collapse="")
	# t2 <- paste(fastav[I,81:160], collapse="")
	# t3 <- paste(fastav[I,161:228], collapse="")

	write(rbind(A,t1), paste(PDB,"_",fr,".fasta",sep=""))
}

system(paste0("mkdir -p out_",PDB,"/fastas"))
system(paste0("mv *.fasta out_",PDB,"/fastas"))
