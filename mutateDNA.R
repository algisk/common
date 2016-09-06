###########################
### 	Environment 	###
###########################
require(foreach)
require(data.table)

args=commandArgs(trailingOnly=TRUE)

# Usage:
# Rscript --vanilla mutateDNA.R arg1 arg2 arg3 arg4
# arg1 - DNA fasta
# arg2 - mutation file
# arg3 - mutation number
# arg4 - codon to change into
###########################
### 		Body 		###
###########################

# Read DNA sequence
DNA <- readLines(args[1])[2]

# Split sequence to codons
DNA3 <- data.table(
	sapply(strsplit(DNA, "(?<=[[:alpha:]]{3})", perl = TRUE),
         function(x) paste(x)))

# Read mutation list
muts <- fread(args[2])
nMut <- args[3]
muts <- muts$Mut[1:nMut]
muts <- as.numeric(gsub("[[:alpha:]]", "", muts))

# Replace AAs and save new fasta
DNA3[muts] <- args[4]
DNA <- paste(DNA3$V1, collapse="")
Name <- readLines(args[1])[1]
sName <- gsub("[[:punct:]]", "", Name)

write(rbind(Name,DNA), paste0(sName,"_",nMut,".fasta"))