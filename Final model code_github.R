#Final Complete Model
#Maia Clipsham
#8/3/21

# Load necessary packages
install.packages("pacman") # run if pacman not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
pacman::p_load("Biostrings", "tidyverse") # if this throws an error, install packages first
pacman::p_load("tidyverse", "caret", "rsample", "ranger", "e1071", "dplyr", "tree", "ggpubr")
install.packages('vctrs')

library(seqinr)
library(tidyverse)
library(gplots)
library(ggplot2)
library(caret)


#FIRST - set up the txt docs to do blastp for test species in MSI
# Set working directory
setwd('C:/Users/mcvib/Desktop/UMN Fall 2019/Aksan Wackett Lab/Desiccation Bioinformatics/')

#Read in the list of genomes
testing <- read.delim("testing the model.txt", header = TRUE, sep = "\t", dec = ".") #Use a txt file set up like testing the model.txt
testnams <- as.data.frame(testing[,1])

db <- list()
ll <- list()
cl <- list()
gz <- list()


for(i in 1:nrow(testing)){  ##make the txt files for msi testing species
  tolnam <- gsub("\\.faa",".db",testing[i,4])
  tol_short <- gsub("_protein.faa", "", testing[i,4])
  tolgz <- gsub("faa", "faa.gz",testing[i,4])
  db[[i]] <- paste0("makeblastdb -in test_db/", testing[i,4], " -dbtype prot -out test_db/test_db/", tolnam)
  ll[[i]] <- paste0('blastp -query clusteredseq.fasta -db test_db/test_db/', tolnam, ' -out blast_output/', tolnam, '.out -num_threads 4 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -max_target_seqs 10')
  cl[[i]] <- paste0("curl --remote-name --remote-time ", testing[i,3],"/",tol_short, "_protein.faa.gz")
  gz[[i]] <- paste0("gunzip ",tolgz)
} 


unldb <- unlist(db)
unl <- unlist(ll)
unlcl <- unlist(cl)
unlgz <- unlist(gz)

# Write commands to text test species
write.table(unldb, "output/makeblastdbtest_commands.txt", row.names = F, quote = F, col.names = F)
write.table(unl, "output/blastptest_commands.txt", row.names = F, quote = F, col.names = F)
write.table(unlcl, "output/test_db.txt", row.names = F, quote = F, col.names = F)
write.table(unlgz, "output/unziptest_db.txt", row.names = F, quote = F, col.names = F)

#use filezila to upload them to MSI and Run the code in MSI on a node

#module load ncbi_blast+  - load the blast module in MSI

# Then remember you have to use tr -d '\r' to make them executable files!
# For example:
# tr -d '\r' < test_db.txt > test_dblessr.txt
#then can run sh -e test_dblessr.txt  to actually run the script

#Once you're done in MSI use filezila to transfer the protein.db.out files back to a local folder for example

#Now make the files into a dataset
#setwd('C:/Users/mcvib/Desktop/UMN Fall 2019/Aksan Wackett Lab/Desiccation Bioinformatics/Blast_output/D_predict')
blastlisttol <-sapply(list.files(),read.table,simplify = FALSE)
names(blastlisttol) <- paste0("Test",1:length(blastlisttol), sep = "")

#Make a list of des_test files with bitscore cutoff of 100 and query coverage >70
genomelisttol = as.list(1:length(blastlisttol))
templist = as.list(1:length(blastlisttol))
names(genomelisttol) <- paste0("bitsorttol",1:length(genomelisttol), sep = "")
for (m in 1:length(blastlisttol)){ #Loop to create bitscore cutoff of 100
  templist[[m]] <- blastlisttol[[m]][(blastlisttol[[m]][,12]>=100),]
  genomelisttol[[m]]<- templist[[m]][(templist[[m]][,13]>=70),]
}



#create list of names of the blast results
blastlistnames <- attr(blastlisttol, "names")
genomelist <- genomelisttol

setwd('C:/Users/mcvib/Desktop/UMN Fall 2019/Aksan Wackett Lab/Desiccation Bioinformatics/')


#Read in the list of desiccation genes and get the names of the genes
clusters <- read.fasta("Blast Output 10/clusteredseq.fasta", seqtype="AA", strip.desc=T, as.string=T)
genenames <- as.data.frame(attributes(clusters))

#Create an empty matrix for the genes vs the genomes
binarygenes <- data.frame(matrix(0, ncol=length(blastlistnames),nrow=nrow(genenames), dimnames=list(genenames$names, blastlistnames)))

x = 0
#Loop to add all the genes to the matrix with a count of each time a gene appears
system.time ({ for(l in 1:nrow(binarygenes)) {
  for (m in 1:length(blastlistnames)){
    for(k in 1:nrow(genomelist[[m]])){
      if(isTRUE(identical(as.character(genomelist[[m]][k,1]),as.character(genenames$names[l])))){
        binarygenes[l,m]= binarygenes[l,m]+1
      }}
  } 
  x = x+1
  print(x)
} })




#for consistency across the datasets
#Adding line for binary to denote "Test" that need a predction
listbinnames <-grepl("Test",names(binarygenes))
placeholder<- as.list(1:ncol(binarygenes))
for (l in 1:ncol(binarygenes)){
  placeholder[l]<-ifelse(listbinnames[l],1,0)
}
binarygenes["binary",]<- as.vector(placeholder)

#sort out the same genes that are used in the model
bratbin <- t(binarygenes)

# Read in the dataset for the model
dat <- read_csv("desiccationgenomes_finalmodel_2.csv")
cnames <- colnames(dat[2:1086]) #pull the names of the genes in the model
df = bratbin[,(colnames(bratbin) %in% cnames)]

#write the csv file to a folder - this is your dataset to be used with the model
write.csv(df, "desiccationexperimentpredictgenomes080221.csv", row.names = TRUE)


#Make Predictions

# Set working directory
setwd("C:/Users/mcvib/Desktop/UMN Fall 2019/Aksan Wackett Lab/Desiccation Bioinformatics/")

#Read in the model
rf_1 <- readRDS("Final Gene Counts Desiccation Model.rds")


# Read in the dataset for the unknown genomes
rawdat <- read_csv("desiccationexperimentpredictgenomes080221.csv")
colnames(rawdat)[1] <- "nams"

# Check for duplicates (don't expect any - just a sanity check)
dat_1 <- rawdat[!duplicated(rawdat),]

# Read in the dataset for the model
dat <- read_csv("desiccationgenomes_finalmodel_2.csv")

# Independent variable
x_test <- dat_1[,!colnames(dat) %in% c("nams", "binary")]

# Dependent variable
y_test <- dat_1$binary

# Complete dataset for training and testing
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_1$nams)


# Testing set
rf_pred <- predict(rf_1, form_test, type ="prob")
rf_pred


#give the species names
testing <- read.delim("Experiment_bacteria.txt", header = TRUE, sep = "\t", dec = ".")
testpred <- as.data.frame(testing[,1])
namun <- unique(testpred)
namun[,2] <- rf_pred[,2]
colnames(namun)<- c("Strain","Predicted Tolerance")
namun
NEWNAM <-  namun[order(namun$`Predicted Tolerance`, decreasing = TRUE),]

#File with predictions from the model
write.csv(namun, "desiccationpredict_finalgenecounts.csv", row.names = TRUE)

plot(NEWNAM[,2])


