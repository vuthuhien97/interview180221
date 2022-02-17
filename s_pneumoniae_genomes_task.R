##Setup
install.packages("seqinr")
install.packages("tidyverse")
installed.packages("data.table")
install.packages("kmer")
install.packages("digest")
library(tidyverse)
library(seqinr)
library(kmer)
library(digest)


#####Reading kmer
setwd("C:/Users/elvin/Desktop/a")
R6 <- read.fasta(file = "s_pneumoniae_genomes/R6.fa",as.string = T, seqtype = c("DNA"))
TIGR4 <- read.fasta(file = "s_pneumoniae_genomes/TIGR4.fa",as.string = T, seqtype = c("DNA"))
contigs82<- read.fasta(file = "s_pneumoniae_genomes/14412_3#82.contigs_velvet.fa",s.string = T, seqtype = c("DNA"))
contig84<- read.fasta(file = "s_pneumoniae_genomes/14412_3#84.contigs_velvet.fa", s.string = T, seqtype = c("DNA"))
 
###Run Though the substrings
short = substr(R6, 1, 30)
kmers = substring(short, seq(1,nchar(short)-13), seq(14,nchar(short)))
kcount 
length(kmers)

splitted = unlist(strsplit(kmers[[1]], ""))
splitted[4]


#### Count the occurences of 14kmers

kmers_R6 = substring(R6, seq(1,nchar(R6)-13), seq(14,nchar(R6)))
nchar(R6)

kmers_TIGR4 = substring(TIGR4, seq(1,nchar(TIGR4)-13), seq(14,nchar(TIGR4)))
nchar(TIGR4)

#define Jaccard Similarity function
jaccard <- function(kmers_R6, kmers_TIGR4) {
  intersection = length(intersect(kmers_R6, kmers_TIGR4))
  union = length(kmers_R6) + length(kmers_TIGR4) - intersection
  return (intersection/union)
}

#find Jaccard Similarity between the two sets 
jaccard(kmers_R6, kmers_TIGR4)




# MinHash Jaccard R6
hash1_R6 <- digest(kmers_R6[4], algo = "murmur32", serialize = F, seed = 0)
hash1_R6

# Run on all the k-mers
hash_run1_R6 <- unlist(lapply(kmers_R6, digest, algo = "murmur32", serialize = F, seed = 0))
hash_run1_R6

hash_run2_R6 <- unlist(lapply(kmers_R6, digest, algo = "murmur32", serialize = F, seed = 0))
all(hash_run1_R6 == hash_run2_R6)


# MinHash Jaccard TIGR4
hash1_TIGR4 <- digest(kmers_TIGR4[4], algo = "murmur32", serialize = F, seed = 0)
hash1_TIGR4


# Run on all the k-mers
hash_run1_TIGR4 <- unlist(lapply(kmers_TIGR4, digest, algo = "murmur32", serialize = F, seed = 0))
hash_run1_TIGR4

hash_run2_TIGR4 <- unlist(lapply(kmers_TIGR4, digest, algo = "murmur32", serialize = F, seed = 0))
all(hash_run1_TIGR4 == hash_run2_TIGR4)
