##SCRIPT TO DISCRETISE NORMALISED MICROARRAY DATA ACCORDING TO RELATIVE GENE EXPRESSION##
##AUTHOR: Tamara Lechon Gomez
##ADDRESS: lechongomezt@cardiff.ac.uk
##DATE CREATED: 24 Aug 2024
##DATE LAST MODIFIED: 26 Sep 2024

# Load required libraries
library(dplyr)

# Read input files 

#Expected input file 1 is a csv file with tab separated columns containing normalised expression values
#Each row corresponds to a probeset/gene; each column is a separate experiment
normalised.input <- read.csv(file.choose(), sep = '\t', header = T)
normalised.df <- normalised.input[2:2374]

#Expected input file 2 is a csv containing a subset of genes to construct Bayesian networks (Locus.ID)
network.nodes.df <- read.csv(file.choose())

#Expected input file 3 is a dictionary linking each probeset to a gene locus ID
#Expected format is a csv file with two fields: Probe.ID and Locus.ID
array.annotation <- read.csv(file.choose())
#array.annotation$Locus.ID <- toupper(array.annotation$Locus.ID) #optional step to ensure Locus.ID is uppercase


# Retrieve probeset for each network node locus ID
network.nodes.probesets <- array.annotation %>% 
  filter(Locus.ID %in% network.nodes.df$Locus.ID)

# Retrieve expression data for each network node
array.input <- normalised.input %>%
  filter(Probesets %in% network.nodes.probesets$Probe.ID)
array.input <- array.input[2:2374]
array.input.means <- rowMeans(array.input)

# Discretise expression data
discretised.output <- data.frame()
for (i in 1:nrow(array.input)) {
  for (j in 1:length(array.input[i,])) {
    if (array.input[i,j]-array.input.means[i] > 1) {
      discretised.output[i,j] <- 2
    } else if (array.input.means[i]-array.input[i,j] > 1){
      discretised.output[i,j] <- 0
    } else {
      discretised.output[i,j] <- 1
    }  
  }
}

discretised.matrix <- data.matrix(discretised.output)
discretised.matrix <- t(discretised.matrix)

# Write output matrix to use as BANJO input matrix
write.table(discretised.matrix, file = file.choose(), #select name and save as .txt file
            sep = '\t', eol = '\r\n', row.names = F, col.names = F)
