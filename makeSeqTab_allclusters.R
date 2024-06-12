#!/bin/Rscript
  
library(tidyverse)
library(Biostrings)
library(dplyr)

#d <- file.path("/path_to_sequences")
#primer1 <- 18S
#primer2 <- 12S
#mincluster <- 100
#inputprimer1 <- ""$primer1"_mapped_polished.fasta"
#inputprimer2 <- ""$primer2"_mapped_polished.fasta"

#outputprimer1 <- ""$primer1"_mapped_polished_min"$mincluster""
#outputprimer2 <- ""$primer2"_mapped_polished_min"$mincluster""


setwd("/path_to_sequences/cdhit")
fasta18S <- readDNAStringSet("18S_mapped_polished.fasta")
fasta12S <- readDNAStringSet("12S_mapped_polished.fasta")
fastaCOI <- readDNAStringSet("COI_mapped_polished.fasta")
fastaBelow <- readDNAStringSet("clustersbelowmin10threshold.fasta")

length(fasta18S)
length(fasta12S)
length(fastaCOI)
length(fastaBelow)

#SeqTable <- c(fasta18S,fasta12S)
SeqTable <- c(fasta18S,fastaCOI,fasta12S,fastaBelow) 
length(SeqTable)

SeqTablenames <- names(SeqTable)
SeqTableseqs <- paste(SeqTable)
SeqTable.df <- data.frame(SeqTablenames,SeqTableseqs)

new.SeqTable <- SeqTable.df %>% mutate(barcode = (word(SeqTablenames, 1,sep = ";")), .before = 1) %>%
mutate(SeqID = (word(SeqTablenames, 2, sep = ";")), .before = 2) %>%
mutate(size = (word(SeqTablenames, 3, sep = ";")), .before = 3)

mv.SeqTable <- new.SeqTable %>% mutate(Barcode = gsub("barcodelabel=","",barcode)) %>%
mutate(UniqueID = gsub("otu=","",SeqID)) %>%
mutate(Abundance = as.numeric(gsub("size=","",size)))
mv.SeqTable$barcode <- NULL
mv.SeqTable$SeqID <- NULL
mv.SeqTable$size <- NULL
mv.SeqTable$SeqTablenames <- NULL
mv.SeqTable$SeqTableseqs <- NULL

tidy.SeqTable <- mv.SeqTable %>%
pivot_wider(names_from = Barcode, values_from = Abundance, values_fill = 0)

write.csv(tidy.SeqTable,"/path_to_sequences/cdhit/18S_12S_COI_allclusters_secondcluster_seqtab.csv")

quit()
--no-save

