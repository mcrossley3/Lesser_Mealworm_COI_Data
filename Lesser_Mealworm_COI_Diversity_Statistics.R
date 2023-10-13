setwd("[path/to/data]/Alignments")
library(pegas)

##Overall Alignment
#Read in Data
alph = read.dna('Alphitobius_Alignment.fasta', format="fasta",as.matrix=TRUE)
#Nucleotide diversity
All.pi = nuc.div(alph,variance = FALSE)
All.pi
#Haplotype diversity
All.H = hap.div(alph,variance = FALSE)
All.H

###USA

##Delaware Farm 1
DE_F1_H1 = read.dna('DE_Farm_1_House_1.fasta', format="fasta",as.matrix=TRUE)
DE_F1_H2 = read.dna('DE_Farm_1_House_2.fasta', format="fasta",as.matrix=TRUE)
DE_F1_H3 = read.dna('DE_Farm_1_House_3.fasta', format="fasta",as.matrix=TRUE)
DE_F1_H4 = read.dna('DE_Farm_1_House_4.fasta', format="fasta",as.matrix=TRUE)
#Nucleotide diversity
DE_F1_H1.pi = nuc.div(DE_F1_H1,variance = FALSE)
DE_F1_H2.pi = nuc.div(DE_F1_H2,variance = FALSE)
DE_F1_H3.pi = nuc.div(DE_F1_H3,variance = FALSE)
DE_F1_H4.pi = nuc.div(DE_F1_H4,variance = FALSE)
#Haplotype diversity
DE_F1_H1.H = hap.div(DE_F1_H1,variance = FALSE)
DE_F1_H2.H = hap.div(DE_F1_H2,variance = FALSE)
DE_F1_H3.H = hap.div(DE_F1_H3,variance = FALSE)
DE_F1_H4.H = hap.div(DE_F1_H4,variance = FALSE)

##Delaware Farm 2
DE_F2_H1 = read.dna('DE_Farm_2_House_1.fasta', format="fasta",as.matrix=TRUE)
DE_F2_H2 = read.dna('DE_Farm_2_House_2.fasta', format="fasta",as.matrix=TRUE)
#Nucleotide diversity
DE_F2_H1.pi = nuc.div(DE_F2_H1,variance = FALSE)
DE_F2_H2.pi = nuc.div(DE_F2_H2,variance = FALSE)
#Haplotype diversity
DE_F2_H1.H = hap.div(DE_F2_H1,variance = FALSE)
DE_F2_H2.H = hap.div(DE_F2_H2,variance = FALSE)

###Alabama
AL_F1_H1 = read.dna('AL_Farm_1_House_1.fasta', format="fasta",as.matrix=TRUE)
AL_F1_H2 = read.dna('AL_Farm_1_House_2.fasta', format="fasta",as.matrix=TRUE)
AL_F2 = read.dna('AL_Farm_2.fasta', format="fasta",as.matrix=TRUE)
AL_F3 = read.dna('AL_Farm_3.fasta', format="fasta",as.matrix=TRUE)
#Nucleotide diversity
AL_F1_H1.pi = nuc.div(Alabama1,variance = FALSE)
AL_F1_H2.pi = nuc.div(Alabama2,variance = FALSE)
AL_F2.pi = nuc.div(Alabama3,variance = FALSE)
AL_F3.pi = nuc.div(Alabama4,variance = FALSE)
#Haplotype diversity
AL_F1_H1.H = hap.div(Alabama1,variance = FALSE)
AL_F1_H2.H = hap.div(Alabama2,variance = FALSE)
AL_F2.H = hap.div(Alabama3,variance = FALSE)
AL_F3.H = hap.div(Alabama4,variance = FALSE)

##Georgia
Geo = read.dna('Georgia.fasta',format='fasta',as.matrix=TRUE)
#Nucleotide diversity
Geo.pi = nuc.div(Geo,variance = FALSE)
#Haplotype diversity
Geo.H = hap.div(Geo,variance = FALSE)

###Maryland

##Maryland Farm 1
MD_F1 = read.dna('MD_Farm_1.fasta',format='fasta',as.matrix=TRUE)
#Nucleotide diversity
MD_F1.pi = nuc.div(MD_F1,variance = FALSE)
#Haplotype diversity
MD_F1.H = hap.div(MD_F1,variance = FALSE)

##Maryland Farm 2
MD_F2 = read.dna('MD_Farm_2.fasta',format='fasta',as.matrix=TRUE)
#Nucleotide diversity
MD_F2.pi = nuc.div(MD_F2,variance = FALSE)
#Haplotype diversity
MD_F2.H = hap.div(MD_F2,variance = FALSE)

###Europe

##Czech
Czech = read.dna('Czech.fasta',format='fasta',as.matrix=TRUE)
#Nucleotide diversity
Cze.pi = nuc.div(Czech,variance = FALSE)
#Haplotype diversity
Cze.H = hap.div(Czech,variance = FALSE)

##Italy
Italy = read.dna('Italy.fasta',format='fasta',as.matrix=TRUE)
#Nucleotide diversity
Ita.pi = nuc.div(Italy,variance = FALSE)
#Haplotype diversity
Ita.H = hap.div(Italy,variance = FALSE)

##Turkey
Turkey = read.dna('Turkey.fasta',format='fasta',as.matrix=TRUE)
#Nucleotide diversity
Tur.pi = nuc.div(Turkey,variance = FALSE)
#Haplotype diversity
Tur.H = hap.div(Turkey,variance = FALSE)

#Consolidate
Site=c('Overall','DE Farm 1 House 1','DE Farm 1 House 2','DE Farm 1 House 3','DE Farm 1 House 4','DE Farm 2 House 1','DE Farm 2 House 2','AL Farm 1 House 1','AL Farm 1 House 2','AL Farm 2',
        'AL Farm 3','Georgia','MD Farm 1','MD Farm 2','Czech','Italy','Turkey')
Pi=c(All.pi,DE_F1_H1.pi,DE_F1_H2.pi,DE_F1_H3.pi,DE_F1_H4.pi,DE_F2_H1.pi,DE_F2_H2.pi,AL_F1_H1.pi,AL_F1_H2.pi,AL_F2.pi,AL_F3.pi,Geo.pi,MD_F1.pi,MD_F2.pi,Cze.pi,Ita.pi,Tur.pi)
H=c(All.H,DE_F1_H1.H,DE_F1_H2.H,DE_F1_H3.H,DE_F1_H4.H,DE_F2_H1.H,DE_F2_H2.H,AL_F1_H1.H,AL_F1_H2.H,AL_F2.H,AL_F3.H,Geo.H,MD_F1.H,MD_F2.H,Cze.H,Ita.H,Tur.H)

cbind(Site,Pi,H)
