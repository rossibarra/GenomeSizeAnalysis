#' ---
#' title: "Maize Genome Size Selection Analyses - Blank Rand GBS RESULTS"
#' author: "Paul Bilinski"
#' date: "May 11th, 2017"
#' ---

setwd("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/SelectionTests/Maize/")

options(scipen = 999)
library(data.table)
library(rrBLUP)
library(emma)
library(OpenMx)

source("../jri_emma_returnbeta.txt")

#GBS with 30% site coverage randomly imputed
#gbs<-read.table("~/Desktop/PaulCS/PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov.3_minTaxaCov.1.RndImp.Endelman.kinship.txt",skip=3,row.names = 1)
#mzfull <-as.matrix(gbs[1:77,1:77])
#A <- as.matrix(mzfull)
#iden <- as.data.frame(substr(row.names(mzfull),1,9))
#colnames(iden)[1] <- "identifier"

#GBS with 60% site coverage mean imputed via RRblup
#gbs<-read.table("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/SelectionTests/PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov0.6_minTaxaCov.1.UsMeanImpute.Endelman.RRBLUP.kinship.txt",skip=0,row.names = 1,header=TRUE)

#GBS with 60% coverage, KNN+mean
#data<-read.table("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/SelectionTests/PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov0.6_minTaxaCov.1.KNNImp.Endelman.kinship.txt",skip=3,row.names = 1)

#GBS with 60% coverage, KNN+rand
#data<-read.table("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/SelectionTests/PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov0.6_minTaxaCov.1.KNNImp.RndImp.Endelman.kinship",skip=3,row.names = 1)

#GBS with 60% coverage, Tassel Mean impute
#data<-read.table("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/SelectionTests/PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov0.6_minTaxaCov.1.MeanImpTassel.Endelman.kinship.txt",skip=3,row.names = 1)

#GBS with 60% site coverage RND imputed
data<-read.table("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/SelectionTests/PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov0.6_minTaxaCov.1.RndImp.Endelman.kinship",skip=3,row.names = 1)

#mzfull <-as.matrix(data[94:170,94:170])
#A <- as.matrix(mzfull)
#iden <- as.data.frame(substr(row.names(mzfull),1,9))
#colnames(iden)[1] <- "identifier"

mzfull <-as.matrix(data[1:77,1:77])
A <- as.matrix(mzfull)
iden <- as.data.frame(substr(row.names(mzfull),1,9))
colnames(iden)[1] <- "identifier"

#taxaID <- as.data.frame(substr(row.names(mzfull),1,9))
#write.csv(taxaID,"~/Desktop/maizetaxa.csv")

#CHECKING CHIPS
#chipgeno <- read.csv("~/Documents/Projects/Genome_Size_Analysis/Github_ParallelGS/SNP_data/Landrace_noSWUS_matrix.csv",header=TRUE,row.names=1)
#colnames(chipgeno)
#chipgeno[,c(4,15,17,19,36,53)] <- NULL
#iden <- as.data.frame(substr(colnames(chipgeno),1,9))
#colnames(iden)[1] <- "identifier"
#dt <- t(chipgeno)
#A <- A.mat(dt)
#A <- as.matrix(A)
#A.chip <- A
#cor(c(mzfull),c(A.chip))
#hundgbs<-read.csv("~/Desktop/Maize_gbs/convertedsnpsgbs.txt")
#hundgbs$X<-NULL
#dt <- t(hundgbs)
#A.2 <- A.mat(dt)
#cor(c(A),c(A.2))
#A <- A.2
#


pheno<- read.csv("Landrace_data.csv")
pheno$X <- NULL
pheno2 <- subset(pheno, pheno$Region!="SWUS")
pheno2$iden <- substr(pheno2$Row.names,1,9)

phenoorder <- merge(iden, pheno2, by.x="identifier", by.y="iden",sort=FALSE)

#Making mean and variance tables
kelly<- phenoorder[,c("identifier","Region","GS_bp","X180knobbp","TR1bp","TotallTebp")]
tapply(kelly$GS_bp,kelly$Region, summary)
tapply(kelly$GS_bp,kelly$Region, var)
tapply(kelly$X180knobbp,kelly$Region, summary)
tapply(kelly$X180knobbp,kelly$Region, var)
tapply(kelly$TR1bp,kelly$Region, summary)
tapply(kelly$TR1bp,kelly$Region, var)
tapply(kelly$TotallTebp,kelly$Region, summary)
tapply(kelly$TotallTebp,kelly$Region, var)
#

phenoorder$Region
mzSA <-as.matrix(data[c(1:14,24:27,29:31,33:39,50:55,59:60,66),c(1:14,24:27,29:31,33:39,50:55,59:60,66)])
A <- as.matrix(mzSA)
iden <- as.data.frame(substr(row.names(mzSA),1,9))
colnames(iden)[1] <- "identifier"

phenoorder <- merge(iden, pheno2, by.x="identifier", by.y="iden",sort=FALSE)
phenoorder$Region

mzM <-as.matrix(data[c(15:23,28,32,40:49,56:58,61:65,67:77),c(15:23,28,32,40:49,56:58,61:65,67:77)])
A <- as.matrix(mzM)
iden <- as.data.frame(substr(row.names(mzM),1,9))
colnames(iden)[1] <- "identifier"

phenoorder <- merge(iden, pheno2, by.x="identifier", by.y="iden",sort=FALSE)
phenoorder$Region


gs <- t ( as.matrix ( phenoorder$GS_bp, ncol = 1 ) )
gsgb <- gs/1000000000
alt <- t ( as.matrix ( phenoorder$Altitude , ncol = 1 ) )
alt <- alt - mean ( alt ) + 0.5

out.gsgb <- emma.REML.t_beta ( gsgb , alt , K = A )
out.gsgb

# pval: 0.00000006426576
# beta: -0.0001327578

knob <- t ( as.matrix ( phenoorder$X180knobbp, ncol = 1 ) )
knobgb <- knob/1000000000
out.knobgb <- emma.REML.t_beta ( knobgb , alt , K = A )

out.knobgb
# pval: 0.0000000442801
# beta: -0.00007538086

#selection on just tr1
tr1 <- t ( as.matrix ( phenoorder$TR1bp, ncol = 1 ) )
tr1gb <- tr1/1000000000
out.tr1gb <- emma.REML.t_beta ( tr1gb , alt , K = A )
out.tr1gb
# pval: 0.000000000490372
# beta: -0.00001248793

#selection on centc
centc <- t ( as.matrix ( phenoorder$CentCbp, ncol = 1 ) )
centcgb <- centc/1000000000
out.centcgb <- emma.REML.t_beta ( centcgb , alt , K = A )

out.centcgb
# pval: 0.08047364
# beta: 0.0000006048351

#selection on overall TE's
tes <- t ( as.matrix ( phenoorder$TotallTebp, ncol = 1 ) )
tesgb <- tes/1000000000 
out.tesgb <- emma.REML.t_beta ( tesgb , alt , K = A )
out.tesgb
#pval: 0.0000004269027
#beta: -0.00007723322

#selection on genome size without TR1 knobs
minustr1 <- t ( as.matrix ( phenoorder$GS_bp - phenoorder$TR1bp , ncol = 1 ) )
minustr1bp <- minustr1/1000000000
out.minustr1bp <- emma.REML.t_beta ( minustr1bp , alt , K = A )
out.minustr1bp

#selection on genome size without 180 and TR1 knobs
minusallknob <- t ( as.matrix ( phenoorder$GS_bp - phenoorder$TR1bp -phenoorder$X180knobbp , ncol = 1 ) )
minusallknob <- minusallknob/1000000000
out.minusallknob <- emma.REML.t_beta ( minusallknob , alt , K = A )
out.minusallknob


#Testing for selection on repeat elements when controlling for genome size
#The beta values are returned in units of GB/M, multiply by 1000 to convert to MB/m
fin.te <- emma.REML.t_beta ( tes , alt , X0 =  cbind ( rep ( 1 , 77 ) , c ( gsgb ) ) , K = A )
fin.te
# pval: 0.7914124, no selection on overall TE content when we account for GS
# beta: -1947.271

#knobgb 180
fin.knob <- emma.REML.t_beta ( knobgb , alt , X0 =  cbind ( rep ( 1 , 77 ) , c ( gsgb ) ) , K = A )
fin.knob
# pval: 0.03202736, no selection on 180knobs when we account for GS
# beta: -0.00002386128

#tr1gb
fin.tr1 <- emma.REML.t_beta ( tr1gb , alt , X0 =  cbind ( rep ( 1 , 77 ) , c ( gsgb ) ) , K = A )
fin.tr1
# pval: 0.0003238046, ?!? SELECTION ON on tr1 when we account for GS
#beta: -0.000006763093

#centcgb
fin.centc <- emma.REML.t_beta ( centcgb , alt , X0 =  cbind ( rep ( 1 , 77 ) , c ( gsgb ) ) , K = A )
fin.centc
#pval:0.1887889
#beta:0.0000005544113

##for mexico separately
#Testing for selection on repeat elements when controlling for genome size IN MEXICO
#The beta values are returned in units of GB/M, multiply by 1000 to convert to MB/m
fin.te <- emma.REML.t_beta ( tes , alt , X0 =  cbind ( rep ( 1 , 40 ) , c ( gsgb ) ) , K = A )
fin.te


#knobgb 180
fin.knob <- emma.REML.t_beta ( knobgb , alt , X0 =  cbind ( rep ( 1 , 40 ) , c ( gsgb ) ) , K = A )
fin.knob


#tr1gb
fin.tr1 <- emma.REML.t_beta ( tr1gb , alt , X0 =  cbind ( rep ( 1 , 40 ) , c ( gsgb ) ) , K = A )
fin.tr1

##FOR SOUTH AMERICA
#Testing for selection on repeat elements when controlling for genome size IN SOUTH AMERICA
#The beta values are returned in units of GB/M, multiply by 1000 to convert to MB/m
fin.te <- emma.REML.t_beta ( tes , alt , X0 =  cbind ( rep ( 1 , 37 ) , c ( gsgb ) ) , K = A )
fin.te


#knobgb 180
fin.knob <- emma.REML.t_beta ( knobgb , alt , X0 =  cbind ( rep ( 1 , 37 ) , c ( gsgb ) ) , K = A )
fin.knob


#tr1gb
fin.tr1 <- emma.REML.t_beta ( tr1gb , alt , X0 =  cbind ( rep ( 1 , 37 ) , c ( gsgb ) ) , K = A )
fin.tr1


#centcgb
fin.centc <- emma.REML.t_beta ( centcgb , alt , X0 =  cbind ( rep ( 1 , 37 ) , c ( gsgb ) ) , K = A )
fin.centc

#analysis of selection on all of the TE subfamilies.  Testing to see if more of them show selection than we would expect by chance.
perte <- phenoorder[,4:1190]
colnames(perte)

#run it with named var
try1 <- t ( as.matrix ( perte$DHH_Hip1_1 , ncol = 1 ) )
tmp <- emma.REML.t ( try1 , alt , X0 =  cbind ( rep ( 1 , 77 ) , c ( gsgb ) ) , K = A )
tmp
#now with column number var
tmp <- emma.REML.t ( perte[,1] , alt , X0 =  cbind ( rep ( 1 , 77 ) , c ( gsgb ) ) , K = A )
tmp
#they return the same
as.numeric(tmp$ps)

#NOW FOR DA WHOLE CABOODLE.
pval <- c()
for(i in 1:1187){
  tmp <- emma.REML.t ( perte[,i] , alt , X0 =  cbind ( rep ( 1 , 77 ) , c ( gsgb ) ) , K = A )
  pval <- c(pval , as.numeric(tmp$ps))
}
#So on the first run of this, it gives me an error where I have my all 0 values.  have to remove those 5 TEs.  I do so by seeing where the pval stops, and thats the next one that i need to NULL out.

perte <- phenoorder[,4:1190]
#we want to go to the original data frame and start nulling out those column IDs, not nulling out one at a time

perte2 <- perte[,-which(colSums(perte)==0)]

#So all of those numbers are a bit off from the original data matrix because im nulling out one at a time.  Below are the indexes in the full data frame maybe?

#also have to 0 out the pval before reruning
pval <- c()
for(i in 1:ncol(perte2)){
  tmp <- emma.REML.t ( perte2[,i] , alt , X0 =  cbind ( rep ( 1 , 77 ) , c ( gsgb ) ) , K = A )
  pval <- c(pval , as.numeric(tmp$ps))
}

length(which(pval< 0.05))

binom.test(46, 1156, p = 0.05)

0.05*ncol(perte2)



#rmarkdown::render("Maize_selectiontests.R", "pdf_document")

