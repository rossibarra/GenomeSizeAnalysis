#' ---
#' title: "Mex Genome Size Selection Analyses - Blank RandImp GBS RESULTS"
#' author: "Paul Bilinski"
#' date: "May 11th, 2017"
#' ---

setwd("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/SelectionTests/Mexicana/")

library(data.table)
library(rrBLUP)
library(emma)
library(OpenMx)

source("../jri_emma_returnbeta.txt")

#data<-read.table("~/Desktop/PaulCS/PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov.3_minTaxaCov.1.RndImp.Endelman.kinship.txt",skip=3,row.names = 1)

#mexfull <-data[78:170,78:170]
#row.names(mexfull)

#remove AM, M, and TZ to get the high elevation material
#mexhigh <- mexfull[c(6:41,51:84),c(6:41,51:84)]
#A <- as.matrix(mexhigh)

#GBS with 60% site coverage mean imputed via RRblup
#data<-read.table("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/SelectionTests/PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov0.6_minTaxaCov.1.UsMeanImpute.Endelman.RRBLUP.kinship.txt",skip=0,row.names = 1,header=TRUE)

#GBS with 60% coverage, KNN+mean
#data<-read.table("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/SelectionTests/PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov0.6_minTaxaCov.1.KNNImp.Endelman.kinship.txt",skip=3,row.names = 1)

#GBS with 60% coverage, KNN+rand
#data<-read.table("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/SelectionTests/PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov0.6_minTaxaCov.1.KNNImp.RndImp.Endelman.kinship",skip=3,row.names = 1)

#GBS with 60% coverage, Tassel Mean impute
#data<-read.table("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/SelectionTests/PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov0.6_minTaxaCov.1.MeanImpTassel.Endelman.kinship.txt",skip=3,row.names = 1)

#GBS with 60% site coverage RND imputed
data<-read.table("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/SelectionTests/PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov0.6_minTaxaCov.1.RndImp.Endelman.kinship",skip=3,row.names = 1)

mexfull <-data[78:170,78:170]
row.names(mexfull)

#taxaID <- as.data.frame(substr(row.names(mexfull),1,3))
#write.csv(taxaID,"~/Desktop/teotaxa.csv")

#remove AM, M, and TZ to get the high elevation material
mexhigh <- mexfull[c(6:41,51:84),c(6:41,51:84)]
A <- as.matrix(mexhigh)

#the five variations on calculating kinship matrices
#PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov0.6_minTaxaCov.1.KNNImp.Endelman.kinship.txt
#PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov0.6_minTaxaCov.1.KNNImp.RndImp.Endelman.kinship
#PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov0.6_minTaxaCov.1.MeanImpTassel.Endelman.kinship.txt
#PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov0.6_minTaxaCov.1.RndImp.Endelman.kinship.txt
#PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov0.6_minTaxaCov.1.UsMeanImpute.Endelman.RRBLUP.kinship.txt

#mexfull <-data[1:93,1:93]
#row.names(mexfull)

#mexhigh <- mexfull[c(6:41,51:84),c(6:41,51:84)]
#A <- as.matrix(mexhigh)

#testing for F statistics variation
popam<- mexfull[c(1:5),c(1:5)]
popm<- mexfull[c(41:49),c(41:49)]
poptz<- mexfull[c(85:93),c(85:93)]
poplows<-mexfull[c(1:5,41:49,85:93),c(1:5,41:49,85:93)]

sum(diag2vec(popam))/length(popam)
sum(diag2vec(popm))/length(popm)
sum(diag2vec(poptz))/length(poptz)
sum(diag2vec(poplows))/length(poplows)
sum(diag2vec(mexfull))/length(mexfull)

t.test(diag2vec(popam),diag2vec(mexfull))
t.test(diag2vec(popm),diag2vec(mexfull))
t.test(diag2vec(poptz),diag2vec(mexfull))
t.test(diag2vec(poplows),diag2vec(mexfull))


pheno <- read.csv("Pheno_alt_threshold_ordered.csv") 

gs <- t ( as.matrix (pheno$GS_bp, ncol=1) )
gsgb <- gs/2000000000

alt <- t ( as.matrix (pheno$Altitude, ncol=1))
alt <- alt - mean(alt) + 0.5

knob180 <- t ( as.matrix (pheno$X180knobbp, ncol=1))
knob180gb <- knob180/2000000000

TR1 <- t ( as.matrix (pheno$TR1bp, ncol=1))
TR1gb <- TR1/2000000000

Te <- t ( as.matrix (pheno$TotalTebp, ncol=1))
Tegb <- Te/2000000000

minustr1 <- t ( as.matrix ( pheno$GS_bp - pheno$TR1bp , ncol = 1 ) )
minustr1bp <- minustr1/2000000000

minusallknob <- t ( as.matrix ( pheno$GS_bp - pheno$TR1bp -pheno$X180knobbp , ncol = 1 ) )
minusallknob <- minusallknob/2000000000

#Testing selection
out.gsgb <- emma.REML.t_beta( gsgb,alt,K=A)
out.gsgb

#pval: 5.456587e-05
#beta: -0.0002716386

out.knob <- emma.REML.t_beta(knob180gb,alt,K=A)
out.knob

#pval: 0.533834
#beta: -4.312718e-05

out.tr1 <- emma.REML.t_beta(TR1gb,alt,K=A)
out.tr1

#pval: 0.008343883
#beta: -7.800644e-06

out.te <- emma.REML.t_beta(Tegb,alt,K=A)
out.te

#pval: 0.0343131
#beta: -0.0001377701

out.minustr1bp <- emma.REML.t_beta ( minustr1bp , alt , K = A )
out.minustr1bp

out.minusallknob <- emma.REML.t_beta ( minusallknob , alt , K = A )
out.minusallknob


#Testing selection on components
fin.gs <- emma.REML.t_beta ( gsgb , alt , X0 =  cbind ( rep ( 1 , 70 ) , c ( TR1gb ) ) , K = A )
fin.gs

#pval: 4.861976e-05
#beta: -0.0002864481

fin.knob <- emma.REML.t_beta ( knob180gb , alt , X0 =  cbind ( rep ( 1 , 70 ) , c ( gsgb ) ) , K = A )
fin.knob

#pval:0.9705134
#beta: -1.510216e-05

fin.tr1 <- emma.REML.t_beta ( TR1gb , alt , X0 =  cbind ( rep ( 1 , 70 ) , c ( gsgb ) ) , K = A )
fin.tr1

#pval: 0.007979344
#beta: -9.098991e-06

fin.te <- emma.REML.t_beta ( Tegb , alt , X0 =  cbind ( rep ( 1 , 70 ) , c ( gsgb ) ) , K = A )
fin.te

#pval: 0.1843356
#beta: 7.281324e-05

#rmarkdown::render("Mexicana_selectiontests.R", "pdf_document")

