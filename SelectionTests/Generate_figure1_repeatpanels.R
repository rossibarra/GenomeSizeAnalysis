#' ---
#' title: "Making Figure 1 for the manuscript"
#' author: "Paul Bilinski"
#' date: "June 2017"
#' ---

library(ggplot2)
library(cowplot)
library(gridExtra)


#Kinship matrix to ensure that everything is accurately captured with regard to the genotyping data
#GBS with 60% site coverage RND imputed
data<-read.table("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/SelectionTests/PaulCellSize_ZeaGBSv27raw_Poly_minSiteCov0.6_minTaxaCov.1.RndImp.Endelman.kinship",skip=3,row.names = 1)

#Maize Pheno
mzfull <-as.matrix(data[1:77,1:77])
A <- as.matrix(mzfull)
iden <- as.data.frame(substr(row.names(mzfull),1,9))
colnames(iden)[1] <- "identifier"
pheno<- read.csv("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/SelectionTests/Maize/Landrace_data.csv")
pheno$X <- NULL
pheno2 <- subset(pheno, pheno$Region!="SWUS")
pheno2$iden <- substr(pheno2$Row.names,1,9)

maizephen <- merge(iden, pheno2, by.x="identifier", by.y="iden",sort=FALSE)
maizephen <- maizephen[,c(1,1198:1231)]
maizephen$GS_bp <- maizephen$GS_bp/1000000
maizephen$X180knobbp <- maizephen$X180knobbp/1000000
maizephen$TR1bp <- maizephen$TR1bp/1000000
maizephen$TotallTebp <- maizephen$TotallTebp/1000000
maizephen$Region<-substr(maizephen$Region,1,1)
maizephen$Region <- gsub("M", "MA", maizephen$Region)
maizephen$Region <- gsub("S", "SA", maizephen$Region)

cor(pheno2$TotallTe.,pheno2$GS_bp)
cor(pheno2$TR1.,pheno2$GS_bp)
cor(pheno2$X180knob.,pheno2$GS_bp)

cor(pheno2$TotallTebp,pheno2$GS_bp)
cor(pheno2$TR1bp,pheno2$GS_bp)
cor(pheno2$X180knobbp,pheno2$GS_bp)

#Mexicana Pheno
mexphen <- read.csv("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/SelectionTests/Mexicana/Pheno_alt_threshold_ordered.csv")
mexphen$X1C_GS <- mexphen$X1C_GS/2
mexphen$GS_bp <- mexphen$GS_bp/2000000
mexphen$X180knobbp <- mexphen$X180knobbp/2000000
mexphen$TR1bp <- mexphen$TR1bp/2000000
mexphen$TotalTebp <- mexphen$TotalTebp/2000000
mexphen$Species <- gsub("mexicana", "Mexicana", mexphen$Species)
plegend <- ggplot(mexphen, aes(Altitude, GS_bp,color=Species)) + geom_point()+ ylab("1C Genome Size (MB)") +geom_smooth(aes(group=Species),method="lm",color="red",linetype="dashed",se=FALSE)+xlab("Altitude (m)")+ theme(axis.text=element_text(size=SIZE))+scale_color_manual(name="Species",values="#E92F07")

SIZE=18

p1 <- ggplot(mexphen, aes(Altitude, GS_bp),) + geom_point(color="red")+ ylab("1C Genome Size (MB)") +geom_smooth(aes(group=Species),method="lm",color="red",linetype="dashed",se=FALSE)+guides(color=FALSE) +xlab("Altitude (m)")+ theme(axis.text=element_text(size=SIZE))
p2 <- ggplot(mexphen, aes(Altitude, TotalTebp),) + geom_point(color="red")+ ylab("TE Content (MB)") +geom_smooth(aes(group=Species),method="lm",color="red",linetype="dashed",se=FALSE)+guides(color=FALSE) +xlab("Altitude (m)")+ theme(axis.text=element_text(size=SIZE))
p3 <- ggplot(mexphen, aes(Altitude, X180knobbp),) + geom_point(color="red")+ ylab("180bp Knob (MB)") +geom_smooth(aes(group=Species),method="lm",color="red",linetype="dashed",se=FALSE)+guides(color=FALSE) +xlab("Altitude (m)")+ theme(axis.text=element_text(size=SIZE))
p4 <- ggplot(mexphen, aes(Altitude, TR1bp),) + geom_point(color="red")+ ylab("TR1 Knob (MB)") +geom_smooth(aes(group=Species),method="lm",color="red",linetype="dashed",se=FALSE)+guides(color=FALSE) +xlab("Altitude (m)")+ theme(axis.text=element_text(size=SIZE))

mlegend<-ggplot(maizephen, aes(Altitude, GS_bp,color=Region)) + geom_point()+ ylab("1C Genome Size (MB)") +geom_smooth(method="lm",linetype="dashed",se = FALSE)+xlab("Altitude (m)")+  scale_fill_manual(name="Maize Region",values=c("#CD853F", "#6ca6cd"), breaks=c("Mesoamerica", "South America"), labels=c("MX", "SA"))+scale_color_manual(values=c("#CD853F", "#6ca6cd")) + theme(axis.text=element_text(size=SIZE))

m1<-ggplot(maizephen, aes(Altitude, GS_bp,color=Region)) + geom_point()+ ylab("1C Genome Size (MB)") +geom_smooth(method="lm",linetype="dashed",se = FALSE)+xlab("Altitude (m)")+  scale_fill_manual(values=c("#CD853F", "#6ca6cd"), name="Region", breaks=c("MA", "SA"), labels=c("MX", "SA"))+scale_color_manual(values=c("#CD853F", "#6ca6cd")) + theme(axis.text=element_text(size=SIZE))+guides(color=FALSE)
m2<-ggplot(maizephen, aes(Altitude, TotallTebp,color=Region)) + geom_point()+ ylab("TE Content (MB)") +geom_smooth(method="lm",linetype="dashed",se = FALSE)+xlab("Altitude (m)")+  scale_fill_manual(values=c("#CD853F", "#6ca6cd"), name="Region", breaks=c("MA", "SA"), labels=c("MX", "SA"))+scale_color_manual(values=c("#CD853F", "#6ca6cd")) + theme(axis.text=element_text(size=SIZE))+guides(color=FALSE)
m3<-ggplot(maizephen, aes(Altitude, X180knobbp,color=Region)) + geom_point()+ ylab("180bp Knob (MB)") +geom_smooth(method="lm",linetype="dashed",se = FALSE)+xlab("Altitude (m)")+  scale_fill_manual(values=c("#CD853F", "#6ca6cd"), name="Region", breaks=c("MA", "SA"), labels=c("MX", "SA"))+scale_color_manual(values=c("#CD853F", "#6ca6cd")) + theme(axis.text=element_text(size=SIZE))+guides(color=FALSE)
m4<-ggplot(maizephen, aes(Altitude, TR1bp,color=Region)) + geom_point()+ ylab("TR1 Knob (MB)") +geom_smooth(method="lm",linetype="dashed",se = FALSE)+xlab("Altitude (m)")+  scale_fill_manual(values=c("#CD853F", "#6ca6cd"), name="Region", breaks=c("MA", "SA"), labels=c("MX", "SA"))+scale_color_manual(values=c("#CD853F", "#6ca6cd")) + theme(axis.text=element_text(size=SIZE))+guides(color=FALSE)

legend <- get_legend(mlegend)
tmp <- plot_grid(m1,m2,m3,m4,labels=c("A", "B","C","D"), ncol = 4, nrow = 1)

p <- plot_grid( tmp, legend, rel_widths = c(3, .3))

legend2 <- get_legend(plegend)

tmp2 <- plot_grid(p1,p2,p3,p4, labels=c("E", "F","G","H"), ncol = 4, nrow = 1)

p2 <- plot_grid( tmp2, legend2, rel_widths = c(3, .3))

#png("~/Desktop/Figure1_repeatpanel.png",width=1500,height=600)
#plot_grid(p,p2,ncol = 1, nrow = 2)
#dev.off()

maizeMpheno <- subset(maizephen,maizephen$Region=="MA")
maizeSpheno <- subset(maizephen,maizephen$Region=="SA")
cor.test(maizeMpheno$Altitude,maizeMpheno$GS_bp)
cor.test(maizeSpheno$Altitude,maizeSpheno$GS_bp)

cor.test(maizeMpheno$Altitude,maizeMpheno$TotallTebp)
cor.test(maizeSpheno$Altitude,maizeSpheno$TotallTebp)

cor.test(maizeMpheno$Altitude,maizeMpheno$X180knobbp)
cor.test(maizeSpheno$Altitude,maizeSpheno$X180knobbp)

cor.test(maizeMpheno$Altitude,maizeMpheno$TR1bp)
cor.test(maizeSpheno$Altitude,maizeSpheno$TR1bp)
