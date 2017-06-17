library(ggplot2)

setwd("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/supp_tevar/")

#When we moved to gbs data, we removed a few individuals from the selection analyses.  Now we have to additionally remove these individuals from this figure.

final77 <- c("RIMMA0388.1","RIMMA0389.1","RIMMA0390.1","RIMMA0392.1","RIMMA0393.1","RIMMA0394.2","RIMMA0395.2","RIMMA0396.1","RIMMA0397.1","RIMMA0398.1","RIMMA0399.1","RIMMA0403.2","RIMMA0404.1","RIMMA0406.1","RIMMA0409.1","RIMMA0416.1","RIMMA0417.1","RIMMA0418.1","RIMMA0421.1","RIMMA0422.1","RIMMA0423.1","RIMMA0424.1","RIMMA0425.1","RIMMA0426.1","RIMMA0428.1","RIMMA0430.1","RIMMA0431.1","RIMMA0433.1","RIMMA0436.1","RIMMA0437.1","RIMMA0438.1","RIMMA0441.1","RIMMA0462.1","RIMMA0464.1","RIMMA0465.1","RIMMA0466.1","RIMMA0467.1","RIMMA0468.1","RIMMA0473.1","RIMMA0614.1","RIMMA0615.1","RIMMA0616.1","RIMMA0619.1","RIMMA0620.1","RIMMA0621.1","RIMMA0623.1","RIMMA0625.1","RIMMA0628.1","RIMMA0630.1","RIMMA0657.1","RIMMA0658.1","RIMMA0661.1","RIMMA0662.1","RIMMA0663.1","RIMMA0667.1","RIMMA0671.1","RIMMA0674.1","RIMMA0680.1","RIMMA0690.1","RIMMA0691.1","RIMMA0696.1","RIMMA0700.1","RIMMA0701.1","RIMMA0702.1","RIMMA0703.1","RIMMA0708.1","RIMMA0709.1","RIMMA0710.1","RIMMA0712.1","RIMMA0716.1","RIMMA0720.1","RIMMA0721.1","RIMMA0727.1","RIMMA0729.1","RIMMA0730.1","RIMMA0731.1","RIMMA0733.1")
final2<-as.data.frame(final77)

png("~/Desktop/supp_tevar_77.png",width=900,height=500)

#RNA elements
listrna <- read.csv("ToptmpRetros2.txt", header=FALSE)
listlandr <- read.csv("Landrace_FTE_aggregated.csv",header=TRUE)
rna50 <- merge(listrna, listlandr, by.x="V1",by.y="FTE_group",sort=FALSE)
rna50$V1 <- NULL
trna50 <- as.data.frame(t(rna50))
gensz <- read.csv("RimmaGS.txt",header=TRUE,row.names="line")
gensz2 <- merge(gensz, final2, by.x="row.names",by.y="final77",sort=FALSE)
tmp <- merge(gensz2, trna50, by.x="Row.names",by.y="row.names",sort=FALSE)
tmp2 <- tmp[order(tmp$size, decreasing=TRUE), ]
sizes <- tmp2$size
tmp2$size <- NULL
accid <- tmp2$Row.names
tmp2$Row.names <- NULL
row.names(tmp2) <- accid
crap <- tmp2 * sizes
maxs <- apply(crap, 2, max)
mins <- apply(crap, 2, min)
scaled.rna15 <- scale(crap, center = mins, scale = maxs - mins )
par(mfrow=c(1,2))
image(scaled.rna15,xlab="Accessions by GS",yaxt='n',xaxt='n',main="Scaled Genomic MB of RNA Elements", col=topo.colors(50))
gencomp <- round(100*colMeans(crap), digits = 2)
umm <- seq(0,1, length.out = 15)
umm2 <- 0:1
axis(2, at=umm, labels=gencomp,las=2)
axis(1,at=umm2,labels=c("Largest","Smallest"))
text(x = -.05, y = 1.15, labels = "A", xpd = NA,cex = 2)

#DNA elements
listdna <- read.csv("ToptmpDNA2.txt", header=FALSE)
listlandr <- read.csv("Landrace_FTE_aggregated.csv",header=TRUE)
dna50 <- merge(listdna, listlandr, by.x="V1",by.y="FTE_group",sort=FALSE)
dna50$V1 <- NULL
tdna50 <- as.data.frame(t(dna50))
gensz <- read.csv("RimmaGS.txt",header=TRUE,row.names="line")
gensz2 <- merge(gensz, final2, by.x="row.names",by.y="final77",sort=FALSE)
tmp <- merge(gensz2, tdna50, by.x="Row.names",by.y="row.names",sort=FALSE)
tmp2 <- tmp[order(tmp$size, decreasing=TRUE), ]
sizes <- tmp2$size
tmp2$size <- NULL
accid <- tmp2$Row.names
tmp2$Row.names <- NULL
row.names(tmp2) <- accid
crap <- tmp2 * sizes
maxs <- apply(crap, 2, max)
mins <- apply(crap, 2, min)
scaled.dna15 <- scale(crap, center = mins, scale = maxs - mins )
image(scaled.dna15,xlab="Accessions by GS",yaxt='n',xaxt='n',main="Scaled Genomic MB of DNA Elements", col=topo.colors(50))
gencomp <- round(100*colMeans(crap), digits = 2)
umm <- seq(0,1, length.out = 15)
umm2 <- 0:1
axis(2, at=umm, labels=gencomp,las=2)
axis(1,at=umm2,labels=c("Largest","Smallest"))
text(x = -.05, y = 1.15, labels = "B", xpd = NA,cex = 2)
dev.off()
