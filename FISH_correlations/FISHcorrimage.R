data <- read.csv("~/Documents/Projects/Genome_Size_Analysis/GenomeSizeZea_Analyses/FISH_correlations/Birchler_FISHcor.csv")

data <- data[1:20,]


library(ggplot2)

#install.packages("cowplot")
library(cowplot)
library(gridExtra)

fishy<-ggplot(data, aes(knob180, X180knobmb)) + geom_point() + ylab("Measured Knob Content (MB)") + xlab("Observed Knob Count")+ geom_smooth(method="lm", se = TRUE, lty=2)+theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))

#png(filename="Fishcor.png",width=480,height=480,units="px")
#fishy
#dev.off()

fishy2<-ggplot(data, aes(tr1, tr1mb)) + geom_point() + ylab("Measured TR1  Content (MB)") + xlab("Observed TR1  Count")+ geom_smooth(method="lm", se = TRUE, lty=2)+theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))

#png(filename="supp_fishtr1count.png",width=480,height=480,units="px")
#fishy2
#dev.off()

grid.arrange(fishy,fishy2,ncol=2)
plot_grid(fishy,fishy2,ncol=2, labels=c("A","B"))

tmp <- plot_grid(m1,m2,m3,m4,labels=c("A", "B","C","D"), ncol = 4, nrow = 1)

cor.test(data$knob180,data$X180knobmb)
cor(data$tr1,data$tr1mb)
