library(rjags)
library(stringr)
library(mcmcplots)

#Master Data frame with all measures is called fgrote, for final grote dataframe.
fgrote <- read.csv("~/Documents/Projects/Genome_Size_Analysis/Github_ParallelGS/Mexicana/Finalgrowthchamberanalysis/grote_Indexmoms_groundedLL_d16NA.csv")

#Make the necessary dataframe inputs into the model.
stoma.positions <- str_detect(names(fgrote),"stom*")
stoma <- fgrote[,names(fgrote)[stoma.positions]]
#Convert microns to centimeters.
stoma <- 0.0001 * stoma 
leaflength.positions <- str_detect(names(fgrote), "leaf3diff*")#change this to leaf3diff*
leaflength <- fgrote[,names(fgrote)[leaflength.positions]]
n.cells <- apply(!is.na(stoma), MAR=1, FUN="sum")

#Data frames needed for figure generation.
gensize.conforming <- matrix(fgrote$genomesize,nrow=nrow(fgrote),ncol=3,byrow=FALSE)
fgrote$leaf3diff10 <- fgrote$leaf3.1 - fgrote$leaf3.0
fgrote$leaf3diff21 <- fgrote$leaf3.2 - fgrote$leaf3.1
fgrote$leaf3diff32 <- fgrote$leaf3.3 - fgrote$leaf3.2
#Empirical derivatives, change per day, necessary for figures.
LEdiffs <- as.matrix(cbind(fgrote$leaf3diff10,fgrote$leaf3diff21,fgrote$leaf3diff32))

#Call JAGS to generate posterior samples. 
DAGLL.jags <- jags.model("DAG_twomediator_Final.txt", 
	list(leaflength=log(leaflength), stoma=log(stoma), 
		n.mothers=max(fgrote$momindex), n.indiv=dim(fgrote)[1],
		n.cells=n.cells, n.times=3, mother.index=fgrote$momindex,
		GS=log(fgrote$genomesize)),
	n.chains=2)
update(DAGLL.jags, n.iter=1000)
DAGLL.samples <- coda.samples(DAGLL.jags, 
	c("beta.GS","tao.GS","GS.contrast",
		"sd.motherStCS","sd.motherLL","sd.indivCS","sd.indivLL","sd.StCS","sd.LL",
		"StCS.coeffs","tao.LL","LE.coeffs"),
	n.iter=5000, thin=50)
#mcmcplot(DAGLL.samples)
chain1 <- as.matrix(DAGLL.samples[[1]])
chain2 <- as.matrix(DAGLL.samples[[2]])
samples.mat <- rbind(chain1, chain2)
write.csv(chain1,"DAGLLfinal.jags.le.6_003_4.8.chain1.csv")
write.csv(chain2,"DAGLLfinal.jags.le.6_003_4.8.chain2.csv")

--------------
hist(table(fgrote$momindex),xlab="Size of Maternal Half-Sib Family",breaks=20,col="peru",ylab="Number of families",main="Histogram of the Size of Maternal Half-Sib Families",cex.main=0.9)
hist(as.factor(fgrote$momindex))
