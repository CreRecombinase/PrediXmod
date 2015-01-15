####by Heather E. Wheeler 20140602####
#Modifications by NWK 01/08/2015 ####
library(peer)
library(preprocessCore)
library(GenABEL)
args <- commandArgs(trailingOnly=T)
date <- Sys.Date()

if(length(args)!=3){
  print("ERROR!:")
  print(length(args))
  stop()
}

###############################################
### Directories & Variables
tissue <- args[1]
my.dir <- args[2]
project.name <- args[3]

Nk <- 15 ##number of peer factors to calculate, recommend 25% of sample size, but no more than 100, GTEx included 15 in pilot analyses
tis <- gsub(" ","",tissue)
 
#tissue <- "Nerve - Tibial" ###check GTEx_Analysis_2014-06-13.SampleTissue.annot for available tissues###

 
setwd(my.dir)

filenames <- paste0(project.name,c(".SampleTissue.annot",
                                   ".RNA-seq.ID.list",
                                   ".RNA-seq.GENE.list",
                                   ".protein-coding.RDS",
                                   ".fam",
                                   ".20genotPCs.txt",
                                   ".protein-coding.anno.RDS"))
if(!all(file.exists(filenames))){
  print("Unable to find files:")
  print(sum(!file.exists(filenames)))
  print(paste(filenames[!file.exists(filenames)]))
  stop()
}
names(filenames) <- c("SampleTissue","RNA-ID","RNA-Genelist","RNA-matrix","fam","PC","ANNO.RDS")
                    
print("Reading SampleTissue")
sam <- read.table(filenames["SampleTissue"],header=T,sep="\t",stringsAsFactors=F,quote='',comment.char='') 
sample <- subset(sam,SMTSD == tissue) ### pull sample list of chosen tissue###

expdata <- t(readRDS(filenames["RNA-matrix"]))
expanno <- readRDS(filenames["ANNO.RDS"])
rownames(expanno) <- expanno$id



expdata <- expdata[intersect(rownames(expdata),sample$SAMPID),] ###pull expression data for chosen tissue###
expdf <- data.frame(id=substr(rownames(expdata),1,10),expdata,stringsAsFactors=F)

print("Reading fam")
fam <- read.table(filenames["fam"],stringsAsFactors=F)
gtsamplelist <- fam$V1
gtsamplelist <- substr(gtsamplelist,1,10) ###to match with exp data###
rownames(fam) <- gtsamplelist
samplelist <- intersect(rownames(fam),expdf$id)
nsample <- length(samplelist)
expdf <- aggregate(.~id,data=expdf,FUN=mean)
rownames(expdf) <- expdf$id
expdata <- data.matrix(expdf[,-1])
expdata <- expdata[samplelist,]
expdata <- expdata[,apply(expdata,2,function(x)mean(x)>0)]

 

###get first 3 PCs from genos, as in GTEx
print("Reading PCs")
pc.matrix<-read.table(filenames["PC"],header=T,stringsAsFactors=F)
rownames(pc.matrix) <- substr(pc.matrix[,1],1,10)
pc.matrix <- data.matrix(pc.matrix[,3:5])
pc.matrix <- pc.matrix[samplelist,]


###pull gender, used as cov in GTEx
fam <- fam[samplelist,]
gender <- fam$V5
names(gender)<-rownames(fam)

###quantile normalize and transform to standard normal exp.w.geno matrix, as in GTEx###
t.expdata <- t(expdata)

rowtable<-function(x) length(table(x))>2 ##function to determine if more than 2 exp levels per gene
nonbin<-apply(t.expdata,1,rowtable) ##apply to matrix
t.expdata <- t.expdata[nonbin,] ##remove binary genes from matrix

qn.t.expdata <- normalize.quantiles(t.expdata) ##quantile normalize
rn.qn.t.expdata <- apply(qn.t.expdata,1,rntransform) ##rank transform to normality & transposes, not sure why?##

###Now we can create the model object, ### from https://github.com/PMBio/peer/wiki/Tutorial

model = PEER()

###set the observed data,

PEER_setPhenoMean(model,data.matrix(rn.qn.t.expdata))

dim(PEER_getPhenoMean(model))

###(NULL response means no error here), say we want to infer K=20 hidden confounders,

PEER_setNk(model,Nk)

PEER_getNk(model)

####and perform the inference. ###for Nk=20 and GTEx-NT, it took 323 iterations, for Nk=15 and GTEx-NT, it took 37 iterations

PEER_update(model)

factors = PEER_getX(model)
rownames(factors) <- rownames(expdata)
filename <- paste0(tis,".",Nk,".PEER.factors.",date,".txt")
write.table(factors,file =filename, quote=F)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)

##png(file= paste0(tis,".",Nk,".PEER.factors.plotmodel.",date,".png"))
##PEER_plotModel(model)
##dev.off()

adj.exp.matrix<-matrix(NA,nrow=nrow(rn.qn.t.expdata),ncol=ncol(rn.qn.t.expdata))

for(i in 1:ncol(rn.qn.t.expdata)){
	res <- summary(lm(rn.qn.t.expdata[,i] ~ factors + pc.matrix + gender, na.action=na.exclude))
	resid <- residuals(res)
	adj.exp.matrix[,i] <- resid
}

colnames(adj.exp.matrix) <- rownames(t.expdata)
rownames(adj.exp.matrix) <- colnames(t.expdata)

outfiles <- paste0(project.name,".",tis,".exp.adj.",Nk,c("PEERfactors.3PCs.gender.IDxGENE",
                                                         "PEERfactors.3PCs.gender.GENE.list",
                                                         "PEERfactors.3PCs.gender.ID.list",
                                                         "PEERfactors.3PCs.gender.ANNO."))
names(outfiles) <- c("IDxGENE",
                     "GENE.list",
                     "ID.list",
                     "ANNO.RDS")
write.table(adj.exp.matrix, file= outfiles["IDxGENE"], quote=F, row.names=F, col.names=F)
saveRDS(adj.exp.matrix,file=paste0(outfiles["IDxGENE"],".RDS"))

write(colnames(adj.exp.matrix), file = outfiles["GENE.list"], ncolumns=1)
write(rownames(adj.exp.matrix), file = outfiles["ID.list"], ncolumns=1)
