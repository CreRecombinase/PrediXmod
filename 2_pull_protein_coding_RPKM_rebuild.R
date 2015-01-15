##Rewrite of 2_pull_protein_coding_RPKM_rebuild.pl in R

args <- commandArgs(trailingOnly=T)

if(length(args)!=4){
  print("Usage: ProjectName expfile genelist outdir")
  stop()
}

projectname <- args[1]
expfile <- args[2]
genelistfile <- args[3]
outdir <- args[4]
setwd(outdir)
outfiles <- paste0(projectname,c(".protein-coding.gct",
                                            ".RNA-seq.ID.list",
                                            ".RNA-seq.GENE.list",
                                            ".RNA-seq.GENExID",
                                            ".RNA-seq.GENEname.list",
                                            ".protein-coding.gct.gz",
                                            ".protein-coding.RDS",
                                            ".protein-coding.anno.RDS"))
names(outfiles) <- c("RPKM","ID","GENE","GENExID","GENENAME","RPKM.gz","EXPRDS","GENEANNO.RDS")


genelist <- read.table(genelistfile,sep="\t",header=F,stringsAsFactors=F)
genelist <- genelist[genelist[,3]=="transcript",]
genename <- sapply(strsplit(sapply(strsplit(genelist[,9],split=";"),"[",5)," "),"[",3)
geneid <- sapply(strsplit(sapply(strsplit(genelist[,9],split=";"),"[",1)," "),"[",2)
proteinCoding <- sapply(strsplit(sapply(strsplit(genelist[,9],split=";"),"[",3)," "),"[",3)
genanno <- data.frame(id=geneid,chrom=genelist[,1],genestart=genelist[,4],genestop=genelist[,4],genename=genename,stringsAsFactors=F)

genanno <- genanno[proteinCodin=="protein_coding",]
geneid <- geneid[proteinCoding=="protein_coding"]

saveRDS(genanno,outfiles["GENEANNO.RDS"])



expfile <- gzfile(expfile)
open(expfile)

firstexp <- readLines(expfile,n=3)
nrows <- as.numeric(strsplit(firstexp[2],"\t")[[1]][1])
ncols <- as.numeric(strsplit(firstexp[2],"\t")[[1]][2])


headerrow <- strsplit(firstexp[3],"\t")[[1]]

write(headerrow[-c(1:2)],file=outfiles["ID"],sep="\n")
write(headerrow,file=outfiles["RPKM"],sep="\t",ncolumns=length(headerrow))

explines <- readLines(expfile,n=nrows)
alldata <- strsplit(explines,"\t")
expdata <- t(sapply(alldata,function(x)as.numeric(x[-c(1:2)])))
colnames(expdata) <- headerrow[-c(1:2)]
rownames(expdata) <- sapply(alldata,"[",1)
expdata <- expdata[rownames(expdata) %in% geneid,]
write(rownames(expdata),file=outfiles["GENE"],sep="\n")
write(sapply(alldata,"[",2),file=outfiles["GENENAME"],sep="\n")
saveRDS(expdata,file=outfiles["EXPRDS"])




