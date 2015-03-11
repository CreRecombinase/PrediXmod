args <- commandArgs(trailingOnly=T)

chromdir <- args[1]
chrompref <- args[2]
dbfile <- args[3]
idnum <- as.integer(args[4])
chunksize <- as.integer(args[5])
outfile <- args[6]

if(length(args)!=6){
    print(paste0("Usage chromdir chrompref dbfile idnum chunksize outfile"))
    stop()
}
allfiles <- dir(chromdir,pattern=chrompref,full.names = T)
refallele <- 4
library(RSQLite)
library(BBmisc)
library(plyr)
library(Matrix)

sqlite    <- dbDriver("SQLite")
con <- dbConnect(sqlite,dbfile)
eqtldf <- dbReadTable(con,"weights")
snplist <- unique(eqtldf$rsid)
genelist <- unique(eqtldf$gene)
genemat <- matrix(data = 0,nrow = length(genelist),ncol = 423)

freadsnps <- function(file,lines,snplist){
  con <- file(file)
  open(con)
  rlines <- 0
  fline <- readLines(con,n=1,warn=F)
  retmatrix <- matrix(unlist(strsplit(fline,split = " ",fixed = T)),nrow = 1,byrow = T)
  while(length(linelist <- readLines(con,n = lines,warn=F))>0){
    tmat <- matrix(unlist(strsplit(linelist,split = " ",fixed=T)),nrow = length(linelist),byrow=T)
    tmat <- tmat[tmat[,2] %in% snplist,]
    if(nrow(tmat)>0){
      retmatrix <- rbind(retmatrix,tmat)
    }
    rlines <- rlines+length(linelist)
    cat(paste0(rlines,"\n"))
  }
  rownames(retmatrix) <- retmatrix[,2]
  return(retmatrix)
}
for(i in 1:length(allfiles)){
  print(paste0("Reading in file ",allfiles[i]))
  dosagemat <- freadsnps(allfiles[i],100000,snplist)
  tempeqtldf <- eqtldf[ eqtldf$rsid %in% rownames(dosagemat),]
  tempeqtldf$dallele <- dosagemat[tempeqtldf$rsid,refallele]
  tempeqtldf$weight[ tempeqtldf$eff_allele != tempeqtldf$dallele] <- -tempeqtldf$weight[ tempeqtldf$eff_allele != tempeqtldf$dallele]
  eqtlmat <- sparseMatrix(i = match(tempeqtldf$gene,genelist),j=match(tempeqtldf$rsid,rownames(dosagemat)),x = tempeqtldf$weight,dims = c(length(genelist),nrow(dosagemat)))
  dmat <- dosagemat[,7:ncol(dosagemat)]
  mode(dmat) <- "numeric"
  genemat <- genemat +(eqtlmat%*%dmat)
}
genemat <- as.matrix(genemat)
genemat <- as.data.frame(genemat)
rownames(genemat) <- genelist
write.table(genemat,file = outfile,col.names=F,row.names=T,quote=F,sep="\t")
