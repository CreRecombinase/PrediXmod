args <- commandArgs(trailingOnly=T)



project.name <- args[1]
my.dir <- args[2]
tissue <- args[3]
tis <- gsub(" ","",tissue)

Nk <- 15

setwd(my.dir)

files <- paste0( project.name,c(paste0(".",tis,".exp.adj.",Nk,"PEERfactors.3PCs.gender.IDxGENE.RDS"),
                                ".SNPxID.RDS",
                                ".SNPanno.RDS",
                                ".protein-coding.anno.RDS"))
names(files) <- c("EXPdata","SNPdata","SNPanno","EXPanno")
if(!all(file.exists(files)[-1])){
    print("File not found:")
    print(files[!file.exists(files)])
    stop()
}

print("Reading Annotation files!")
snpdata <- readRDS(files["SNPdata"])
snpanno <- readRDS(files["SNPanno"])
expanno <- readRDS(files["EXPanno"])

rownames(snpanno) <- snpanno$rsid
rownames(expanno) <- expanno$id

setwd(tis)

print(paste0("Reading ",tis," specific files"))

expdata <- readRDS(files["EXPdata"])




sizesave <- function(snpanno,geneanno,snpdata,expdata,i,j,s,project.name,tis){
    if(ncol(expdata)<300){
        outfiles <- paste0(project.name,".",tis,c(".SNPanno.",
                                                  ".EXPanno.",
                                                  ".IDxSNP.",
                                                  ".IDxGENE."),i,".",j,".RDS")
        names(outfiles) <- c("SNPANNO","EXPANNO","SNPDATA","EXPDATA")
        print(paste(i,j,sep="."))
        saveRDS(snpanno,outfiles["SNPANNO"])
        saveRDS(geneanno,outfiles["EXPANNO"])
        saveRDS(snpdata,outfiles["SNPDATA"])
        saveRDS(expdata,outfiles["EXPDATA"])
    }
    else{
        sizesave(snpanno,geneanno[1:300,],snpdata,expdata[,1:300],i,j,s,project.name,tis)
        s <- c(s,length(s)+2)
        sizesave(snpanno,geneanno[301:nrow(geneanno),],snpdata,expdata[,301:ncol(expdata)],i,s[length(s)],s,project.name,tis)
    }
}

rownames(expdata) <- sapply(strsplit(rownames(expdata),"-",fixed=T),function(x)paste0(x[1],"-",x[2]))
snpdata <- t(snpdata)
snpdata <- snpdata[rownames(snpdata) %in% rownames(expdata),]
expdata <- expdata[rownames(expdata) %in% rownames(snpdata),]

snpanno <- snpanno[rownames(snpanno) %in% colnames(snpdata),]
expanno <- expanno[rownames(expanno) %in% colnames(expdata),]


snpl <- split(snpanno,snpanno$chrom)
genel <- split(expanno,expanno$chrom)
rm(snpanno)
rm(expanno)
gc()
x <- list()
for(j in names(snpl)){
    print(j)
    tgenel <- genel[[j]]
    tsnpl <- snpl[[j]]
    x[[j]]<- numeric()
    for(i in seq(1e06,max(c(max(tgenel$genestop),max(tsnpl$pos))),1e06)){
        lessgenepos <- max(c(tgenel$genestop[tgenel$genestop<i],0),na.rm=T)
        lesssnppos <- max(c(tsnpl$pos[tsnpl$pos<i],0),na.rm=T)
        moregenepos <- min(c(tgenel$genestart[tgenel$genestart>i],max(tgenel$genestart)),na.rm = T)
        moresnppos <- min(c(tsnpl$pos[tsnpl$pos>i],max(tsnpl$pos)),na.rm=T)
        if(((i-lesssnppos)>1e06&&(i-lessgenepos)>1e06)&&((moregenepos-i)>1e06&&(moresnppos-i)>1e06)){
            x[[j]] <- c(x[[j]],i)
        }
    }
    gc()   
}

for(i in names(x)){
    print(i)
    tsnpl <- snpl[[i]]
    tgenel <- genel[[i]]
    s <- seq(length.out=length(x[[i]]))
    for(j in s){
        print(j)
        if(j==1){
            indexSNP <- tsnpl$pos<x[[i]][j]
            indexEXP <- tgenel$genestop<x[[i]][j]
        }
        else{
            indexSNP <- (tsnpl$pos<x[[i]][j])&(tsnpl$pos>x[[i]][j-1])
            indexEXP <- (tgenel$genestop<x[[i]][j])&(tgenel$genestart>x[[i]][j-1])
        }
        if((sum(indexSNP)>0)&&(sum(indexEXP)>0)){
            tsnpanno <- tsnpl[indexSNP,]
            tgeneanno <- tgenel[indexEXP,]
            tsnpdata <- snpdata[,rownames(tsnpanno),drop=F]
            texpdata <- expdata[,rownames(tgeneanno),drop=F]
            s <- sizesave(tsnpanno,tgeneanno,tsnpdata,texpdata,i,j,s,project.name,tis)
            
        }
        else{
            if(j!=length(x[[i]])){
                print(paste(i,j,"empty!"))
            }
        }
        if(j==length(x[[i]])){
            indexSNP <- tsnpl$pos>x[[i]][j]
            indexEXP <- tgenel$genestart>x[[i]][j]
            if((sum(indexSNP)>0)&&(sum(indexEXP)>0)){
                tsnpanno <- tsnpl[indexSNP,]
                tgeneanno <- tgenel[indexEXP,]
                tsnpdata <- snpdata[,rownames(tsnpanno),drop=F]
                texpdata <- expdata[,rownames(tgeneanno),drop=F]
                s <- sizesave(tsnpanno,tgeneanno,tsnpdata,texpdata,i,j+1,s,project.name,tis)
            }else{
                print(paste(i,j,"empty!"))
            }
        }        
    }
}
            
        
       





                                
    
