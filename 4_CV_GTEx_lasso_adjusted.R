####script file found in /nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/
####by Heather E. Wheeler 20140602####
##Modified by NWK ##
args <- commandArgs(trailingOnly=T)
###############################################
### Directories & Variables

library(glmnet)
#library(doMC)
library(plyr)
#registerDoMC(cores=2)
start.time <- Sys.time()
project.name <- args[1]
my.dir <- args[2]
tissue <- args[3]
ch <- args[4]
seg <- args[5]


tis <- gsub(" ","",tissue)


setwd(my.dir)
setwd(tis)
filenames <- paste0(project.name,".",tis,c(".SNPanno.",
                                           ".EXPanno.",
                                           ".IDxSNP.",
                                           ".IDxGENE."),ch,".",seg,".","RDS")
if(!all(file.exists(filenames))){
    print("Files not found:")
    print(filenames[!file.exists(filenames)])
    stop()
}
names(filenames) <- c("SNPANNO","EXPANNO","IDxSNP","IDxGENE")

################################################
###input adjusted expression data###
expdata <- readRDS(filenames["IDxGENE"])
snpdata <- readRDS(filenames["IDxSNP"])
snpanno <- readRDS(filenames["SNPANNO"])
expanno <- readRDS(filenames["EXPANNO"])


k <- 10 ### k-fold CV
n <- 3 #number of k-fold CV replicates
 
################################################
### Functions & Libraries

SNP_Exp <- function(gene,snpdata=snpdata,expdata=expdata,snpanno=snpanno,expanno=expanno){
    start <- max(expanno[gene,"genestart"]-1e6,0)
    end <- expanno[gene,"genestop"]+1e6
    cissnpdata <- snpdata[,snpanno$rsid[(snpanno$pos>=start & snpanno$pos<=end)],drop=F]
    cissnpdata[cissnpdata>=3] <- NA
    cissnpdata <- cissnpdata[,apply(cissnpdata,2,function(x)mean(x,na.rm=T)>0),drop=F]
    cissnpdata <- scale(cissnpdata,center=T,scale=T)
    cissnpdata[is.na(cissnpdata)] <- 0
    exppheno <- expdata[,gene]
    exppheno <- scale(exppheno,center=T,scale=T)
    exppheno[is.na(exppheno)] <- 0
    rownames(exppheno) <- rownames(expdata)
    genename <- expanno[gene,"genename"]
    return(list(exp=exppheno,snp=cissnpdata,genename=genename,ensid=gene))
}

glmnet_wrapper <- function(listel,alpha,nrep.set,nfold.set,isParallel=F){
    exp <- listel[["exp"]]
    snp <- listel[["snp"]]
    genename <- listel[["genename"]]
    if(ncol(snp)<=2){
        print(paste0("No SNPs for gene: ",genename))
        return(NULL)
    }
    else{
        print("starting CV")
        return(replicate(cv.glmnet(snp,exp,nfolds=nfold.set,alpha=alpha,keep=T,parallel=isParallel,standardize=F),n=nrep.set,simplify=F))
    }
}

list_indexer <- function(l,name,simplify=T){
  return(sapply(l,function(x)x[[name]],simplify=simplify))
}

results_parser <- function(cv.list,listel){
    snp <- listel[["snp"]]
    exp <- listel[["exp"]]
    ensid <- listel[["ensid"]]
    genename <- listel[["genename"]]
    if(is.null(cv.list)){
        return(NULL)
    }

    cvm.mins <- apply(list_indexer(cv.list,"cvm","array"),2,min)
    cvm.avg <- mean(cvm.mins)
    cvm.var <- var(cvm.mins)
    each.min <- apply(list_indexer(cv.list,"cvm","array"),2,which.min)
                                        #Find the average of the predictors for the best lambda in each cv.glmnet call
                                        #This line simply takes the columns of fit.preval corresponding to the best lambda for each cv.glmnet, puts them in a matrix, and takes the average of the rows
                                        #Find the most popular lambda (sorta)
    
    nrow.max <- round(mean(each.min))
    best.lambda <- cv.list[[1]][["lambda"]][nrow.max]
    lambdavar <- var(list_indexer(cv.list,"lambda.min"))
    lambda.frac.diff <- sum(list_indexer(cv.list,"lambda.min")!=best.lambda)/length(cv.list)
    pred.avg <- predict(cv.list[[1]],newx=snp,s=best.lambda)
                                        #Find the betas for our favorite lambda
    betadf <- data.frame(beta=cv.list[[1]][["glmnet.fit"]][["beta"]][,nrow.max],gene=genename)
    betadf <- betadf[betadf$beta>0,,drop=F]
    
    res <- summary(lm(exp~pred.avg))

    resultsrow <- data.frame(gene=genename,
                             ensid=ensid,
                             mean.cvm=cvm.avg,
                             var.cvm =cvm.var,
                             lambda.var=lambdavar,
                             lambda.frac.diff=lambda.frac.diff,
                             mean.lambda.iteration=nrow.max,
                             lambda.min=best.lambda,
                             n.snps=nrow(betadf),
                             R2=res$r.squared,
                             alpha=alpha,
                             pval=ifelse(nrow(res$coefficients)>1,
                                 res$coefficients[2,4],
                                 1),
                             stringsAsFactors=F)
    rownames(resultsrow) <- gene
    return(list(resultsrow=resultsrow,betadf=betadf))
}
    

###run LASSO CV

set.seed(1001)

beta.list <- list()
j <- 0
firstresult <- T
firstbeta<-T
resultsname <- paste0(project.name,".",tis,".ResultsArray.",ch,".",seg,".txt")
outfile <- paste0(project.name,".",tis,".BetaTable.",ch,".",seg,".txt")
finishedfile <- paste0(project.name,".",tis,".Finished.",ch,".",seg,".txt")
for(i in 1:ncol(expdata)){
    gene <- colnames(expdata)[i]
    cat(i,"/",ncol(expdata),"\n")
                                        #Pull out the SNP matrix and expression vector corresponding to our gene 
    SNP_EXPlist <- SNP_Exp(gene,snpdata=snpdata,expdata=expdata,snpanno=snpanno,expanno=expanno)
                                        #cv.list will be a data frame where each column is a cv.glmnet result, and each row is a return field
    for (alpha in c(1,0.5,0.95,0.05)){
        cv.list <- glmnet_wrapper(SNP_EXPlist,alpha=alpha,nrep.set=n,nfold.set=k)
        
        print("CV finished")
        parsed_results <- results_parser(cv.list,SNP_EXPlist)

        
        results.df <- parsed_results[["resultsrow"]]
        write.table(results.df,file=resultsname,sep="\t",col.names=firstresult,append=(!firstresult),quote=F,row.names=F)
        firstresult <- F
                                        #cv.list is null if there aren't enough SNPs which satisfy our criteria
        bestbetainfo <- parsed_results[["betadf"]]
        if(!is.null(bestbetainfo)){
            if(nrow(bestbetainfo)>0){
                bestbetainfo <- data.frame(bestbetainfo,snpanno[rownames(bestbetainfo),],alpha=alpha,stringsAsFactors=F)
                write.table(bestbetainfo,file=outfile,sep="\t",quote=F,row.names=F,col.names=firstbeta,append=!firstbeta)
                firstbeta <- F
            }
        }
    }
}

time.stop <- Sys.time()
elapsed <- time.stop-start.time
write(elapsed,file=finishedfile,append=F)
print(elapsed)



