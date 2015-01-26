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
n <- 10 #number of k-fold CV replicates
 
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
    return(list(exp=exppheno,snp=cissnpdata,genename=genename))
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
        return(replicate(cv.glmnet(snp,exp,nfolds=nfold.set,alpha=alpha,keep=T,parallel=isParallel,standardize=F),n=nrep.set))
    }
}

results_parser <- function(cv.list,exp,genename,ensid){
    if(is.null(cv.list)){
        return(NULL)
    }
    cvm.avg <- mean(sapply(cv.list["cvm",],min))
                                        #Find the index of the lambda with the lowest cvm for each cv.glmnet call
    each.min <- sapply(cv.list["cvm",],which.min)

                                        #Find the average of the predictors for the best lambda in each cv.glmnet call
                                        #This line simply takes the columns of fit.preval corresponding to the best lambda for each cv.glmnet, puts them in a matrix, and takes the average of the rows
    pred.avg <- rowMeans(mapply(
        function(x,y){
            x[,y]
        },
        cv.list["fit.preval",],
        each.min
    )
                         )
                                        #Find the most popular lambda (sorta)
    nrow.max <- as.integer(mean(each.min))        
    best.lambda <- cv.list["lambda",1][[1]][nrow.max]
    
                                        #Find the betas for our favorite lambda
    betadf <- data.frame(beta=cv.list["glmnet.fit",1][[1]][["beta"]][,nrow.max],gene=genename)
    betadf <- betadf[betadf$beta>0,,drop=F]
    
    res <- summary(lm(exp~pred.avg)) 
    resultsrow <- data.frame(gene=genename,
                             ensid=ensid,
                             mean.cvm=cvm.avg,
                             mean.lambda.iteration=nrow.max,
                             lambda.min=best.lambda,
                             n.snps=nrow(betadf),
                             R2=res$r.squared,
                             alpha=alpha,
                             pval=res$coefficients[2,4],stringsAsFactors=F)
    rownames(resultsrow) <- gene
    return(list(resultsrow=resultsrow,betadf=betadf))
}
    

###run LASSO CV

set.seed(1001)
results.list <- list()
j <- 0
for(i in 1:ncol(expdata)){
    gene <- colnames(expdata)[i]
    cat(i,"/",ncol(expdata),"\n")
                                        #Pull out the SNP matrix and expression vector corresponding to our gene 
    SNP_EXPlist <- SNP_Exp(gene,snpdata=snpdata,expdata=expdata,snpanno=snpanno,expanno=expanno)
                                        #cv.list will be a data frame where each column is a cv.glmnet result, and each row is a return field
    for (alpha in c(1,0.5,0.95)){
        j <- j+1
        cv.list <- glmnet_wrapper(SNP_EXPlist,alpha=alpha,nrep.set=n,nfold.set=k)
        
        print("CV finished")
        parsed_results <- results_parser(cv.list,SNP_EXPlist[["exp"]],SNP_EXPlist[["genename"]],gene)
        results.list[[j]] <- parsed_results[["resultsrow"]]
                                        #cv.list is null if there aren't enough SNPs which satisfy our criteria

    
        bestbetainfo <- parsed_results[["betadf"]]
        if(!is.null(bestbetainfo)){
            if(nrow(bestbetainfo)>0){
                bestbetainfo <- data.frame(bestbetainfo,snpanno[rownames(bestbetainfo),],alpha=alpha,stringsAsFactors=F)
                
                outfile <- paste0(project.name,".",alpha,".",tis,".BetaTable.",ch,".",seg,".txt")
                write.table(bestbetainfo, file=outfile,quote=F,row.names=F,sep="\t",col.names=(j==1),append=(j!=1))
            }
        }
    }
}
resultsarray <- rbind.fill(results.list)
time.stop <- Sys.time()
elapsed <- time.stop-start.time
print(elapsed)
saveRDS(resultsarray,paste0(project.name,".",alpha,".",tis,".ResultsArray.",ch,".",seg,".RDS"))



