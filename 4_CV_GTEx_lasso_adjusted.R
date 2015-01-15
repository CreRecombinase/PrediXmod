####script file found in /nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/
####by Heather E. Wheeler 20140602####
args <- commandArgs(trailingOnly=T)
date <- Sys.Date() 
###############################################
### Directories & Variables


project.name <- args[1]
my.dir <- args[2]
tissue <- args[3]
ch <- args[4]
seg <- args[5]

tis <- gsub(" ","",tissue)


setwd(my.dir)
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
  

library(glmnet)
library(doMC)
registerDoMC(cores=2)


stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(x))
lower <- function(x) quantile(x,0.025,na.rm=TRUE)
upper <- function(x) quantile(x,0.975,na.rm=TRUE)

## convenience function to select best lambda over cv bootstraps for model linear by Keston edited by Heather to get predicted values
glmnet.select <- function(response, covariates, nrep.set = 10, nfold.set = 10, alpha.set = 1, ...) {
    require(glmnet)
    best.lam.sim = vector()
    best.cvm.sim = vector()
    pred.matrix = matrix(0,nrow=nrow(covariates),ncol=nrep.set)
    for (i in 1:nrep.set) {
        glmnet.fit <- cv.glmnet(covariates, response, nfolds = nfold.set, alpha = alpha.set, keep = TRUE,parallel=T)
        best.lam.sim[i] <- which.min(glmnet.fit$lambda)
        best.cvm.sim[i] <- min(glmnet.fit$cvm)
        pred.matrix[,i] <- glmnet.fit$fit.preval[,best.lam.sim[i]]  #Betas from iteration of cv.glmnet
    }
    cvm.avg = mean(best.cvm.sim) # average cvm
    nrow.max = as.integer(mean(best.lam.sim)) # best lambda over cv bootstraps
    ret <- as.data.frame(glmnet.fit$glmnet.fit$beta[,nrow.max])
    ret[ret == 0.0] <- NA
    ret.vec = as.vector(ret[!is.na(ret),]) # vector of non-zero betas
    names(ret.vec) = rownames(ret)[!is.na(ret)]
    min.lambda <- glmnet.fit$glmnet.fit$lambda[nrow.max]
    pred.avg <- rowMeans(pred.matrix)
    output = list(ret.vec=ret.vec, cvm.avg=cvm.avg, nrow.max=nrow.max, min.lambda=min.lambda, pred.avg=pred.avg)
    return(output)
}

 




###input genotype data###
 
###create results array
resultsarray <- array(0,c(ncol(expdata),7))
dimnames(resultsarray)[[1]] <- colnames(expdata)
dimnames(resultsarray)[[2]] <- c("gene","mean.cvm","mean.lambda.iteration","lambda.min","n.snps","R2","pval")

###run LASSO CV

set.seed(1001)

for(i in 1:ncol(expdata)){
    gene <- colnames(expdata)[i]
    cat(i,"/",ncol(expdata),"\n")

    start <- max(expanno[gene,"genestart"]-1e6,0)
    end <- expanno[gene,"genestop"]+1e6
    cissnpdata <- snpdata[,snpanno$rsid[(snpanno$pos>=start&snpanno$pos<=end)]]
    cissnpdata[cissnpdata >= 3] <- NA ###.gds files have missing=3
                                        #    if(is.null(dim(cisgenos)) | dim(cisgenos)[2] == 0){###effectively skips genes with <2 cis-SNPs
    if(is.null(dim(cissnpdata))){
        bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
    }else if(ncol(cissnpdata) == 0){
        bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
    }else{
        cissnpdata <- cissnpdata[,apply(cissnpdata,2,function(x)mean(x,na.rm=T)>0)]
        cissnpdata <- scale(cissnpdata,center=T,scale=T)
        cissnpdata[is.na(cissnpdata)] <- 0
	if(is.null(dim(cissnpdata))){
            bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
        }else if(ncol(cissnpdata) == 0){
            bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
        }else{
            exppheno <- expdata[,gene]
            exppheno <- scale(exppheno, center=T, scale=T)  ###need to scale for fastLmPure to work properly
            exppheno[is.na(exppheno)] <- 0
            rownames(exppheno) <- rownames(expdata)
            print("starting CV")
            
            cv <- glmnet.select(exppheno,cissnpdata,nrep.set=n,nfold.set=k,alpha.set=1) ###run lasso k-fold CV n times to determine best lambda & betas
            print("CV finished")

            bestbetas <- cv[["ret.vec"]] ###how many SNPs in best predictor?
        }
    }
    if(length(bestbetas) > 0){
        pred.lasso <- cv[["pred.avg"]] ###mean k-fold CV predictions from n reps

        ### calculate correlation between predicted and observed expression
        res <- summary(lm(exppheno~pred.lasso))
        genename <- expanno[gene,"id"]
        resultsarray[gene,1] <- genename
        resultsarray[gene,2] <- cv[["cvm.avg"]] 
        resultsarray[gene,3] <- cv[["nrow.max"]] ###add mean of best lambda iteration to results
        resultsarray[gene,4] <- cv[["min.lambda"]] ###add best lambda to results
        resultsarray[gene,5] <- length(bestbetas) ###add #snps in prediction to results
        resultsarray[gene,6] <- res$r.squared ###lm R2
        resultsarray[gene,7] <- res$coefficients[2,4] ###lm p-value

        ### output bestbetas for PrediXcan
        bestbetalist <- names(bestbetas)
        bestbetainfo <- snpanno[bestbetalist,]
        bestbetainfo$beta <- bestbetas
        bestbetainfo$gene <- gene
        outfile <- paste0(project.name,".",tis,".BetaTable.",ch,".",seg,".txt")
        cat("Writing file")
        write.table(bestbetainfo, file=outfile,quote=F,row.names=F,sep="\t",col.names=(i==1),append=(i!=1))
    }else{
	genename <- gene
        resultsarray[gene,1] <- genename
        resultsarray[gene,2:7] <- c(NA,NA,NA,0,NA,NA)
    }
}

saveRDS(resultsarray,paste0(project.name,".",tis,".ResultsArray.",ch,".",seg,".RDS"))



