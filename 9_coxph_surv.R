library(survival)
library(riskRegression)
library(dplyr)



chrs <- 1:22
alphas <- c("0.5")
Tissues <- c("Nerve-Tibial","Thyroid")
impdir <- "/group/dolan-lab/nknoblauch/CALGB40101/imputedResults/"
outdir <- "/group/dolan-lab/nknoblauch/CALGB40101/survdir/"
phenofile <- "/group/dolan-lab/nknoblauch/CALGB40101/phenodata/CALGB40101.G2.coxph.phen"
phenodat <- read.table(phenofile,header=T,stringsAsFactors=F,sep="\t")

badphenos <- scan("/group/dolan-lab/nknoblauch/CALGB40101/phenodata/badphenotypes.txt",what=character())


rownames(phenodat) <- phenodat[,1]
goodphenodata <- phenodat[!rownames(phenodat) %in% badphenos,]
msurv <- Surv(phenodat$dose,phenodat$event)
gmsurv <- Surv(goodphenodata$dose,goodphenodata$event)


for(i in chrs){
    for(j in alphas){
        for(k in Tissues){
            print(i)
            print(j)
            print(k)
            texpfile <- paste0(impdir,"chr",i,".",k,".",j,".txt")
            if(!file.exists(texpfile)){
                print(paste0("input file not found!: ",texpfile))
                stop()
            }
            toutfile <- paste0(outdir,"chr",i,".",k,".",j,"QC_cox.txt")
            expdat <- read.table(texpfile,header=F,stringsAsFactors=F,sep="\t")
            rownames(expdat) <- expdat[,1]
            expdat <- expdat[,-1]
            expdat <- expdat[,!rownames(phenodat) %in% badphenos]
            expdat <- data.matrix(expdat)
            acoeffs <- t(apply(expdat,1,function(x)coef(summary(coxph(gmsurv~x)))))
            colnames(acoeffs) <- c("coef","exp(coef)","se(coef)","z","Pr(>|z|)")
            dcoeffs <- data.frame(gene=rownames(acoeffs),acoeffs)
            write.table(dcoeffs,toutfile,col.names=T,row.names=T,quote=F,sep="\t")
        }
    }
}
    

