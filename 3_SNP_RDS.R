
args <- commandArgs(trailingOnly=T)
project.name <- args[1]
my.dir <- args[2]

setwd(my.dir)
files <- paste0(project.name,c(".SNP.list",
                               ".SNPxID",
                               ".ID.list",
                               ".bim"))
names(files) <- c("SNPs","SNPMatrix","IDs","BIM")

bimdata <- read.table(files["BIM"],header=F,stringsAsFactors=F,sep="\t")
colnames(bimdata) <- c("chrom","rsid","z","pos","ref","alt")
saveRDS(bimdata,paste0(project.name,".SNPanno.RDS"))

snps <- scan(files["SNPs"],what="character",sep="\n")
IDs <- scan(files["IDs"],what="character",sep="\n")

snpdata <- matrix(scan(files["SNPMatrix"],what=numeric()),ncol=length(IDs),byrow=T)
print(dim(snpdata))
print(length(snps))
print(length(IDs))
rownames(snpdata) <- snps
colnames(snpdata) <- IDs
saveRDS(snpdata,paste0(project.name,".SNPxID.RDS"))

