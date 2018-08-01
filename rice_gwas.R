###Rice_GWAS
###Genotype data
# install.packages("BGLR")
library(BGLR)
out<- read_ped("/Users/zyang35/Documents/yanglab/work/GWAS/RiceDiversity_44K_Genotypes_PLINK/sativas413.ped")
p=out$p
n=out$n
out=out$x
#Recode snp to 0,1,2 format using allele 1
# 0 --> 0
# 1 --> 1
# 2 --> NA
# 3 --> 2
out[out==2]=NA
out[out==3]=2
W <- matrix(out, nrow=p, ncol=n, byrow=TRUE)
W <- t(W) 
dim(W) # # 413 x 36901

###accession ID
fam <-read.table("/Users/zyang35/Documents/yanglab/work/GWAS/RiceDiversity_44K_Genotypes_PLINK/sativas413.fam", header = FALSE, stringsAsFactors = FALSE)  
head(fam)
rownames(W) <- fam$V2 # 413 x 36901
#dim(fam)
#tail(fam)


### phenotypes
rice.pheno <- read.table("http://www.ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt", header=TRUE, stringsAsFactors = FALSE, sep = "\t")
table(rownames(W) == rice.pheno$NSFTVID)
y <- matrix(rice.pheno$Flowering.time.at.Arkansas) # # use the first trait 
rownames(y) <- rice.pheno$NSFTVID
index <- !is.na(y)

y <- y[index, 1, drop=FALSE] # 374
W <- W[index, ] # 374 x 36901
table(rownames(W) == rownames(y))


###Population structure

# PC plots
gp <-read.csv("http://ricediversity.org/data/sets/44kgwas/RiceDiversity.44K.germplasm.csv", header = TRUE, skip = 1,  stringsAsFactors = FALSE)   # 431(? 413?) x 12
gp2 <- gp[match(rownames(y), gp$NSFTV.ID), ]
table(rownames(y) == gp2$NSFTV.ID)

plot(gp2$PC1, gp2$PC2, xlab="PC1", ylab="PC2", col=c(1:6)[factor(gp2$Sub.population)])
legend(x="topleft", legend = levels(factor(gp2$Sub.population)), col=c(1:6), pch=1, cex=0.6)

#ver=c(1:6)[factor(gp2$Sub.population)]
#verr=factor(gp2$Sub.population)
#verrr=levels(factor(gp2$Sub.population))
#ver
#verr
#verrr

###Genotype imputation

for (j in 1:ncol(W)){
  W[,j] <- ifelse(is.na(W[,j]), mean(W[,j], na.rm=TRUE), W[,j])
}

###Compute allele frequencies for all SNPs
p <- colSums(W) / (2 * nrow(W)) # or colMean(X) / 2

###Minor allele frequency (MAF)
maf <- ifelse(p > 0.5, 1-p, p)
maf.index <- which(maf < 0.05)
W2 <- W[, -maf.index]  # 374 x 33558 



###GWAS
library(rrBLUP)
map <- read.table("/Users/zyang35/Documents/yanglab/work/GWAS/RiceDiversity_44K_Genotypes_PLINK/sativas413.map", header = FALSE, stringsAsFactors = FALSE)
my.geno <- data.frame(marker=map[,2], chrom=map[,1], pos=map[,4], t(W-1), check.names = FALSE) # W = \in{-1, 0, 1}
my.pheno <- data.frame(NSFTV_ID=rownames(y), y=y) 

rel <- GWAS(my.pheno, my.geno, min.MAF=0.05, P3D=TRUE, plot=FALSE)
head(rel$y)
tail(rel$y)

library(qqman)
manhattan(x = rel, chr = "chrom", bp = "pos", p = "y", snp = "marker", col = c("blue4", "orange3"), logp = FALSE)

qq(rel$y, main = "Q-Q plot of GWAS p-values", xlim = c(0, 7), ylim = c(0, 7), pch = 18, col = "blue4", cex = 1.5, las = 1)