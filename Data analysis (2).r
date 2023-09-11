library("openxlsx")

file <- "/D/Cameroon-standard format1.xlsx"
data <- read.xlsx(file,1)

# data from years 2001 and 2002
head(data)
ncol(data)
# (1:51)[-2]
data1 <- data[is.element(data$year,c("01","02")),(1:51)[-2]]
data1 <- as.data.frame(data1)
head(data1)

markers <- 2:ncol(data1)

## Transform data from 2001/2002 to binary representation per marker
data1.tranf <- data.format(data1,markers)
head(data1.tranf[[1]])


        ### Warmup
        data.comp <- apply(data1.tranf[[1]],1,function(x) paste(x,collapse="-"))
        Nx <- table(data.comp)
        Nx.names <- names(Nx)
        X <- t(sapply(Nx.names, function(x) unlist(strsplit(x,"-")) ))
        rownames(X) <- NULL
        X2 <- as.numeric(X)  # do not use as.integer this causes numerical problems 

        ### Now X2 ans Nx are in the right format

## We need to consider the SNP drug-resistance conferring loci + 1 microsatellite locus; or only the resistance conferrring loci
## Simplest case: we want dhfr and dhps separately!

### dhfr  

DATA <- data1.tranf[[1]][,1:4]
num.alleles <- data1.tranf[[3]][1:4]
list.alleles <- data1.tranf[[2]][1:4]
data.comp <- apply(DATA,1,function(x) paste(x,collapse="-"))
Nx <- table(data.comp)
Nx.names <- names(Nx)
X <- t(sapply(Nx.names, function(x) unlist(strsplit(x,"-")) ))
rownames(X) <- NULL
X2 <- array(as.numeric(X),dim(X))  # do not use as.integer

# Now select only rows without missing data , i.e., no entry =0
sel <- rowSums(X2==0)==0  
Nx1 <- Nx[sel]
X1 <- X2[sel,]
X1
Nx1

estimate <- est(X1,Nx1,num.alleles)
freq.est <- estimate[[2]] ## frequency estimate of haplotypes

hapl.names(freq.est,list.alleles) ## frequency estimates with proper haplotype names

#### dhps 


### 2004/2005

### dhfr + 1 misatellite marker


colnames(data1.tranf[[1]])
DATA <- data1.tranf[[1]][,c(1:4,10)]
num.alleles <- data1.tranf[[3]][c(1:4,10)]
list.alleles <- data1.tranf[[2]][c(1:4,10)]
data.comp <- apply(DATA,1,function(x) paste(x,collapse="-"))
Nx <- table(data.comp)
Nx.names <- names(Nx)
X <- t(sapply(Nx.names, function(x) unlist(strsplit(x,"-")) ))
rownames(X) <- NULL
X2 <- array(as.numeric(X),dim(X))  # do not use as.integer

# Now select only rows without missing data , i.e., no entry =0
sel <- rowSums(X2==0)==0  
Nx1 <- Nx[sel]
X1 <- X2[sel,]
X1
Nx1

estimate <- est(X1,Nx1,num.alleles)
freq.est <- estimate[[2]] ## frequency estimate of haplotypes

freq.est1 <-hapl.names(freq.est,list.alleles) ## frequency estimates with proper haplotype names
rownames(freq.est1)
names1 <- t(sapply(rownames(freq.est1), function(x){
                                        y <- unlist(strsplit(x,"-"))
                                        c(paste(y[-length(y)],collapse="-"),y[length(y)])
                                    }
                             ))

freq.est2 <- as.data.frame(names1)

freq.est2$freq <- freq.est1

pp <- freq.est2[freq.est2[,1]=="ASN-ARG-ASN-ILE",3]
pp <- pp/sum(pp)
1-sum(pp^2)   ### heterozygosity

### For next time:
###         1. get frequency distribution of dhfr haplotypes for 2004/2005 
###         2. get frequency distribution of dhps haplotypes for 2001/2002 and 2004/2005
###         3. get herterozgosity among the dhfr Wildtype, 51/108, 59/108, 51/59/108 mutations for all microsatellite loci
###         4. get herterozgosity among the dhps Wildtype, 436, 437, 436/437 mutations for all microsatellite loci
###         5. Think of how to get bootstrap confidence intervalls for heterozygosity

table(as.data.frame(names1))
library(questionr)
wtd.table(names1[,1],names1[,2],freq.est1)
as.integer(unlist(strsplit("1-2-3","-")))


paste(c(1,2,3),collapse="-")
paste(c(1,2,3),3,4,3,"dfvfev",sep="-")

# data from years 2004 and 2005
head(data)
ncol(data)
(1:51)[-2]
data2 <- data[is.element(data$year,c("04","05")),(1:51)[-2]]
data2 <- as.data.frame(data2)

