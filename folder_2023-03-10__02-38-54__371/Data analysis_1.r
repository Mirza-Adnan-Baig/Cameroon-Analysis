library("openxlsx")
library(ggplot2)
file <- "/Users/kristanschneider/Meine Ablage (mathmalaria@gmail.com)/Lehre/Research in Applications/Cameroon-standard format.xlsx"
path <- "/Users/kristanschneider/Meine Ablage (mathmalaria@gmail.com)/Lehre/Research in Applications/Plots/"
data <- read.xlsx(file,1)
colnames(data) <- colnames(read.xlsx(file,2))

# data from years 2001 and 2002
head(data)
ncol(data)
# (1:51)[-2]
data1 <- data[is.element(data$year,c("01","02")),(1:51)[-2]]
data1 <- as.data.frame(data1)
head(data1)

markers <- 2:ncol(data1)

## relative allele frequencies per table
allele.freqs <- apply(data1[,-1],2,function(x) prop.table(table(x)))
### visualistions 

barplot(allele.freqs[[15]])  ## this is not nice -> package ggplot



allele.freqs1 <- lapply(allele.freqs,data.frame)
for(k in 1:length(allele.freqs1)){
        cbPalette <- c( "#0072B2", "#E69F00" , "#009E73","#56B4F9", "#CC79A7")
        p <- ggplot(data=allele.freqs1[[k]], aes(x=as.factor(x), y=Freq)) +
            geom_bar(stat="identity",position=position_dodge(), colour="black",fill= "#0072B2") 
        p <- p + theme(panel.background = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(colour='black',fill=NA),
                        axis.text = element_text(size = rel(1.3),color='black'),
                        axis.text.x = element_text(angle=75,vjust=0.5),
                        axis.title = element_text(size = rel(1.3)),
                        plot.title = element_text(size = rel(1.4),color='black',hjust=0.5),
                        legend.text = element_text(size = rel(1.3)),
                        legend.title = element_text(size = rel(1.3)))
        p <- p + scale_fill_manual(values=cbPalette,name="Years") + ylim(0,1)
        p <- p + labs(x="alleles",y="frequencies",title=parse(text=names(allele.freqs1)[[k]]))
        p
        outpath <- paste(path,"Allele_freqs_",names(allele.freqs1)[[k]],".pdf",sep="")  ## path where to put file
        pdf(outpath,height=5)
        print(p)
        dev.off()
}



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


### dhfr   + dhps 

DATA <- data1.tranf[[1]][,c(1:3,5,6,8,9)]  ## dhfr 51, 59, 108, dhps 436, 437, 581, 613
num.alleles <- data1.tranf[[3]][c(1:3,5,6,8,9)]
list.alleles <- data1.tranf[[2]][c(1:3,5,6,8,9)]
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

hapl.freqs <- hapl.names(freq.est,list.alleles) ## frequency estimates with proper haplotype names
rownames(hapl.freqs)[hapl.freqs>0.05]




#### dhfr + dhps + 1 microsatellite marker

# We want the heterozygosity among those drug-resistant haplotypes which have frequency > 1%
hapllist <- rownames(hapl.freqs)[hapl.freqs>0.01]
hapllist

# 5 dhfr/dhps haplotypes have a frequency > 1% 
# We want heterozygosity for each marker on the backgroung of each of the 5 haplotypes

SNPs <- c(1:3,5,6,8,9)
STR.markers <- 10:49
out.het <- array(,c(length(hapllist),length(STR.markers)))  # matrix to store the heterozygosity values rows - haplotype backgroung, columns are the microsatellite markers
rownames(out.het) <- hapllist

for(kk in 1:length(STR.markers)){ # we go first marker by marker
    k <- STR.markers[kk]  # the column in the original data
    colnames(data1.tranf[[1]])
    DATA <- data1.tranf[[1]][,c(SNPs,k)]  ## the SNP markers at dhfr/dhps + one microsatellite marker
    num.alleles <- data1.tranf[[3]][c(SNPs,k)]
    list.alleles <- data1.tranf[[2]][c(SNPs,k)]
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

    freq.est1 <- hapl.names(freq.est,list.alleles) ## frequency estimates with proper haplotype names
    rownames(freq.est1)
    # next remove the microsatellite from the haplotype names
    names1 <- t(sapply(rownames(freq.est1), function(x){
                                            y <- unlist(strsplit(x,"-"))
                                            c(paste(y[-length(y)],collapse="-"),y[length(y)])
                                        }
                                ))

    freq.est2 <- as.data.frame(names1)  # haplotype freqeuncy estimates

    freq.est2$freq <- freq.est1 
    for(h in hapllist){ # we calculate the heterozygosity for conditioned on each of the 5 dhfr/dhps haplotype of interest
        pp <- freq.est2[freq.est2[,1]==h,3] # allele frequencies at the microsat locus given dhfr/dhps haplotype
        pp <- pp/sum(pp) # conditional relative frequencies
        out.het[h,kk] <- 1-sum(pp^2)   ### heterozygosity
    }
    
    

}

plot(out.het[1,]/out.het[5,],type="l")

plot(out.het[11,],type="l")
colnames(dat0)


### dhfr + 1 misatellite marker

hapllist <- c("ASN-ARG-ASN-ILE","ASN-CYS-SER-ILE","ILE-ARG-ASN-ILE")
STR.markers <- 10:49
out.het <- array(,c(length(hapllist),length(STR.markers)))
rownames(out.het) <- hapllist
colnames(out.het) <- colnames(data1.tranf[[1]])[STR.markers]

for(kk in 1:length(STR.markers)){
    k <- STR.markers[kk]
    colnames(data1.tranf[[1]])
    DATA <- data1.tranf[[1]][,c(1:4,k)]
    num.alleles <- data1.tranf[[3]][c(1:4,k)]
    list.alleles <- data1.tranf[[2]][c(1:4,k)]
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

    freq.est1 <- hapl.names(freq.est,list.alleles) ## frequency estimates with proper haplotype names
    rownames(freq.est1)
    names1 <- t(sapply(rownames(freq.est1), function(x){
                                            y <- unlist(strsplit(x,"-"))
                                            c(paste(y[-length(y)],collapse="-"),y[length(y)])
                                        }
                                ))

    freq.est2 <- as.data.frame(names1)

    freq.est2$freq <- freq.est1
    for(h in hapllist){
        pp <- freq.est2[freq.est2[,1]==h,3]
        pp <- pp/sum(pp)
        out.het[h,kk] <- 1-sum(pp^2)   ### heterozygosity
    }
    
    

}

plot(out.het[1,1:18]/out.het[2,1:18],type="l")  # het of 59/108 mutant relqtiv eto wildtype
plot(out.het[3,1:18]/out.het[2,1:18],type="l")  # het of 51/59/108 mutant relqtiv eto wildtype

summary(data.frame(data1,stringsAsFactors = FALSE))
ddd <- data.frame(data1,stringsAsFactors = TRUE)
class(ddd[,15])
is.data.frame(data1)
summary(data1)

data2 <- data1
data2[,] <- lapply(data2,factor)
summary(data2)

DATA <- data1.tranf[[1]][,]
num.alleles <- data1.tranf[[3]][1:4]
list.alleles <- data1.tranf[[2]][1:4]

library(ISwR)
summary(lung)
class(lung[,2])
is.data.frame(lung[,2])
summary(data1)
colnames(dat0)



### dhfr + 2 misatellite markers

hapllist <- c("ASN-ARG-ASN-ILE","ASN-CYS-SER-ILE","ILE-ARG-ASN-ILE")
STR.markers <- 10:49
m <- length(STR.markers)
out.LD <- array(,c(length(hapllist),m*(m-1)/2))
rownames(out.LD) <- hapllist
#colnames(out.het) <- colnames(data1.tranf[[1]])[STR.markers]

m <- 4

for(kk in 2:m){
    for(jj in 1:(kk-1) )
    k <- STR.markers[kk]
    j <- STR.markers[jj]

    k <- 11
    j <- 10
    colnames(data1.tranf[[1]])
    DATA <- data1.tranf[[1]][,c(1:4,j,k)]
    num.alleles <- data1.tranf[[3]][c(1:4,j,k)]
    list.alleles <- data1.tranf[[2]][c(1:4,j,k)]
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

    freq.est1 <- hapl.names(freq.est,list.alleles) ## frequency estimates with proper haplotype names
    rownames(freq.est1)
    names1 <- t(sapply(rownames(freq.est1), function(x){
                                            y <- unlist(strsplit(x,"-"))
                                            c(paste(y[1:(length(y)-2)],collapse="-"),y[length(y)-1],y[length(y)])
                                        }
                                ))

    freq.est2 <- as.data.frame(names1)
    print(freq.est2)
    freq.est2$freq <- freq.est1
    #for(h in hapllist){
    #    pp <- freq.est2[freq.est2[,1]==h,3]
    #    pp <- pp/sum(pp)
    #    out.het[h,kk] <- 1-sum(pp^2)   ### heterozygosity
    #}
    
    

}

plot(out.het[1,1:18]/out.het[2,1:18],type="l")  # het of 59/108 mutant relqtiv eto wildtype
plot(out.het[3,1:18]/out.het[2,1:18],type="l")  # het of 51/59/108 mutant relqtiv eto wildtype

summary(data.frame(data1,stringsAsFactors = FALSE))
ddd <- data.frame(data1,stringsAsFactors = TRUE)
class(ddd[,15])
is.data.frame(data1)
summary(data1)

data2 <- data1
data2[,] <- lapply(data2,factor)
summary(data2)

DATA <- data1.tranf[[1]][,]
num.alleles <- data1.tranf[[3]][1:4]
list.alleles <- data1.tranf[[2]][1:4]

library(ISwR)
summary(lung)
class(lung[,2])
is.data.frame(lung[,2])
summary(data1)
colnames(dat0)
### For next time:
###         1. get frequency distribution of dhfr haplotypes for 2004/2005 
###         2. get frequency distribution of dhps haplotypes for 2001/2002 and 2004/2005
###         3. get herterozgosity among the dhfr Wildtype, 51/108, 59/108, 51/59/108 mutations for all microsatellite loci
###         4. get herterozgosity among the dhps Wildtype, 436, 437, 436/437 mutations for all microsatellite loci
###         5. Think of how to get bootstrap confidence intervalls for heterozygosity





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


### 2004/2005
