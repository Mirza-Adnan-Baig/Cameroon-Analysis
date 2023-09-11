library("openxlsx")
#library(bignum)
library(Rmpfr)
file <- "/Users/kristanschneider/Meine Ablage (mathmalaria@gmail.com)/Lehre/Research in Applications/Cameroon-standard format.xlsx"
dat0 <- read.xlsx(file,1)


head(dat0) ##
dim(dat0)
markers <- 3:ncol(dat0)
head(dat0) ##
 #markers <- c(3:5,7:8,10:15) # ncol(dat)
 #unique(dat0[,9])
dat <- dat0

### function data.format outputs a list containing as first element  the data in new format (1 row per sample, 1 colum by marker) 
### entries are integerers. If transformed into binary numbers 0-1 vectors they indicate absence/presence of alleles 
### e.g., 0.... no allele present, 1... first allele present, 2... 2nd allele present, 3.... 1st and 2nd allele present ect. 
### The second element gives the order of the alleles pper locus
### The third element gives th number of alleles per locucs

data.format <- function(dat,markers){
  ### dat... is the input data set in standard format of package MLMOI, 1 st column contains smaple IDs
  ### markers ... vector of columms containing markers to be included
  
  ### list of alleles per marker###############
  allele.list <- apply(dat[,markers],2, function(x){ 
                          y=sort(unique(x))
                          y[!is.na(y)]
                        } 
                      )
  ###number of alleles per marker########
  allele.num <- unlist(lapply(allele.list,length))  
  test <- min(allele.num < 17) ## determines whether we need high precision arithmetics
  #### split data by sample ID
  dat.split <- split(dat[,markers],dat[,1])  #splits data by sample IDs
  l <- length(markers)
  if(test){ ### no high precision needed
    #### Binary representation of allele being absent and present
    samples.coded <- t(sapply(dat.split, function(x){mapply(function(x,y,z){as.integer(is.element(y,x)) %*% 2^(0:(z-1))}, x , allele.list, allele.num)}))
  }
  else{ #high precision is needed
    #### Binary representation of allele being absent and present
    N <- length(dat.split)
    m <- sum(allele.num)   
    ## we need 2^m numbers  2^m = 10^x m log 2 = x log 10 -> x = m* log(2)/log(10) -> x = ceiling (m* log(2)/log(10))
    ##m1 <- ceiling(m * log(2)/log(10))  ## precision we need
    samples.coded <- mpfrArray(NA,m,c(N,l))
    for(n in 1:N){
      for(k in 1:l){
        samples.coded[n,k] <- as.integer(is.element(allele.list[[k]],unlist(dat.split[[n]][k]))) %*%  mpfr(2^(0:(allele.num[[k]]-1)),m)   # binary coding
      }
    }
  }
  list(samples.coded,allele.list,allele.num)
  
}

### This function outputs a list
#### 1st gives the compact notation for the data
#### 2nd element base factos 1 g1 g1g2. ....,g1*...g(l-1)
#### 3rd element has the number of alleles in the prope order per locus
#### 4th element is the number of alleles per locus
data.formal.AL <- function(data){
  l <- length(data[[3]])   ## number of loci
  m <- sum(data[[3]])
  if(m <16){ # We need no high precision
    basefactors <- cumprod(c(as_biginteger(1),2^(as_biginteger(data[[3]][-l]))))
    M <- data[[1]] %*% basefactors
    print("low precision")
  }else{  # We need  high precision
    print("high precision")
    #m1 <- m*log(2)/log(10)
    basefactors <- cumprod(c(mpfr(1,m),2^(mpfr(data[[3]][-l],m))))
    M <- data[[1]] %*% basefactors
    
  }
  list(M, basefactors,data[[2]],data[[3]])
  
}

### Going back
div <- function(X,ff){
  x <- NULL
  for(k in length(ff):1){
    x <- c(X %/% ff[k],x) 
    X <- X %% ff[k]
  }
  x
}

reconst.data <- function(data){
  if(class(data[[1]])!="mpfrMatrix"){ # low precision
    newdat <- t(sapply(data[[1]], function(x) div(as_biginteger(x),data[[2]]) ) )
  }
  else{ # high precision is necessary
    N <- dim(data[[1]])[1]  # sample size
    L <- length(dat[[2]])   # number of loci
    m <- getPrec(dat[[1]][1]) # precision
    newdat <-mpfrArray(0,m,c(N,L))
    for(n in 1:N){
      newdat[n,] <- div(data[[1]][n],data[[2]])
    }
  }
  list(newdat,data[[3]],data[[4]])
}

lapply(1:10,function(x) list())
rec[[3]]
rec[[1]][1,1]
div(rec[[1]][1,1],2^(0:1))==1
reconst.data.long <- function(data){
  if(class(dat[[1]])!="mpfrMatrix"){
    newdat <- t(sapply(data[[1]], function(x) div(as_biginteger(x),data[[2]]) ) )
  }
  else{
    d <-  dim(data[[1]])
    N <- d[1]  # sample size
    L <- d[2]   # number of loci
    outlist <- lapply(1:N,function(x) list())
    
    for(l in 1:L){
      base <- 2^(0:(data[[3]][l]-1))
      for(n in 1:N){
        outlist[[n]][[l]] <- data[[2]][div(data[[1]][n,l],base)==1]
      }
    }
  }
  list(outlist,data[[2]],data[[3]])
}
reconst.data.long(rec)
# Try

datanew <- data.format(dat0,markers)
head(datanew[[1]])
is.list(datanew[[1]])
unlist(datanew[[1]])

dat <- data.formal.AL(datanew)
rec <- reconst.data(dat)
as.integer(rec[[1]][1,])
array(round(reconst.data(dat)[[1]]),dim(rec[[1]]))


head(rec[[1]])


