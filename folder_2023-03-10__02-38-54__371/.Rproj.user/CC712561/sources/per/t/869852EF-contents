library("openxlsx")

file <- "/Users/kristanschneider/Meine Ablage (mathmalaria@gmail.com)/Lehre/Research in Applications/Cameroon-standard format.xlsx"
dat0 <- read.xlsx(file,1)


head(dat0) ##
markers <- 3:ncol(dat0)
head(dat0) ##
markers <- 3:20 # ncol(dat)

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
    #### split data b sample ID
    dat.split <- split(dat[,markers],dat[,1])
    #### Binary representation of allele being absent and present
    samples.coded <- t(sapply(dat.split, function(x){mapply(function(x,y,z){as.integer(is.element(y,x)) %*% 2^(0:(z-1))}, x , allele.list, allele.num)}))
    list(samples.coded,allele.list,allele.num)
}



### This function outputs a list
#### 1st gives the compact notation for the data
#### 2nd element base factos 1 g1 g1g2. ....,g1*...g(l-1)
#### 3rd element has the number of alleles in the prope order per locus
#### 4th element is the number of alleles per locus
dafa.formal.AL <- function(data){
        l <- length(data[[3]])   ## number of loci
        basefactors <- cumprod(c(1,2^(data[[3]][-l])))
        list(data[[1]] %*% basefactors, basefactors,data[[2]],data[[3]])
}



### Going back


div <- function(X,f){
    x <- NULL
    for(f in rev(dat[[2]])){
        #print(c(X,f))
        x <- c(X %/% f,x) 
        X <- X%%f
    }
    x
}

reconst.data <- function(data){
    newdat <- t(sapply(data[[1]], function(x) div(x,data[[2]]) ) )
    list(newdat,data[[3]],data[[4]])
}

# Try


datanew <- data.format(dat0,markers)
head(datanew[[1]])
datanew[[3]]

l <- length(datanew[[3]])   ## number of loci

basefactors <-  cumprod(c(1,2^(datanew[[3]][-l])))

datanew[[1]] %*% basefactors



dat <- dafa.formal.AL(datanew)



reconst.data(dat)
dat
out <-table(dat[[1]])
inpt <- list(as.integer(names(out)),dat[[2]],dat[[3]],dat[[4]])
reconst.data(inpt)



head(t(sapply(dat[[1]], function(x) div(x) ) ))
head(datanew[[1]])
dat[[2]]

datanew[[1]][1:2,] %*% cumprod(c(1,2^datanew[[3]][1:7]))
div <- function(X,f){
    x <- NULL
    for(f in rev(dat[[2]])){
        #print(c(X,f))
        x <- c(X %/% f,x) 
        X <- X%%f
    }
    x
}
cumprod(c(1,2^datanew[[3]][1:7]))
div( 6870,cumprod(c(1,2^datanew[[3]][1:7])))

6870 %/% 4096 
6870 %% 4096 
2774 %/% 2048
2774 %% 2048

720%%512

#R MPFR - Multiple Precision Floating-Point Reliable
65 %/% 7
65 %% 7##

X <-dat[[1]][1]
x <- NULL
for(f in rev(dat[[2]])){
    #print(c(X,f))
    x <- c(X %/% f,x) 
    X <- X%%f
}

 x

 datanew[[1]][1,]
#################### .   ended on the cutting room floor - backup


library("openxlsx")

file <- "/Users/kristanschneider/Meine Ablage (mathmalaria@gmail.com)/Lehre/Research in Applications/Cameroon-standard format.xlsx"
dat <- read.xlsx(file,1)

head(dat) 


### list of alleles per marker###############
markers <- 3:ncol(dat)
allele.list <- apply(dat[,markers],2, function(x){ 
    y=sort(unique(x))
    y[!is.na(y)]
    } )

allele.num <- unlist(lapply(allele.list,length))  ###number of alleles per marker########

dat.split <- split(dat[,markers],dat[,1])

samples.coded <- t(sapply(dat.split, function(x){mapply(function(x,y,z){as.integer(is.element(y,x)) %*% 2^(0:(z-1))}, x , allele.list, allele.num)}))

head(samples.coded)


as.matrix(unlist(samples.coded,recursive=TRUE,use.names = FALSE),nrows= length(markers),ncol= length(unique(dat[,1])))



unlist(samples.coded,recursive=TRUE,use.names = FALSE)


length(c(unlist(samples.coded)))
head(t(array(c(unlist(samples.coded)),c(50,331))))
allele.num[[1]]

is.element(NA,1:10)


lapply(dat.split, function(x){dim(x)})


mapply(function(x,y,z){
    as.integer(is.element(y,x) %*% 2^(0:(z-1)))
}, dat.split , allele.list, allele.num)



dat.split <- split(dat,dat[,1])[[300]]


mapply(function(x,y,z){is.element(y,x) %*% 2^(0:(z-1))}, dat.split , allele.list, allele.num)


mapply(function(x,y,z){as.integer(is.element(y,x)) %*%  2^(0:(z-1))}, dat.split , allele.list, allele.num)

mapply(function(x,y,z){as.integer(as.integer(is.element(y,x)) %*%  2^(0:(z-1)))}, dat.split , allele.list, allele.num)


as.integer(is.element(allele.list[[1]],dat.split[,1])) %*% 2^(0:(allele.num[[1]]-1))

2^(0:(allele.num[[27]]-1))

dat.split[,27]
allele.list[[27]]
