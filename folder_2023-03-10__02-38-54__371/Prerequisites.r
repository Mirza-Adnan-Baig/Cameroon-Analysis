

## To transform data



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



####



#---------------------------------------
# help functions 

varsets2 <- function(l){   #calculate all var sets
  # n number of loci
  # l number of alleles per locus
  n <- length(l)
  B <- array(1,c(prod(l),n))
  B[1:l[1],1] <- 1:l[1]
  lkmo <- l[1]
  if(n>1){
    for(k in 2:n){
      if(l[k]>1){
        lk <- lkmo*l[k]
        pick1 <- (lkmo+1):lk
        B[pick1,] <- B[rep(1:lkmo,l[k]-1),]
        B[pick1,k] <- rep(2:l[k],each=lkmo)
        lkmo <- lk
      }
                              
    }
  }
  B
}

varsets1 <- function(l){   #calculate all var sets
  # n number of loci
  # l number of alleles per locus
  n <- length(l)
  B <- array(0,c(prod(l),n))
  B[1:l[1],1] <- 0:(l[1]-1)
  lkmo <- l[1]
  if(n>1){
    for(k in 2:n){
      if(l[k]>1){
        lk <- lkmo*l[k]
        pick1 <- (lkmo+1):lk
        B[pick1,] <- B[rep(1:lkmo,l[k]-1),]
        B[pick1,k] <- rep(1:(l[k]-1),each=lkmo)
        lkmo <- lk   
      }
                           
    }
  }
  B
}

varsets <- function(l,n){   #calculate all var sets
  # n number of loci
  # l number of alleles per locus
  B <- array(0,c(l^n,n))
  B[1:l,1] <- 0:(l-1)
  lkmo <- l
  if(n>1){
    for(k in 2:n){
      lk <- lkmo*l
      pick1 <- (lkmo+1):lk
      B[pick1,] <- B[rep(1:lkmo,l-1),]
      B[pick1,k] <- rep(1:(l-1),each=lkmo)
      lkmo <- lk                        
    }
  }
  B
}

gead <- function(x,l,n){   ## calculates geadic expression of each element of vectorx 
  
  l <- rep(l,n)

  out <- array(0,c(length(x),n))
  div <- c(1,cumprod(l[1:(n-1)]))
  for(k in n:1){
    r <- x%%div[k]
    #print(c(r,div[k]))
    out[,k] <- (x-r)/div[k]
    x <- r
  }
  out
}

gead1 <- function(x,l){   ## calculates general geadic expression of each element of vector x
  n <- length(l)
  out <- array(0,c(length(x),n))
  div <- c(1,cumprod(l[1:(n-1)]))
  for(k in n:1){
    r <- x%%div[k]
    out[,k] <- (x-r)/div[k]
    x <- r
  }
  out
}



#_____________________________
# main function
est <- function(X,Nx,l){
  eps <- 10^-8
  N <- sum(Nx)
  nn <- nrow(X) 
  n <- ncol(X)
  x <- X
  
  #_____________________________
  # this calculates the list to pick proper haplotype freuencies, i.e. sets Ay
  
  #allconf <- function(x,l,n){  # x array
  hapll <- list()
  if(length(l)==1){
    l <- rep(l,n)
  }else{
    n <- length(l)
  }
  ggead <- c(1,cumprod(l[1:(n-1)]))
  Hx <- list()
  Ax <- list()
  alnum <- list()
  bin2num <- list()
  for(k in 1:n){
    alnum[[k]] <- 1:l[[k]]
    bin2num[[k]] <- 2^(0:(l[[k]]-1))
  }
  alcnt <- array(0,n)
  for(u in 1:nn){
    Hx[[u]] <- list(array(0,n),list(),list(),list(),list(),array(0,n))
    Ax[[u]] <- list(list(),list(),list(),list())
    for(k in 1:n){
      temp <- gead(x[u,k],2,l[[k]])
      temp1 <- varsets1(temp+1)[-1,]
      Hx[[u]][[1]][k] <- sum(temp)
      Hx[[u]][[2]][[k]] <- (alnum[[k]])[temp*alnum[[k]]]
      Hx[[u]][[3]][[k]] <- temp1
      Hx[[u]][[4]][[k]] <- temp1%*%bin2num[[k]] # subsets of alleles
      Hx[[u]][[5]][[k]] <- Hx[[u]][[3]][[k]]%*%rep(1,l[k]) # number of alleles in subset at locus k
      Hx[[u]][[6]][k] <- length(temp1%*%bin2num[[k]])
    }
    vz1 <- sum(Hx[[u]][[1]])
    temp2 <- prod(Hx[[u]][[6]])
    Ax[[u]][[1]] <- temp2
    Ax[[u]][[2]] <- varsets2(Hx[[u]][[6]])
    for(k in 1:n){
      Ax[[u]][[2]][,k] <- Hx[[u]][[4]][[k]][Ax[[u]][[2]][,k]]
    }
    for(j in 1:temp2){
      Ax[[u]][[3]][[j]] <- list() 
      for(k in 1:n){
        temp <- gead(Ax[[u]][[2]][j,k],2,l[[k]])
        temp1 <- (alnum[[k]])[temp*alnum[[k]]]
        alcnt[k] <- length(temp1)
        Ax[[u]][[3]][[j]][[k]] <- temp1
      }
      Ax[[u]][[4]][[j]] <- varsets2(alcnt)
      for(k in 1:n){
        Ax[[u]][[4]][[j]][,k] <- Ax[[u]][[3]][[j]][[k]][Ax[[u]][[4]][[j]][,k]]
      }    
      Ax[[u]][[4]][[j]] <- as.character((Ax[[u]][[4]][[j]]-1)%*%ggead+1)
      Ax[[u]][[3]][[j]] <- (-1)^(vz1+sum(alcnt))
    }
    hapll[[u]] <- Ax[[u]][[4]][[temp2]]
    
  }  
  
  hapl1 <- unique(unlist(hapll))
  
  #---------------------------------------
  # initialize parameters
  H <- length(hapl1)
  pp <- array(rep(1/H,H),c(H,1))
  rownames(pp) <- hapl1
  
  #initial list#
  num0 <- pp*0
  cond1 <- 1  ## condition to stop EM alg! 
  la <- 2
  num <- num0
  rownames(num) <- hapl1
  Bcoeff <- num0
  rownames(Bcoeff) <- hapl1
  
  while(cond1>eps){
    Ccoeff <- 0
    Bcoeff <- num0 #reset B coefficients to 0 in next iteration
    num <- num0  #reset numerator to 0 in next iteration
    for(u in 1:nn){
      denom <- 0
      num <- num0
      CC <- 0
      for(k in 1:Ax[[u]][[1]]){
        p <- sum(pp[Ax[[u]][[4]][[k]],])
        vz <- Ax[[u]][[3]][[k]]
        lap <- la*p
        exlap <- vz*exp(lap)
        denom <- denom + exlap-vz  ##   = (1-)^(Nx-Ny)*(Exp(lambda*sum p)-1) = (1-)^(Nx-Ny)*G(sum p)
        num[Ax[[u]][[4]][[k]],] <- num[Ax[[u]][[4]][[k]],]+ exlap#*pp[Ax[[u]][[1]][[k]],]
        ## exlap =  (1-)^(Nx-Ny) G'(sum p)   --- denominator of generating functions cancels out!
        CC <- CC + exlap*p
      }
      num <- num*pp
      denom <- Nx[u]/denom
      denom <- la*denom
      Ccoeff <- Ccoeff + CC*denom
      Bcoeff <- Bcoeff + num*denom
    }
    Ccoeff <- Ccoeff/N
    ppn <- Bcoeff/(sum(Bcoeff))
    
    ### Newton step
    cond2 <- 1
    xt <- Ccoeff   ### good initial condition
    while(cond2 > eps){
      ex <- exp(-xt)
      xtn <- xt + (1-ex)*(xt + Ccoeff*ex - Ccoeff)/(ex*xt+ex-1)
      cond2 <- abs(xtn-xt)
      xt <- xtn
    }
    cond1 <- abs(xt-la)+sqrt(sum((pp-ppn)^2))
    la <- xt
    pp <- ppn
    #print("________")
    #print(cond1)
    #print(c(la,pp))
    
  }
  list(la,pp)
}

#_____________________________
# likelihood function
likeGen <- function(pp,la,Nx,N,Ax){
  logli <- 0
  num0 <- 0
  Bcoeff <- num0 #reset B coefficients to 0 in next iteration
  num <- num0  #reset numerator to 0 in next iteration
  nn <- (length(Ax))
  for(u in 1:nn){
    denom <- 0
    num <- num0
    CC <- 0
    for(k in 1:Ax[[u]][[1]]){
      p <- sum(pp[Ax[[u]][[4]][[k]],])
      vz <- Ax[[u]][[3]][[k]]
      lap <- la*p
      exlap <- vz*exp(lap)
      denom <- denom + (exlap-vz) 
    }
    
    logli <-  logli+ Nx[u]*log(denom/(exp(la)-1))
  }
 logli 
  
    
}  


#_____________________________
# labels output frequencies correctly
hapl.names <- function(pp,allist){ ## allist list with alleles per locos
  n <- length(allist)
  allnum <- unlist(lapply(allist,length))
  hapl <- gead1(as.numeric(rownames(pp))-1,allnum)+1
  hapl1 <- array(,dim(hapl))
  for(l in 1: ncol(hapl)){
    for(m in 1: nrow(hapl)){
      hapl1[m,l] <- allist[[l]][hapl[m,l]]
    }
  }
  hapl1 <- apply(hapl1,1,function(x) paste(x,sep="",collapse="-"))
  rownames(pp) <- hapl1
  pp
}
#_____________________________
# example data
l <- c(3,4,2)
X <- array(c(3,1,7,2,2,1),c(2,3))
Nx <- c(2,3)
est(X,Nx,l)

(1,1,0)
(1,1,1,0)
(0,1,0,0)
