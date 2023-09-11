library('xlsx')
library('ggplot2')
library('stringr')
library('Rmpfr')

#************************************************************************************
#This function imports data
#************************************************************************************

DatImp <- function(path){
  if(substring(path,nchar(path)-3,nchar(path))==".xls"){
    dat <- read.xlsx(path,1)
  }
  else{
    if(substring(path,nchar(path)-4,nchar(path))==".xlsx"){
      dat <- read.xlsx(path,1)
    }
    else{
      if(substring(path,nchar(path)-3,nchar(path))==".txt"){
        dat <- read.table(path,header=TRUE, sep="\t")
      }
      else{
        if(substring(path,nchar(path)-3,nchar(path))==".csv"){
          dat <- read.csv(path,header=TRUE,sep=";")
        }
      }
    }  
  }
  dat
}  

#************************************************************************************
#This function calculate Nk
#************************************************************************************    

Nk <- function(dat){
  for(k in 1:nrow(dat)){
    if(dat[k,1]==""||is.na(dat[k,1])){
      dat[k,1] <- dat[k-1,1]
    }
  }
  dat <- dat[!is.na(dat[,2]),]
  N <- length(unique(dat[,1]))
  Nknum <- length(unique(dat[,2]))
  dat <- dat[!duplicated(dat),]
  list(N,t(as.matrix(table(dat[,2]))))
}

#************************************************************************************
#This function calculates the ML estimate
#************************************************************************************

MLE<-function(Nk,N){
  out <- list(NA, NA,NA,NA,NA)
  if(sum(Nk)>N && min(N-Nk)>0 && N>0){
    sel <- Nk
    Nk <- sel[sel>0]
    nk <- Nk/N
    l1 <- 2.5         # initial value
    la <- 2.5
    l0 <- 0
    eps <- 10^(-8)       # precision 
    out <- list(NA, NA,NA,NA,NA)
    k <- 1
    while(abs(l0-l1)>eps && k<50 && l1>0){
      k <- k+1
      l0 <- l1
      l1 <- l0-(l0+sum(log(1-nk*(1-exp(-l0)))))/(1-sum(nk/(exp(l0)*(1-nk)+nk)))
    }
    if(k==50 || l1<0){
      # print(c(l0,l1,Nk))
      for(st in 1:10){
        #print(st)
        l1 <- st
        l0 <- l1+1
        k <- 1
        while(abs(l0-l1)>eps && k<100 && l1>0){
          k <- k+1
          l0 <- l1
          l1 <- l0-(l0+sum(log(1-nk*(1-exp(-l0)))))/(1-sum(nk/(exp(l0)*(1-nk)+nk)))
        }
        if(abs(l0-l1)<eps){
          break
        }
      }
      if(abs(l0-l1)>eps){      # if numerical problems occur, calculations are performed with higher precision
        print(l0-l1)
        l1 <- mpfr(10*la,precBits=100)
        l0 <- l1+1
        eps <- mpfr(eps,precBits=100)
        while(abs(l0-l1)>eps){
          l0 <- l1
          l1=l0-(l0+sum(log(1-nk*(1-exp(-l0)))))/(1-sum(nk/(exp(l0)*(1-nk)+nk)))
          print(l1)
        }
      }        
    }
    pk <- -1/l1*log(1-nk*(1-exp(-l1)))   
    ml <- (-N)*log(exp(l1)-1)+sum(Nk*log(exp(l1*pk)-1))	 #maximum log-likelihood 
    pk1 <- array(0,length(sel))  
    pk1[sel>0] <- pk    
    out <- list(ml,l1,pk1,N,sel)
  }else{
    if(N==0){
      warning("Sample size N needs to be positive")
    }else{
      if(min(N-Nk)==0){
        warning("Nk=N for one k, the MLE does not exist")
      }else{
        if(sum(Nk)==N){
          warning("Only one lineage per sample, the MLE is degenerate")
          sel <- Nk
          Nk <- sel[sel>0]
          pk <- Nk/N
          ml <- sum(Nk*log(Nk))
          pk1 <- array(0,length(sel)) 
          pk1[sel>0] <- pk    
          l1 <- 0
          out <- list(ml,l1,pk1,N,sel)
        }
        
      }
    }
  }
  out	
}

#************************************************************************************
#This function provides conditinally Poisson distributed numbers
#************************************************************************************

cpoiss<-function(lambda,n){
  # lambda...Poisson parameter, n...number of conditional Poisson numbers to be generated
  
  m <- 100 # to accelerate computation it is assumed that m<100 is generically drawn
  out <- rep(0,n)
  x <- runif(n,min=0,max=1)
  p0 <- ppois(0,lambda)
  nc <- 1/(1-exp(-lambda))
  pvec <- (ppois(1:m,lambda)-p0)*nc
  pvec <- c(pvec,1) 
  for (i in 1:n){
    k <- 1
    while(x[i] > pvec[k]){
      k <- k+1
    }
    if(k==m){ # if a m>=100 is drawn this is executed
      k <- k+1
      a <- dpois(k,lambda)*nc
      b <- pvec[m]+a
      while(x[i]>b){
        k <- k+1
        a <- a*lambda/k
        b <- b+a
      }
    }
    out[i] <- k
  }
  out
}

#************************************************************************************
#This function returns multinomial samples
#************************************************************************************

mnom<-function(M,dist){#M=(m_1,...,m_N), dist=(p_1,...,p_n)...frequency distribution  
  # function return N multinomial vectors, where the kth vector is distributiod Mult(m_k;p_1,...,p_n) 
  N <-length(M)
  out<-array(0,dim=c(N,length(dist)))
  for(k in 1:N){
    out[k,]=rmultinom(1,M[k],dist)
  }
  out
}

#************************************************************************************
#This function calculates differences of frequency estimates and true parameters 
#************************************************************************************ 

pstat <- function(p,logp,hp){
	dp <- p-hp
	adp <- abs(dp)
	n1 <- sum(adp)
	n2 <- sqrt(sum(adp^2))
	ninf <- max(adp)
	sp <- hp+p
	JSdiv <- 2*log(2)+hp%*%log(hp/sp)+p%*%log(p/sp)
	Hdiv <- sqrt(sum((sqrt(p)-sqrt(hp))^2)/2)
	Bdiv <- -log(sum(sqrt(p*hp)))
	c(n1,n2,ninf,JSdiv,Hdiv,Bdiv)	
}

#************************************************************************************
#This is the important loop
#************************************************************************************ 

run1 <- function(lambda,p,N,K,path){
  logp <- log(p)
  n <- length(p)
  out <- array(0,c(K,7+length(p))) 
  out1 <- array(0,c(K,7+length(p))) 
  kk <- 0
  k <- 0
  while(k<K){
    kk <- kk+1
    m <- cpoiss(lambda,N)
    Nk <- colSums(sign(mnom(m,p)))
    if(sum(Nk)>N && min(N-Nk)>0){
      k <- k+1
      data <- MLE(Nk,N)
      out[k,] <- c(data[[2]],data[[3]],pstat(p,logp,data[[3]]))
    }
  }
  la <- out[,1] 
  out[,1]  <- la/(1-exp(-la))
  a <- rbind(apply(out,2,function(x) as.numeric(summary(x[!is.na(x)],digits=7))),apply(out,2,function(x) var(x[!is.na(x)])))
  out1 <- t(c(c(lambda,p,N,K,kk),a))
  write.table(round(out1,digits=7), paste(path,toString(n),".txt",sep=""), sep="\t", append = TRUE, col.names=FALSE, row.names = FALSE)
  
}

#************************************************************************************
#This function runs the simulations
#************************************************************************************ 

runsim <- function(lamin,lamax,lainc,pmat,Nvec,K=10000,path){
	for(k in 1: dim(pmat)[1]){
		p <- pmat[k,]
		#print(p)
		for(n in 1:length(Nvec)){
			N <- Nvec[n]
			#print(N)
			for(la in seq(lamin,lamax,lainc)){
        #print(la)
				run1(la,p,N,K,path)
			}
		}
	}
}

#************************************************************************************
#This function plots output
#************************************************************************************ 

Glpe <- function(inpath,n,quant,sta,outpath){
  warn <- 0
  if(!is.element(sta,c("min","Q1","me","E","Q2","max","var","CV"))){
    warn <- 1
    warning("4th input argument must be one of the following: min, Q1, me, E, Q2, max, var, CV" )
  }
  if(!is.element(quant,c("MLE",paste(c("p"), 1:n, sep=""),"1-norm", "2-norm", "sup-norm", "JS-div", "H-div", "B-div" ))){
    warn <- 1
    warning("3rd input argument must be one of the following: MLE, p1, ... , pn, 1-norm, 2-norm, sup-norm, JS-div, H-div, B-div" )
  }
  if(sta=="CV"){
    if(!is.element(quant,c("MLE",paste(c("p"), 1:n, sep="")))){
      warn <- 1
      warning("3rd input argument must be one of the following: MLE, p1, ... , pn" )
    }
  }
  dat <- read.table(inpath, sep="\t")
  data <- label(dat,n)
  pmat <- unique(data[,2:(n+1)])
  K <- dim(pmat)[1]
  pmat <- t(pmat)
  attributes(pmat) <- NULL
  pmat <- array(pmat,c(n,K))
  qu <- str_replace_all(quant,pattern=" ",repl="_")
  a <- 1:n
  data1 <- data
  data1[,"lam"] <- data1[,"lam"]/(1-exp(-data1[,"lam"]))
  xax <- expression(frac(lambda,1-e^- lambda))#label for x axis
  yax <- paste(sta,quant,sep=" ")
  if(sta=="E"){
    aaa <- quant
    if(quant=="JS-div"){
      #print(data1[,"E JS-div"])
      aaa <- "Jensen-Shannon div."
    }
    if(quant=="H-div"){
      #print(data1[,"E JS-div"])
      aaa <- "Hellinger distance"
    }
    yax <- paste("mean",aaa,sep=" ")
  }
  if(sta=="E" && quant=="MLE"){
    data1[,"E MLE"] <- (data1[,"E MLE"]/data1[,"lam"]-1)*100
    yax <- "bias in %"
  }
  if(sta=="me"&& quant=="MLE"){
    data1[,"me MLE"] <- (data1[,"me MLE"]/data1[,"lam"]-1)*100
    yax <- "median bias in %"
  }
  if(sta=="CV"){
    yax <- "CV in %"
    pp <- paste("p",1:n,sep="")
    ppCR <- paste("CR p",1:n,sep="")
    el <- exp(data[,"lam"])
    pmat1 <- data[,pp]
    elpk <- exp(data[,"lam"]*pmat1)
    A <- rowSums(elpk-1)
    CRlam <- el*(el-1-data[,"lam"])^2/(data[,"N"]*(el-1)^2)*A/(el-1-A)
    CRpp <- (el-1)^2/(data[,"N"]*el*data[,"lam"]^2)*((elpk-1)/(el-1)+(pmat1^2*A-2*pmat1*(elpk-1)+(elpk-1)^2/(el-1))/(el-1-A))
    #print(CRpp)
    for(j in 1:n){
      data1[,ppCR[j]] <- sqrt(CRpp[,j])/pmat1[,j] *100
    }
    data1[,"CR MLE"] <- sqrt(CRlam)*(1-1/el)/data[,"lam"]*100
    
    data1[,paste("CV",quant,sep=" ")] <- sqrt(data1[,paste("var",quant,sep=" ")])/data1[,paste("E",quant,sep=" ")]*100
  }
  for(k in 1:K){
    p <- pmat[,k]
    dattemp1 <- data1[colSums(abs(t(data[,2:(n+1)])-p))==0,]
    title <- paste(",\", p\"[",a,"],\"=",p,"\"",sep="",collapse="")
    title <- paste("paste(\"",substring(title,4,nchar(title)),")",collapse="")
    if(n==10){
      pnew <- array("",c(1,n))
      temp <- unique(p)
      for(l in 1:length(temp)){
        pnew[max((1:n)[p==temp[l]])] <- paste(toString(temp[l]),", ",sep="",collapse="")
        title <- paste(",\" p\"[",a,"],\"=",pnew,"\"",collapse="",sep="")
        title <- paste("paste(\"",substring(title,4,nchar(title)-3),"\")",collapse="")
      }
    }
    if(warn==0){
    out <- paste(outpath,sta,"_",qu,toString(n),"_",toString(k),".pdf",sep="")
    pdf(out,height=5)
    print(Lpe(dattemp1,quant,sta,"N",title,xax,yax))
    dev.off()
    }
  }
}

Lpe <- function(data,quant,quant1,gp,title,xax,yax){
  cbPalette <- c("#000000", "#1B9E77",	"#D95F02",	"#7570B3",	"#E7298A",	"#66A61E",	"#E6AB02",	"#A6761D",	"#666666", "#1f78b4", "#33a02c", "#999999", "#332288")	
  E <- paste(quant1,quant,sep=" ")
  data1 <- data[,c(gp,"lam",gp,E)]
  data1[,3] <- data1[,4]
  colnames(data1) <- c("N","lam","bias","mean")
  data1 <- data.frame(N=as.factor(data1[,1]),lam=data1[,2],bias=data1[,3],mean=data1[,4])
  yll <- min(-0.01,round(min(data1[,"bias"]),digits=2)-0.01)
  yul <- max(0.04,round(max(data1[,"bias"]),digits=2)+0.01)
  if(quant1=="CV"){
    data1 <- data[,c(gp,"lam",gp,E,paste("CR",quant,sep=" "))]
    data1[,3] <- data1[,4]
    colnames(data1) <- c("N","lam","bias","mean","CR")
    data1 <- data.frame(N=as.factor(data1[,1]),lam=data1[,2],bias=data1[,3],mean=data1[,4],CR=data1[,5])
    yll <- min(-0.01,round(min(c(data1[,"bias"],data1[,5])),digits=2)-0.01)
    yul <- max(0.04,round(max(c(data1[,"bias"],data1[,5])),digits=2)+0.01)
  }
  
  
  p <- ggplot(data=data1,aes(x=lam,y=bias,group=N,colour=N))
  p <- p + geom_line()
  if(quant1=="CV"){
    p <- p + geom_line(data=data1,aes(x=lam,y=CR,group=N,colour=N),linetype=5)
  }
  p <- p + theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
  p <- p + theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank())
  p <- p + theme(panel.background = element_rect(colour='black',fill='white'))
  p <- p + theme(legend.key=element_rect(colour='white', fill='white'))
  p <- p + theme(axis.text = element_text(colour='black',size = rel(1.4)),title = element_text(size = rel(1.3)))
  p <- p + theme(axis.ticks = element_line(color = "black"))
  p <- p + theme(axis.title = element_text(size = rel(1.3)))
  p <- p + theme(legend.text = element_text(size = rel(1.3)))
  p <- p + theme(legend.title = element_text(size = rel(1.1)))
  p <- p + ylim(yll,yul)
  p <- p + theme(plot.title = element_text(hjust = 0.5))
  p <- p + labs(title=parse(text=title))
  p <- p + labs(x=xax,y=yax)
  p <- p + scale_colour_manual(values=cbPalette)
  p
  
}

label <- function(dat,n){
  cn<- c("MLE",paste(c("p"), 1:n, sep=""),"1-norm", "2-norm", "sup-norm", "JS-div", "H-div", "B-div" )
  dn <- c("min","Q1","me","E","Q2","max","var")
  lab <- paste(dn,rep(cn,each=7),sep=" ")
  en <- c('lam',paste(c("p"), 1:n, sep=""),'N','K','K1')
  lab <- c(en,lab)
  colnames(dat) <- lab
  dat
} 

#############################################
#############################################
#END
#############################################
#############################################