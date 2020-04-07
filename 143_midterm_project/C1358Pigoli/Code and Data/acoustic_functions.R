library(abind)



sep.cov<-function(Data){ 
  N<-dim(Data)[1]
  d1 <- dim(Data)[2]
  d2 <- dim(Data)[3]
  
  # the two marginal covariance estimates
  C1<-matrix(0,dim(Data)[2],dim(Data)[2])
  C2<-matrix(0,dim(Data)[3],dim(Data)[3])
  
  M<-apply(Data,c(2,3),mean)
  
  ## remove the mean from the data
  for (k in 1:N){ 
    Data[k,,]<-Data[k,,]-M
  }
  
  ## compute marginal covariance estimate
  tmp.mat=matrix(c(aperm(Data, c(1,3,2))), ncol=d1)
  C1=cov(tmp.mat)*(N*d2-1)/N
  rm(tmp.mat)
  
  ## compute marginal covariance estimate
  ## this is faster
  tmp.mat=matrix(c(Data), ncol=d2)
  C2=cov(tmp.mat)*(N*d1-1)/N
  rm(tmp.mat)
  
  
  C1<-C1/sqrt(sum(diag(C1)))
  C2<-C2/sqrt(sum(diag(C2)))
  
  
  
  out<-list(C1,C2)
  out
}



sqrtM<-function(A){
  
  auto<-eigen(A)
  val<-pmax(auto$values,rep(0,dim(A)[1]))
  d<-sqrt(val)
  D<-diag(d)
  V<-auto$vectors
  L<-V%*%D%*%t(V)
  L
}

Procdist<-function(A,B){
  
  L1<-sqrtM(A)
  L2<-sqrtM(B)
  C<-t(L2)%*%L1
  SVD<-svd(C)
  R<-SVD$u%*%t(SVD$v)
  Obj<-L1-L2%*%R
  out<-sqrt(sum( rowSums( t(Obj) * t(Obj))))
  out
}


ProcMean<-function(G,toll=10^(-6),maxiter=10000){
  
  n<-dim(G)[1]
  p<-dim(G)[2]
  
  
  G0<-matrix(0,dim(G)[2],dim(G)[3])
  G1<-array(0,dim(G))
  
  for (k in 1:n){
    G1[k,,]<-sqrtM(G[k,,])
    G0<-G0+G1[k,,]
  }
  G0next<-G0/n
  iter<-0
  while(Procdist(G0,G0next)>toll & iter<maxiter){   
    G0<-G0next
    for (k in 1:n){
      C<-t(G1[k,,])%*%G0
      SVD<-svd(C)
      R<-SVD$u%*%t(SVD$v)
      G1[k,,]<-G1[k,,]%*%R
    }
    G0next<-matrix(0,dim(G1)[2],dim(G1)[3])
    for (k in 1:n){
      G0next<-G0next+G1[k,,]
    }
    G0next<-G0next/n
    iter<-iter+1
  }
  
  out<-G0next%*%t(G0next)
  
  out
}


Means_stat<-function(D,digit){
  
  I1<-which(digit==1)
  I2<-which(digit==2)
  I3<-which(digit==3)
  I4<-which(digit==4)
  I5<-which(digit==5)
  I6<-which(digit==6)
  I7<-which(digit==7)
  I8<-which(digit==8)
  I9<-which(digit==9)
  I10<-which(digit==10)
  
  M1<-apply(D[I1,,],c(2,3),mean)
  
  M2<-apply(D[I2,,],c(2,3),mean)
  M3<-apply(D[I3,,],c(2,3),mean)
  M4<-apply(D[I4,,],c(2,3),mean)
  M5<-apply(D[I5,,],c(2,3),mean)
  
  M6<-apply(D[I6,,],c(2,3),mean)
  M7<-apply(D[I7,,],c(2,3),mean)
  M8<-apply(D[I8,,],c(2,3),mean)
  M9<-apply(D[I9,,],c(2,3),mean)
  M10<-apply(D[I10,,],c(2,3),mean)

  M<-(1/10)*(M1+M2+M3+M4+M5+M6+M7+M8+M9+M10)
  
  f<-sum((M1-M)^2)+sum((M2-M)^2)+sum((M3-M)^2)+sum((M4-M)^2)+sum((M5-M)^2)+sum((M6-M)^2)+sum((M7-M)^2)+sum((M8-M)^2)+sum((M9-M)^2)+sum((M10-M)^2)
  
  return(f)
}





Means_test<-function(D1,d1,M){
  f0<-Means_stat(D1,d1)
  f<-c()
  for (k in 1:M){
    sub<-sample(d1,length(d1),replace=FALSE)
    f<-c(f,Means_stat(D1,sub))
  }
  pv<-length(which(f>f0))/M
  pv
}


Cov_stat<-function(D,digit,dir,dist=Procdist,mean=ProcMean){
  
  I1<-which(digit==1)
  I2<-which(digit==2)
  I3<-which(digit==3)
  I4<-which(digit==4)
  I5<-which(digit==5)
  I6<-which(digit==6)
  I7<-which(digit==7)
  I8<-which(digit==8)
  I9<-which(digit==9)
  I10<-which(digit==10)
  
  C1<-sep.cov(D[I1,,])[[dir]]
  
  C2<-sep.cov(D[I2,,])[[dir]]
  
  
  C3<-sep.cov(D[I3,,])[[dir]]
  
  
  C4<-sep.cov(D[I4,,])[[dir]]
  C5<-sep.cov(D[I5,,])[[dir]]
  C6<-sep.cov(D[I6,,])[[dir]]
  C7<-sep.cov(D[I7,,])[[dir]]
  C8<-sep.cov(D[I8,,])[[dir]]
  C9<-sep.cov(D[I9,,])[[dir]]
  C10<-sep.cov(D[I10,,])[[dir]]
  
  G<-abind(C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,along=0)
  C<-mean(G)
  f<-dist(C,C1)^2+dist(C,C2)^2+dist(C,C3)^2+dist(C,C4)^2+dist(C,C5)^2+dist(C,C6)^2+dist(C,C7)^2+dist(C,C8)^2+dist(C,C9)^2+dist(C,C10)^2
  
  f
}



Cov_test<-function(D1,d1,dir,M){
  f0<-Cov_stat(D1,d1,dir)
  f<-c()
  for (k in 1:M){
    sub<-sample(d1,length(d1),replace=FALSE)
    f<-c(f,Cov_stat(D1,sub,dir))
  }
  pv<-length(which(f>f0))/M
  pv
}

require(audio)


sound_from_spectrogram <- function(spec.logamplitude, spec.phase, rate, bit, wn)
{
  require(tuneR)
  A.grid<-sqrt(exp(spec.logamplitude/10))
  P<-spec.phase
  
  spec<-A.grid*(complex(real = cos(P), imaginary = sin(P)))
  
  x<-inverse_specgram(spec,wn)
  # return(audioSample(10*x,rate=rate))
  return(normalize(Wave(left=c(10*x), right=numeric(0), samp.rate=rate, bit=bit), as.character(bit)))
}


inverse_specgram<-function(stft, wn,h=ceiling(wn/2)){ 
  
  
  # estimate the length of the signal
  coln = dim(stft)[2]
  xlen = wn + (coln-1)*h
  x = rep(0, xlen)
  
  # form a periodic hamming window
  win = gausswin(wn)
  
  steps<-seq(0,h*(coln-1),by=h)
  # perform IFFT and weighted-OLA
  for (b in steps){ 
    # extract FFT points
    X <- stft[, 1 + b/h]
    X <- c(X, Conj(X[(length(X)):1]))
    
    # IFFT
    xprim = Re(fft(X,inverse=TRUE)/ length(X))
    
    # weighted-OLA
    x[(b+1):(b+wn)] = x[(b+1):(b+wn)] + (xprim*win)
  }
  
  
  W0 = sum(win^2)                   # find W0
  x = x*h/W0                        # scale the weighted-OLA
  
  x
}


