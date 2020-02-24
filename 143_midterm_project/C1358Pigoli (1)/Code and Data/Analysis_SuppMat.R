library(R.matlab)
library(fields)
library(rgl)
library(fda)
library(signal)

#############################################################################
########   This script implements the analysis of acoustic phonetic data ####
########           described in Pigoli et al., 2015                      ####
########                                                                 ####
#############################################################################

source('../Code and Data/acoustic_functions.R')


########## Preprocessing ########

Sounds<-readMat('../Code and Data/PrunedLanguageDataset.mat')



wn<-160
rate<-16000
#plot(Sounds[[136]],type='l',main=names(Sounds)[136],col=4)
plot(Sounds[[100]],type='l',main=names(Sounds)[100],col=4)


## extract sounds and compute spectrograms (takes some time) 

language<-rep(NA, length(Sounds))
digit<-rep(NA, length(Sounds))
gender<-rep(NA, length(Sounds))
sounds.spec <- vector("list",length(Sounds))
for (k in 1:length(Sounds)){ 
  S<-specgram(Sounds[[k]], n=wn, Fs=16000, window=gausswin(wn))
  A<-abs(S$S)^2
  A<-A/max(A)
  A = pmax(A, 10^(-16))
  sounds.spec[[k]]<-list(power=10*log(A), phase=Arg(S$S))
  header<-unlist(strsplit(names(Sounds)[k], "[.]"))
  language[k]<-header[1]
  digit[k]<-as.double(header[3])
  gender[k]<-header[4]
}
Speaker_info<-data.frame(language=language,digit=digit,gender=gender)


nsounds <- length(sounds.spec)
nfreq<-dim(sounds.spec[[1]]$power)[1]

######## interpolate on common grid ##########

eq.grid<-seq(0,1,len=100)

# Data.grid<-array(NA, c(nsounds, nfreq, 100))
# for (k in 1:nsounds){
#   print(paste0(k, "/", nsounds))
#   for (i in 1:nfreq){ 
#     Data.grid[k,i,]<-approx(seq(0,1,len=dim(sounds.spec[[k]]$power)[2]),sounds.spec[[k]]$power[i,],eq.grid)$y
#   }
# }

# example 
# example<-smooth.surface(sounds.spec[[147]]$power,nfreq,100,b=20)
# 
# brk<-seq(-250,0,len=100)
# color = tim.colors(length(brk)-1)
# image.plot(sounds.spec[[147]]$power,breaks = brk,col=color)
# image.plot(example,breaks = brk,col=color)


Data.grid<-array(NA, c(nsounds, nfreq, 100))
for  (k in 1:nsounds){
  print(paste0(k, "/", nsounds))
  Data.grid[k,,]<-smooth.surface(sounds.spec[[k]]$power,nfreq,100,b=20)
}

save(Speaker_info,Data.grid,file="SmoothPDS.RData")

load("SmoothPDS.RData")
attach(Speaker_info)

Data2<-array(NA,dim(Data.grid))
for(num in 1:10){
  cat(paste0("\n\n\n", 
             "----------------------------------------",
             "\n-------- ", num, "--------\n",
             "----------------------------------------",
             "\n\n\n") )
  subset.i <- which( digit == num )
  local.data.grid <- Data.grid[subset.i,,]
  
  library(fdasrvf)
  xseq<-seq(0,1,len=dim(local.data.grid)[3])
  
  npc=3
  #FV<-align_fPCA(t(local.data.grid[,1,]),xseq)
  FV<-align_fPCA(t(local.data.grid[,5,]),xseq, num=npc, parallel=TRUE, cores=2, showplot=FALSE); ## align according to frequency 5, where the peak is most often
  
  
  local.aligned.sounds<-array(NA,dim(local.data.grid))
  
  for (i in 1:dim(local.aligned.sounds)[1]){
    for (j in 1:dim(local.aligned.sounds)[2]){
      local.aligned.sounds[i,j,]<-approx(xseq,local.data.grid[i,j,],xout=FV$gam[,i],rule=2)$y
    }
  }
  
  Data2[subset.i,,]<-local.aligned.sounds
}

save(Speaker_info,Data2,file="SVRF_WarpedPSD.RData")

load("SVRF_WarpedPSD_SuppMat.RData")
attach(Speaker_info)


# Estimation of word means

MFr<-array(NA,dim(Data2))
MI<-array(NA,dim(Data2))
MP<-array(NA,dim(Data2))
MSA<-array(NA,dim(Data2))
MSI<-array(NA,dim(Data2))


for ( i in 1:10){ 
  
  Fr1<-which(language=='FR'& digit==i)
  I1<-which(language=='IT'& digit==i)
  P1<-which(language=='PO'& digit==i)
  SA1<-which(language=='SA'& digit==i)
  SI1<-which(language=='SI'& digit==i)
  
  
  
  MFr[i,,]<-apply(Data2[Fr1,,],c(2,3),mean)
  
  MI[i,,]<-apply(Data2[I1,,],c(2,3),mean)
  
  
  MP[i,,]<-apply(Data2[P1,,],c(2,3),mean)
  
  MSA[i,,]<-apply(Data2[SA1,,],c(2,3),mean)
  
  MSI[i,,]<-apply(Data2[SI1,,],c(2,3),mean)
  
}


### Estimation of the separable language covariances

D<-array(0,dim(Data2))

for ( i in 1:10){ 
  
  Fr1<-which(language=='FR'& digit==i)
  I1<-which(language=='IT'& digit==i)
  P1<-which(language=='PO'& digit==i)
  SA1<-which(language=='SA'& digit==i)
  SI1<-which(language=='SI'& digit==i)
  
  
  
  M<-apply(Data2[Fr1,,],c(2,3),mean)
  for (k in Fr1){
    D[k,,]<-Data2[k,,]-M
  }
  
  M<-apply(Data2[I1,,],c(2,3),mean)
  for (k in I1){
    D[k,,]<-Data2[k,,]-M
  }
  
  M<-apply(Data2[P1,,],c(2,3),mean)
  for (k in P1){
    D[k,,]<-Data2[k,,]-M
  }
  
  M<-apply(Data2[SA1,,],c(2,3),mean)
  for (k in SA1){
    D[k,,]<-Data2[k,,]-M
  }
  M<-apply(Data2[SI1,,],c(2,3),mean)
  for (k in SI1){
    D[k,,]<-Data2[k,,]-M
  }
  
}


Fr<-which(language=='FR')
I<-which(language=='IT')
P<-which(language=='PO')
SA<-which(language=='SA')
SI<-which(language=='SI')


CFT<-sep.cov(D[Fr,,])[[1]]

CPT<-sep.cov(D[P,,])[[1]]


CIT<-sep.cov(D[I,,])[[1]]


CSAT<-sep.cov(D[SA,,])[[1]]

CSIT<-sep.cov(D[SI,,])[[1]]

CFW<-sep.cov(D[Fr,,])[[2]]

CPW<-sep.cov(D[P,,])[[2]]



CIW<-sep.cov(D[I,,])[[2]]


CSAW<-sep.cov(D[SA,,])[[2]]

CSIW<-sep.cov(D[SI,,])[[2]]



#########   Permutations tests for the means ########### 

lg<-unique(language)

# P-vlues of the permutation tests for the equality of the means across the words
# for each language

pv_means<-rep(NA,length(lg))

for (k in 1:length(lg)){ 
  lang<-which(language==lg[k])  
  pv_means[k]<-Means_test(Data2[lang,,],digit[lang],M=1000)
}

pv_means

# PO IT SI FR SA
# 0.962 0.014 0.002 0.000 0.000



#########   Permutations tests for the covariance ########### 

lg<-unique(language)

# Centering with respect to the word means

D<-array(NA,dim(Data2))

for ( i in 1:10){ 
  
  d1<-which(digit==i)
  
  M<-apply(Data2[d1,,],c(2,3),mean)
  for (k in d1){
    D[k,,]<-Data2[k,,]-M
  }
  
  
}




# P-vlues of the permutation tests for the equality of the frequency covariance across the words
# for each language



pv_freq<-rep(NA,length(lg))

for (k in 1:length(lg)){ 
  lang<-which(language==lg[k])  
  pv_freq[k]<-Cov_test(D[lang,,],digit[lang],dir=1,M=100)
}

pv_freq

# PO IT SI FR SA
# 0.95 1.00 1.00 0.38 0.83

# P-values of the permutation tests for the equality of the time covariance across the words
# for each language


pv_time<-rep(NA,length(lg))

for (k in 1:length(lg)){ 
  lang<-which(language==lg[k])  
  pv_time[k]<-Cov_test(D[lang,,],digit[lang],dir=2,M=100)
}

pv_time

# 0.74 1.00 0.20 0.91 0.71



###### figures #######

P1<-which(language=='FR'& digit==1)

GP<-Data2[P1[1],,]
#GP<-MFr[1,,]
brk<-seq(min(GP),max(GP),len=60)
#brk<-seq(-0.01,0.05,len=60)
color = tim.colors(length(brk)-1)

#M<-par3d("userMatrix")
#rgl.viewpoint(userMatrix=M) 
par3d(cex=2)
zcol  = cut(GP, brk)
time<-seq(0,1,len=dim(GP)[2])
freq<-seq(0,8,len=80)
persp3d(freq,time,GP,xlab="",ylab="",zlab="",cex.lab=1.5,lwd=2,cex.axis=1.5,col=color[zcol])

#rgl.postscript("spectro1.pdf", fmt="pdf", drawText=TRUE )
rgl.postscript("spectroF1_SM.eps", fmt="eps", drawText=TRUE )

#rgl.snapshot("spectroF1_R.png", fmt="png" )

image.plot(time,freq,t(GP),main="French speaker, word 'un'",xlab="Time (standardized)",ylab="Frequency (kHz)",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)

par(mar=c(5,6,4,2))
image.plot(freq,freq,t(CFT),main="French",xlab="Frequency (kHz)",ylab="Frequency (kHz)",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
image.plot(freq,freq,t(CPT),main="Portuguese",xlab="Frequency (kHz)",ylab="Frequency (kHz)",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
image.plot(freq,freq,t(CIT),main="Italian",xlab="Frequency (kHz)",ylab="Frequency (kHz)",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
image.plot(freq,freq,t(CSAT),main="American Spanish",xlab="Frequency (kHz)",ylab="Frequency (kHz)",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
image.plot(freq,freq,t(CSIT),main="Castilian Spanish",xlab="Frequency (kHz)",ylab="Frequency (kHz)",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)

image.plot(time,time,t(CFW),main="French",xlab="Time",ylab="Time",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
image.plot(time,time,t(CPW),main="Portuguese",xlab="Time",ylab="Time",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
image.plot(time,time,t(CIW),main="Italian",xlab="Time",ylab="Time",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
image.plot(time,time,t(CSAW),main="American Spanish",xlab="Time",ylab="Time",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
image.plot(time,time,t(CSIW),main="Castilian Spanish",xlab="Time",ylab="Time",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)


#################################################################################################



##### Transformation of a speaker log-spectrograms from French to Portuguese #############

l1<-"FR"
l2<-"PO"

# Select the parameters of the language

C1<-CFT
C2<-CPT
T1<-CFW
T2<-CPW


# Choice of the word

dig<-1

# Computation of the means

L1<-which(language==l1 & digit==dig)
L2<-which(language==l2 & digit==dig)


M1<-apply(Data2[L1,,],c(2,3),mean)
M2<-apply(Data2[L2,,],c(2,3),mean)

# Residual log-spectrogram 

Res<-Data2[L1[1],,]-M1

# Transformation:

PS<-(sqrtM(C2)%*%solve(sqrtM(C1)))%*%Res%*%(solve(sqrtM(T1))%*%sqrtM(T2))+M2

# save it in a matlab file for sound reconstruction
#writeMat("FR2I_proj.mat",F2P=PS)


image.plot(time,freq,t(PS),main="",xlab="Time (standardized)",ylab="Frequency (kHz)",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)

image.plot(time,freq,t(Data2[L2[1],,]),main="",xlab="Time (standardized)",ylab="Frequency (kHz)",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)


############ Interpolation and extrapolation of sound surfaces ################

l1<-"FR"
l2<-"PO"

# Select the parameters of the language

C1<-CFT
C2<-CPT
T1<-CFW
T2<-CPW

# Choice of the word

dig<-5

# Computation of the means

L1<-which(language==l1 & digit==dig)
L2<-which(language==l2 & digit==dig)

M1<-apply(Data2[L1,,],c(2,3),mean)
M2<-apply(Data2[L2,,],c(2,3),mean)


# Grid for the interpolation (values between 0 and 1) or extrapolation (<0 or >1)

weights<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)


######### Interpolation/extrapolation between two (different) recorded speakers ########

# Speakers

spk1<-1
spk2<-1

Res1<-Data2[L1[spk1],,]-M1
Res2<-Data2[L2[spk2],,]-M1


Spect_path<-list()

for ( h in 1:length(weights)){
  
  w<-weights[h]
  
  
  ############## Frequency covariance ###################
  
  SI<-sqrtM(C1)
  SA<-sqrtM(C2)
  C<-t(SA)%*%SI
  SVD<-svd(C)
  R<-SVD$u%*%t(SVD$v)
  C3<-(SI+w*(SA-SI%*%R))%*%t(SI+w*(SA-SI%*%R))
  
  ############## Time covariance ###################
  
  SI<-sqrtM(T1)
  SA<-sqrtM(T2)
  C<-t(SA)%*%SI
  SVD<-svd(C)
  R<-SVD$u%*%t(SVD$v)
  T3<-(SI+w*(SA-SI%*%R))%*%t(SI+w*(SA-SI%*%R))
  
  ############# Means interpolation ##############  
  
  M<-(1-w)*M1+w*M2
  
  
  ############ Residual standardization #########
  
  Res1x<-solve(sqrtM(C1))%*%Res1%*%(solve(sqrtM(T1)))
  Res2x<-solve(sqrtM(C2))%*%Res2%*%(solve(sqrtM(T2)))
  
  Resx<-(1-w)*Res1x+w*Res2x
  ################################################  
  
  # Interpolation 
  
  Y<-(sqrtM(C3))%*%Resx%*%(sqrtM(T3))+M
  
  Spect_path<-c(Spect_path,list(Y))
  
}

# save the interpolated path in a matlab file for sound reconstruction.

#writeMat("Spect_path.mat",Spect_path1=Spect_path[[1]],Spect_path2=Spect_path[[2]],Spect_path3=Spect_path[[3]],Spect_path4=Spect_path[[4]],Spect_path5=Spect_path[[5]],Spect_path6=Spect_path[[6]],Spect_path7=Spect_path[[7]],Spect_path8=Spect_path[[8]],Spect_path9=Spect_path[[9]],Spect_path10=Spect_path[[10]],Spect_path11=Spect_path[[11]])


#### Interpolation between a recorded speaker and the transformation in a different language ########


# Speaker
spk<-1

# Residual for the speaker

Res1<-Data2[L1[spk],,]-M1



Spect_proj<-list()

for ( h in 1:length(weights)){
  
  w<-weights[h]
  
  ############## Frequency covariance ###################
  
  
  SI<-sqrtM(C1)
  SA<-sqrtM(C2)
  C<-t(SA)%*%SI
  SVD<-svd(C)
  R<-SVD$u%*%t(SVD$v)
  C3<-(SI+w*(SA-SI%*%R))%*%t(SI+w*(SA-SI%*%R))
  
  ############## Time covariance ###################
  
  SI<-sqrtM(T1)
  SA<-sqrtM(T2)
  C<-t(SA)%*%SI
  SVD<-svd(C)
  R<-SVD$u%*%t(SVD$v)
  T3<-(SI+w*(SA-SI%*%R))%*%t(SI+w*(SA-SI%*%R))
  
  
  ############# Means interpolation ##############  
  
  
  M<-(1-w)*M1+w*M2
  
  
  ################################################  
  
  # Interpolation 
  
  Y<-(sqrtM(C3)%*%solve(sqrtM(C1)))%*%Res1%*%(solve(sqrtM(T1))%*%sqrtM(T3))+M
  
  Spect_proj<-c(Spect_proj,list(Y))
  
}

#writeMat("Spect_proj.mat",Spect_proj1=Spect_proj[[1]],Spect_proj2=Spect_proj[[2]],Spect_proj3=Spect_proj[[3]],Spect_proj4=Spect_proj[[4]],Spect_proj5=Spect_proj[[5]],Spect_proj6=Spect_proj[[6]],Spect_proj7=Spect_proj[[7]],Spect_proj8=Spect_proj[[8]],Spect_proj9=Spect_proj[[9]],Spect_proj10=Spect_proj[[10]],Spect_proj11=Spect_proj[[11]])


######## Sound reconstruction ############ 

# Spectrogram sequence to be transformed

#Spectro<-Spect_proj
Spectro<-Spect_path

### recover phase info ####

#Sounds<-readMat('../Code and Data/PrunedLanguageDataset.mat')
wn<-160
rate<-16000
S<-specgram(Sounds$FR.F.5.F.c, n=wn, Fs=16000, window=gausswin(wn))
phase<-Arg(S$S)
output<-c()

for (k in 1:length(Spectro)){
  PSD<-matrix(NA,dim(phase)[1],dim(phase)[2])
  for (j in 1:dim(Spectro[[k]])[1]){
    PSD[j,]<-approx(time,Spectro[[k]][j,],xout=seq(0,1,len=dim(phase)[2]))$y
  }
  out<-sound_from_spectrogram(PSD,phase,rate,bit=16,wn)
  output<-c(output,out@left)
}
output<-Wave(output,right = numeric(0), samp.rate = rate, bit = 16)
writeWave(output,filename ="F2P_spk_digit5_V2.wav")
#writeWave(output,filename ="F2P_proj_digit5_V2.wav")


######## Sound reconstruction: example with one spectrogram ############ 

# Spectrogram sequence to be transformed

#Spectro<-Spect_proj
Spectro<-Data2[1,,]

### recover phase info ####

#Sounds<-readMat('../Code and Data/PrunedLanguageDataset.mat')
wn<-160
rate<-16000
S<-specgram(Sounds$FR.F.5.F.c, n=wn, Fs=16000, window=gausswin(wn))
phase<-Arg(S$S)

PSD<-matrix(NA,dim(phase)[1],dim(phase)[2])
for (j in 1:dim(Spectro)[1]){
  PSD[j,]<-approx(time,Spectro[j,],xout=seq(0,1,len=dim(phase)[2]))$y
}
out<-sound_from_spectrogram(PSD,phase,rate,bit=16,wn)
  
writeWave(out,filename ="sound_example.wav")
#writeWave(output,filename ="F2P_proj_digit5_V2.wav")



# to listen to original sounds:

output<-Wave(Sounds[[1]]*20000,right = numeric(0), samp.rate = rate, bit = 16)
writeWave(output,filename ="sound_example.wav")
