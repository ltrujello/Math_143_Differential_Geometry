library(R.matlab)
library(fields)
library(rgl)

#############################################################################
########   This script implements the analysis of acoustic phonetic data ####
########           described in Pigoli et al., 2015                      ####
########                                                                 ####
#############################################################################

source('acoustic_functions.R')

# Import data from Matlab save files

#Data<-readMat('WarpedPSD.mat')

# Import data from R workspace

load('WarpedPSD.RData')

Data2<-array(0,c(219,81,100))

for (i in 1:219){
  Data2[i,,]<-Data[[1]][i][[1]]
}

# Import speaker information

Speaker_info<-read.table('SpeakerInfo.txt',header=T)
attach(Speaker_info)



# Estimation of word means

MFr<-array(NA,dim(Data2))
MI<-array(NA,dim(Data2))
MP<-array(NA,dim(Data2))
MSA<-array(NA,dim(Data2))
MSI<-array(NA,dim(Data2))


for ( i in 1:10){ 
  
  Fr1<-which(language=='French'& digit==i)
  I1<-which(language=='Italian'& digit==i)
  P1<-which(language=='Portuguese'& digit==i)
  SA1<-which(language=='American_Spanish'& digit==i)
  SI1<-which(language=='Iberian_Spanish'& digit==i)
  
  
  
  MFr[i,,]<-apply(Data2[Fr1,,],c(2,3),mean)
  
  MI[i,,]<-apply(Data2[I1,,],c(2,3),mean)
  
  
  MP[i,,]<-apply(Data2[P1,,],c(2,3),mean)
  
  MSA[i,,]<-apply(Data2[SA1,,],c(2,3),mean)

  MSI[i,,]<-apply(Data2[SI1,,],c(2,3),mean)
  
}


### Estimation of the separable language covariances

D<-array(0,dim(Data2))

for ( i in 1:10){ 
  
  Fr1<-which(language=='French'& digit==i)
  I1<-which(language=='Italian'& digit==i)
  P1<-which(language=='Portuguese'& digit==i)
  SA1<-which(language=='American_Spanish'& digit==i)
  SI1<-which(language=='Iberian_Spanish'& digit==i)
  
  
  
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


Fr<-which(language=='French')
I<-which(language=='Italian')
P<-which(language=='Portuguese')
SA<-which(language=='American_Spanish')
SI<-which(language=='Iberian_Spanish')


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
  pv_freq[k]<-Cov_test(D[lang,,],digit[lang],dir=1,M=1000)
}

pv_freq

# P-vlues of the permutation tests for the equality of the time covariance across the words
# for each language


pv_time<-rep(NA,length(lg))

for (k in 1:length(lg)){ 
  lang<-which(language==lg[k])  
  pv_time[k]<-Cov_test(D[lang,,],digit[lang],dir=2,M=1000)
}

pv_time





#################################################################################################



##### Transformation of a speaker log-spectrograms from French to Portuguese #############

l1<-"French"
l2<-"Portuguese"

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

# Residual log-spectrogram 

Res<-Data2[L1[1],,]-M1

# Transformation:
  
PS<-(sqrtM(C2)%*%solve(sqrtM(C1)))%*%Res%*%(solve(sqrtM(T1))%*%sqrtM(T2))+M2

# save it in a matlab file for sound reconstruction
writeMat("FR2I_proj.mat",F2P=PS)





############ Interpolation and extrapolation of sound surfaces ################

l1<-"French"
l2<-"Portuguese"

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
spk2<-2

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

writeMat("Spect_path.mat",Spect_path1=Spect_path[[1]],Spect_path2=Spect_path[[2]],Spect_path3=Spect_path[[3]],Spect_path4=Spect_path[[4]],Spect_path5=Spect_path[[5]],Spect_path6=Spect_path[[6]],Spect_path7=Spect_path[[7]],Spect_path8=Spect_path[[8]],Spect_path9=Spect_path[[9]],Spect_path10=Spect_path[[10]],Spect_path11=Spect_path[[11]])


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

writeMat("Spect_proj.mat",Spect_proj1=Spect_proj[[1]],Spect_proj2=Spect_proj[[2]],Spect_proj3=Spect_proj[[3]],Spect_proj4=Spect_proj[[4]],Spect_proj5=Spect_proj[[5]],Spect_proj6=Spect_proj[[6]],Spect_proj7=Spect_proj[[7]],Spect_proj8=Spect_proj[[8]],Spect_proj9=Spect_proj[[9]],Spect_proj10=Spect_proj[[10]],Spect_proj11=Spect_proj[[11]])






