#########SFPCA########################################################################
#Equally spaced grid from 0 to 100 (stepsize is the analog of time here)
T <- seq(1, 100, 1) #stepsize
p <- length(T)
L1=67 #67 subjects; I
J=L2=9 #9 thresholds within each subject; J
L3=72 #72 directions within each threshold; K

I_index=rep(1:L1, each=L2*L3) #Index for each subject; numeric vector of 67 subjects, each repeated 72*9 times
n=length(I_index) #samp size; n=43416 curves
IJ_index=rep(1:(L1*L2), each=L3) #IJ_index = subject-threshold index: 1 to (67*9) with each number duplicated 72 times

I=max(I_index); I #I=67; # of subj
n_I0=table(I_index); n_I0 #648 threshold-direction pairs in each of 67 subj; #samp size (N=(72*9)) for each I 
IJ=max(IJ_index); IJ #603 subject-thresh pairs
n_IJ=table(IJ_index); n_IJ #72 direction in each subj-thresh pair; #samp size (N=72) for each subj-thresh index (67*9)

#Standardize final matrix
final_matrix_stdized <- scale(final_matrix, center=apply(final_matrix, 2, mean), scale=F) #demean (columnn-wise: subtract mean from each element), but do not scale
head(final_matrix_stdized)[,1:5] 

smooth.y <- TRUE
set.seed(66)
### smooth the raw data w/ spline
if(smooth.y){
  Y_new <-list()
  Y_new<-lapply(1:n, function(i) {m<-smooth.spline(as.numeric(final_matrix_stdized[,i]));  return(predict(m)$y)})
  Y <-Reduce(cbind,Y_new)
}

dim(Y)
#Y is 100x43416 matrix
#Each column represents an EC curve as a function of stepsize for a specific direction, threshold, and subject

#Save Y as Rdata file
Y_matrix_SFPCA <- Y
dim(Y_matrix_SFPCA)
setwd("/Users/garyzhou/Downloads/Dropbox/hippo_tau7")
save(Y_matrix_SFPCA, file = "Y_matrix_SFPCA.rda")

#To avoid waiting for Y, load "Y_matrix_SFPCA.rda" directly and proceed:
Y <- Y_matrix_SFPCA


#Covariance estimation
k1=sum(n_IJ^2) 
Y_IJ <- t(rowsum(t(Y),IJ_index))
k2=sum(n_I0^2)
Y_I <- t(rowsum(t(Y),I_index))

##############
#Might have to allocate more memory to the calculation
#memory.limit(size=8000)
library(usethis) 
usethis::edit_r_environ()
##############

#Original Functions
#Corresponds directly to formulas given in manuscript
HW <- (Y%*%diag(n_IJ[IJ_index])%*%t(Y))-(Y_IJ%*%t(Y_IJ))*2/(k1-n) 
HU <- (Y%*%diag(n_I0[I_index])%*%t(Y)-(Y_I)%*%t(Y_I)-(k1-n)/2*HW)*2/(k2-k1)
HX <- (n*Y%*%t(Y)-rowSums(Y)%*%t(rowSums(Y))-(k1-n)/2*HW-(k2-k1)/2*HU)*2/(n^2-k2)

#Memory-efficient crossprod of transpose
HW <- (Y%*%tcrossprod(diag(n_IJ[IJ_index]), Y))-(tcrossprod(Y_IJ, Y_IJ))*2/(k1-n)
HU <- (Y%*%tcrossprod(diag(n_I0[I_index]),t(Y))-(tcrossprod(Y_I, Y_I))-(k1-n)/2*HW)*2/(k2-k1)
HX <- (n*tcrossprod(Y,Y)-tcrossprod(rowSums(Y), rowSums(Y))-(k1-n)/2*HW-(k2-k1)/2*HU)*2/(n^2-k2) #no diagonalization so takes least amt of time

KW <- HW/2 #Covar operator for direction-level
KU <- (HU - HW)/2 #Covar operator for threshold-level
KX <- (HX - HU)/2 #Covar operator for subj-level
#Outputs are all 100x100 covar matrices

#To avoid waiting, load "KW.standardized.twogroup.rda" and the other two RDA files directly and proceed:

#Save
KW.standardized <- KW
KU.standardized <- KU
KX.standardized <- KX
setwd("/Users/garyzhou/Downloads/Dropbox/hippo_tau7")
save(KW.standardized, file = "KW.standardized.twogroup.rda")
save(KU.standardized, file = "KU.standardized.twogroup.rda")
save(KX.standardized, file = "KX.standardized.twogroup.rda")

#Number of principal components to include; choose # to account for 95% of variability. Here we use two just as example (covers less than 90%)
N1 <- N2 <- N3 <- 2

#The eigenfunctions and the eigenvalues then can be estimated.
phi3e <- eigen(KW.standardized)$vectors[,1:N3] #direction; 100 eigenvectors for each column
phi2e<- eigen(KU.standardized)$vectors[,1:N2] #threshold
phi1e <- eigen(KX.standardized)$vectors[,1:N1] #subject

#The first principal component for layer level over time (stepsize) can be plotted. 
T_new <- seq(1, 100, 1)

par(mfrow=c(2,3))
plot(T_new, -phi3e[,1], type="l", xlab="Stepsize", ylab="Eigenvector", main="PC1 for Direction")
plot(T_new, -phi2e[,1], type="l", xlab="Stepsize", ylab="", main="PC1 for Threshold")
plot(T_new, -phi1e[,1], type="l", xlab="Stepsize", ylab="", main="PC1 for Subject")

#2nd FPC
plot(T_new, -phi3e[,2], type="l", xlab="Stepsize", ylab="Eigenvector", main="PC2 for Direction")
plot(T_new, -phi2e[,2], type="l", xlab="Stepsize", ylab="", main="PC2 for Threshold")
plot(T_new, -phi1e[,2], type="l", xlab="Stepsize", ylab="", main="PC2 for Subject")

#Eigenvalues
### Direction ###
prop_eigenval_KW <- eigen(KW.standardized)$values/sum(eigen(KW.standardized)$values)
prop_eigenval_KW #% of level-specific variability explained by FPC: 53.21% 1st FPC; 12.92% 2nd FPC; etc.
cumsum(prop_eigenval_KW) #91% of level-specific variability explained by first 8 FPC

#Cum variance explained by FPC plot
plot(
  x = seq(1:length(prop_eigenval_KW)), y = prop_eigenval_KW,
  type = "o",
  xlab = "FPC for Direction", ylab = "Variance")

#Percent of variability explained by Direction-Level - 44.51%
sum(eigen(KW.standardized)$values)/(abs(sum(eigen(KW.standardized)$values)) + abs(sum(eigen(KU.standardized)$values)) + abs(sum(eigen(KX.standardized)$values)))

### Threshold ###
prop_eigenval_KU <- eigen(KU.standardized)$values/sum(eigen(KU.standardized)$values)
prop_eigenval_KU
plot(
  x = seq(1:length(prop_eigenval_KU)), y = prop_eigenval_KU,
  type = "o",
  xlab = "FPC for Threshold", ylab = "Variance")

#Percent of variability explained by Threshold-Level - 50%
abs(sum(eigen(KU.standardized)$values))/(abs(sum(eigen(KW.standardized)$values)) + abs(sum(eigen(KU.standardized)$values)) + abs(sum(eigen(KX.standardized)$values)))

### Subject ###
prop_eigenval_KX <- eigen(KX.standardized)$values/sum(eigen(KX.standardized)$values)
prop_eigenval_KX
plot(
  x = seq(1:length(prop_eigenval_KX)), y = prop_eigenval_KX,
  type = "o",
  xlab = "FPC for Subject", ylab = "Variance")

#Percent of variability explained by Subject-Level - 5.49%
abs(sum(eigen(KX.standardized)$values))/(abs(sum(eigen(KW.standardized)$values)) + abs(sum(eigen(KU.standardized)$values)) + abs(sum(eigen(KX.standardized)$values)))

# Project the data onto the principal component space and get FPC scores
### Subject ###
fpca_scores <- t(phi1e) %*% Y_matrix_SFPCA
# Print the FPCA scores
dim(fpca_scores)
plot(fpca_scores[1,], fpca_scores[2,], xlab = "FPCA 1 Scores", ylab = "FPCA 2 Scores", main = "FPCA Scores Plot Subject")
# Color code the points based on column index; first 25272 columns correspond to CN group; the rest to AD group
plot(fpca_scores[1, 1:25272], fpca_scores[2, 1:25272], col = "red", xlab = "FPCA 1 Scores", ylab = "FPCA 2 Scores", main = "FPCA Scores Plot Subj")
points(fpca_scores[1, 25273:43416], fpca_scores[2, 25273:43416], col = "blue")

#Threshold
fpca_scores_thres <- t(phi2e) %*% Y_matrix_SFPCA
plot(fpca_scores_thres[1, 1:25272], fpca_scores_thres[2, 1:25272], col = "red", xlab = "FPCA 1 Scores", ylab = "FPCA 2 Scores", main = "FPCA Scores Plot Threshold")
points(fpca_scores_thres[1, 25273:43416], fpca_scores_thres[2, 25273:43416], col = "blue")

#Direction
fpca_scores_dir <- t(phi3e) %*% Y_matrix_SFPCA
plot(fpca_scores_dir[1, 1:25272], fpca_scores_dir[2, 1:25272], col = "red", xlab = "FPCA 1 Scores", ylab = "FPCA 2 Scores", main = "FPCA Scores Plot Direction")
points(fpca_scores_dir[1, 25273:43416], fpca_scores_dir[2, 25273:43416], col = "blue")


