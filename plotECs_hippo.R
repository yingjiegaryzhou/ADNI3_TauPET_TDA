### Plot ECs by diagnosis
### Run clinical.R before this to obtain the clinical markers

#Tau7
#To loop over all thresholds and iteratively got SEC matrices for each:
load_and_process <- function(threshold, clinical_data) {
  # File path for each thresh
  file_path = paste0("/Users/garyzhou/Downloads/Dropbox/hippo_tau7/Q", threshold, "/MRIECs_hippo.RData")
  
  # Load file
  load(file_path)
  
  # Obtain dimensions
  nrot = ncol(MRI_list[[1]]$EC)
  stepsize = nrow(MRI_list[[1]]$EC)
  
  # Generate matrix with k observations and nrot*stepsize column
  ECs = matrix(nrow = length(MRI_list), ncol = nrot*stepsize)
  rownames(ECs) = 1:nrow(ECs)
  
  # Place the curves in an nxp matrix with patient tags as the row names
  for(i in 1:nrow(ECs)){
    ECs[i,] = c(MRI_list[[i]]$EC)
    rownames(ECs)[i] = MRI_list[[i]]$name
  }
  
  # Remove the vectors of zeros where a new MRI slice begins
  ECs = ECs[,-seq(101,ncol(ECs),by=101)]
  
  ECs.dat = as.data.frame(ECs)
  ECs.dat$PTID = subj.names
  
  # Merge the clinical data and EC curves
  plot.dat = merge(clinical_data, ECs.dat, by = "PTID")
  
  return(plot.dat)
}

# Iterate over each threshold and store the results
results_list = list()
for (i in 1:length(threshold_val)) {
  results_list[[i]] = load_and_process(threshold_val[i], clinical)
}

#results_list contains all 'plot.dat' dfs in a list

##########################################################################################
##Demonstration of plotting EC curves for 1 threshold

#Q1
load("/Users/garyzhou/Downloads/Dropbox/hippo_tau7/Q1.04/MRIECs_hippo.RData")

#MRI_list output has two components for each of the k (159) patients: name and EC. 
#EC is a matrix with 101 observ (stepsize with 1 extra to denote start of new slice) and 72 cols (nrot)
nrot = ncol(MRI_list[[1]]$EC); stepsize = nrow(MRI_list[[1]]$EC)
ECs = matrix(nrow = length(MRI_list), ncol = nrot*stepsize) #generate matrix with k observ and nrot*stepsize column
rownames(ECs) = 1:nrow(ECs)
dim(ECs) #159x7272

### Place the Curves in an nxp Matrix with Patient Tags as the Row Names ###
for(i in 1:nrow(ECs)){
  ECs[i,] = c(MRI_list[[i]]$EC) #reduce matrix to length 7272 vector going down 
  #(first 101 elements correspond to first column (first rot angle))
  rownames(ECs)[i] = MRI_list[[i]]$name
}

### Remove the Vectors of Zeros where a New MRI Slice Begins ###
ECs = ECs[,-seq(101,ncol(ECs),by=101)]
ECs.dat = as.data.frame(ECs)
ECs.dat$PTID = subj.names
dim(ECs.dat)
dim(clinical) 

### Merge the clinical data and EC curves
plot.dat = merge(clinical, ECs.dat, by = "PTID") 
#Reduced to 102 patients with first 5 col being var from clinical dataset; 102x7205
dim(plot.dat)

### Add averages (median) for all directions
sub.cn_tot = plot.dat[plot.dat$DIAGNOSIS == 1,6:7205]
sub.cn_tot = sub.cn_tot[complete.cases(sub.cn_tot),] #39 complete cases
ec.cn_tot = apply(sub.cn_tot, 2, median)
sub.mci_tot = plot.dat[plot.dat$DIAGNOSIS == 2,6:7205]
sub.mci_tot = sub.mci_tot[complete.cases(sub.mci_tot),]
ec.mci_tot = apply(sub.mci_tot, 2, median) #35
sub.dem_tot = plot.dat[plot.dat$DIAGNOSIS == 3,6:7205]
sub.dem_tot = sub.dem_tot[complete.cases(sub.dem_tot),]
ec.dem_tot = apply(sub.dem_tot, 2, median) #28

#####Mean#####
### Add averages (mean) for all directions
sub.cn_tot = plot.dat[plot.dat$DIAGNOSIS == 1,6:7205]
sub.cn_tot = sub.cn_tot[complete.cases(sub.cn_tot),] #39 complete cases
ec.cn_tot = apply(sub.cn_tot, 2, mean)
sub.mci_tot = plot.dat[plot.dat$DIAGNOSIS == 2,6:7205]
sub.mci_tot = sub.mci_tot[complete.cases(sub.mci_tot),]
ec.mci_tot = apply(sub.mci_tot, 2, mean) #35
sub.dem_tot = plot.dat[plot.dat$DIAGNOSIS == 3,6:7205]
sub.dem_tot = sub.dem_tot[complete.cases(sub.dem_tot),]
ec.dem_tot = apply(sub.dem_tot, 2, mean) #28

##############################################################################
#Direction-Step Combo Example
cols = c("#5ab4ac", "#d8b365", "#b2182b")
#Plot for 2 groups
par(mfrow=c(3,3))
plot(1:7200, ec.dem_tot, col = cols[3], type = "l",
     lwd = 1, bty = "n",ylab = "SECT Value",
     xlab = "Direction-Step Combination", main="a=1.04", ylim = c(0,4)) 
lines(1:7200, ec.dem_tot, col = cols[3], lwd = 2)
lines(1:7200, ec.cn_tot, col = cols[1], lwd = 2)

###############################################

#Example of one curve (for poster)
plot(1:100, plot.dat[5,6:105], col = "red", type = "l",
     lwd = 1, bty = "n",ylab = "SEC Value",
     xlab = "Stepsize", ylim = c(0,4)) 

cols = c("#5ab4ac", "#d8b365", "#b2182b")
#Plot of EC curves for all 3 groups
plot(1:7200, ec.dem_tot, col = cols[3], type = "l",
     lwd = 1, bty = "n",ylab = "SECT Values",
     xlab = "Direction-Step Combination", ylim = c(0,4)) 
lines(1:7200, ec.dem_tot, col = cols[3], lwd = 1)
lines(1:7200, ec.mci_tot, col = cols[2], lwd = 1)
lines(1:7200, ec.cn_tot, col = cols[1], lwd = 1)
legend("bottom",legend = c("CN","MCI","Dementia"), col = cols, lwd = c(1, 1, 1), horiz = TRUE,bty = "n")

cols = c("blue", "red", "green")
#Plot for 2 groups
plot(1:7200, ec.dem_tot, col = cols[3], type = "l",
     lwd = 1, bty = "n",ylab = "Smooth Euler Characteristic Transform (SECT)",
     xlab = "Sublevel Sets", ylim = c(0,4)) 
lines(1:7200, ec.dem_tot, col = cols[3], lwd = 2)
lines(1:7200, ec.cn_tot, col = cols[1], lwd = 2)

#Plot all curves for 2 groups
matplot(t(sub.dem_tot), type = "l", lty = 1, col = "green", xlab = "Column Index", ylab = "SECT Value")
matlines(t(sub.cn_tot), type = "l", lty = 1, col = "blue")
matlines(t(sub.mci_tot), type = "l", lty = 1, col = "red")
legend("topright", legend = c("AD", "CN"), lty = 1, col = c("green", "blue"), title = "Diagnosis Group", cex=0.4)

plot(1:7200, sub.dem_tot, col = cols[3], type = "l",
     lwd = 1, bty = "n",ylab = "Smooth Euler Characteristic Transform (SECT)",
     xlab = "Sublevel Sets", ylim = c(0,4)) 
lines(1:7200, ec.dem_tot, col = cols[3], lwd = 2)
lines(1:7200, ec.cn_tot, col = cols[1], lwd = 2)

#Q1.04
setwd("/Users/garyzhou/Downloads/Dropbox/hippo_tau7/Q1.04")
save(ec.cn_tot, file = "CN_SECT_Values.rda")
save(ec.mci_tot, file = "MCI_SECT_Values.rda")
save(ec.dem_tot, file = "DEM_SECT_Values.rda")
plot.dat.1.04_hippo <- plot.dat
save(plot.dat.1.04_hippo, file = "plotdat1.04.rda")


##########################################################################################
##### Plot all EC curves colored by DX for angle 1 #####
cols = c("blue", "red", "green")
tmp = matrix(plot.dat[1,6:105], 1, 100)
plot(1:100, tmp, col = cols[plot.dat$DIAGNOSIS[1]], type = "l",
     lwd = 2, bty = "n",ylab = "Smooth Euler Characteristic Transform (SECT)",
     xlab = "Sublevel Sets", ylim = c(1.5,4)) 

for(i in 2:dim(plot.dat)[1]){
  tmp = matrix(plot.dat[i,6:105], 1, 100)
  lines(1:100, tmp,col = cols[plot.dat$DIAGNOSIS[i]], lty = 1, lwd = 2)
}

##### Plot median EC curves colored by DX for angle 1 #####
sub.cn = plot.dat[plot.dat$DIAGNOSIS == 1,6:105]
sub.cn = sub.cn[complete.cases(sub.cn),] #39 complete cases
ec.cn = apply(sub.cn, 2, median)
sub.mci = plot.dat[plot.dat$DIAGNOSIS == 2,6:105]
sub.mci = sub.mci[complete.cases(sub.mci),]
ec.mci = apply(sub.mci, 2, median) #35
sub.dem = plot.dat[plot.dat$DIAGNOSIS == 3,6:105]
sub.dem = sub.dem[complete.cases(sub.dem),]
ec.dem = apply(sub.dem, 2, median) #28

plot(1:100, ec.dem, col = cols[3], type = "l",
     lwd = 3, bty = "n",ylab = "Smooth Euler Characteristic Transform (SECT)",
     xlab = "Sublevel Sets", ylim = c(1.5,4)) 
lines(1:100, ec.dem, col = cols[3], lwd = 3)
lines(1:100, ec.mci, col = cols[2], lwd = 3)
lines(1:100, ec.cn, col = cols[1], lwd = 3)
legend("top",legend = c("CN","MCI","Dementia"),col = cols,lty = c(1,2,4),lwd = 2,horiz = TRUE,bty = "n")

##################################################################################################################################
### Add variance for all directions
#Q1.04
sub.cn_tot_var = plot.dat.1.04[plot.dat.1.04$DIAGNOSIS == 1,6:7205]
sub.cn_tot_var = sub.cn_tot_var[complete.cases(sub.cn_tot_var),] #39 complete cases
ec.cn_tot_var = apply(sub.cn_tot_var, 2, var)
sub.mci_tot_var = plot.dat.1.04[plot.dat.1.04$DIAGNOSIS == 2,6:7205]
sub.mci_tot_var = sub.mci_tot_var[complete.cases(sub.mci_tot_var),]
ec.mci_tot_var = apply(sub.mci_tot_var, 2, var) #35
sub.dem_tot_var = plot.dat.1.04[plot.dat.1.04$DIAGNOSIS == 3,6:7205]
sub.dem_tot_var = sub.dem_tot_var[complete.cases(sub.dem_tot_var),]
ec.dem_tot_var = apply(sub.dem_tot_var, 2, var) #28

plot(1:7200, ec.dem_tot_var, col = cols[3], type = "l",
     lwd = 3, bty = "n",ylab = "Smooth Euler Characteristic Transform (SECT)",
     xlab = "Sublevel Sets", ylim = c(0,0.5)) 
lines(1:7200, ec.dem_tot_var, col = cols[3], lwd = 3)
lines(1:7200, ec.mci_tot_var, col = cols[2], lwd = 3)
lines(1:7200, ec.cn_tot_var, col = cols[1], lwd = 3)

legend("top",legend = c("CN","MCI","Dementia"),col = cols,lty = c(1,2,4),lwd = 2,horiz = TRUE,bty = "n")

#Q1.00
setwd("/Users/garyzhou/Downloads/Dropbox/hippocampus_tau7/Q1.04")
save(ec.cn_tot_var, file = "CN_SECT_Values_Var.rda")
save(ec.mci_tot_var, file = "MCI_SECT_Values_Var.rda")
save(ec.dem_tot_var, file = "DEM_SECT_Values_Var.rda")




