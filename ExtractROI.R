### Load packages
library(oro.nifti)

### Read in Tau 7 files
dat.dir <- "/Users/garyzhou/Downloads/Dropbox/MRI_TAU7"

### Load AAL mask
aal <- readNIfTI("/Users/garyzhou/Downloads/Dropbox/atlas/AAL2.nii")
dim(aal)

### Extract the ROI 4101+4102 from the AAL atlas corresponding to hippocampus
ind <- which(aal == 4101, arr.ind = TRUE)
ind <- rbind(ind, which(aal == 4102, arr.ind = TRUE))

### Create a 3D box where the ROI will be saved
ROI <- array(0, c(max(ind[,1])-min(ind[,1])+2, 
                  max(ind[,2])-min(ind[,2])+2,
                  max(ind[,3])-min(ind[,3])+2))

### Indices in the ROI 3D box corresponding to those in native space
ind1 <- ind
ind1[,1] <- ind[,1] - min(ind[,1])+1
ind1[,2] <- ind[,2] - min(ind[,2])+1
ind1[,3] <- ind[,3] - min(ind[,3])+1

# Obtain the paths to the PET scan for each subject
subj.dat <- dir(dat.dir, pattern = "*", full.names = TRUE)

# Obtain the subject IDs from the file names
subj.names <- dir(dat.dir, pattern = "*", full.names = FALSE)

# A list that will contain the 3D boxes of the ROI extracted from all subjects
roi.list <- c()

### Extract the ROI for each subject
for(subj in 1:length(subj.dat)){
  # Obtain the file name from the subject folder
  file.name <- dir(subj.dat[subj], pattern = "*", full.names = TRUE)[1]
  
  # Load the dataset for a subject
  dat <- readNIfTI(file.name)
  
  # Extract the region
  ROI[ind1] <- dat[ind]
  roi.list[[subj]] <- ROI
  print(subj)
}

###############################################################################
### Save ROI list and subject names
setwd("/Users/garyzhou/Downloads/Dropbox/hippo_tau7")
save(roi.list, file = "ROIlist_hippocampus.rda")
save(subj.names, file = "SubjectNames_hippocampus.rda")
