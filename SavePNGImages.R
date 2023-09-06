# Load Tau 7 ROI list and subject names
load("/Users/garyzhou/Downloads/Dropbox/hippo_tau7/ROIlist_hippocampus.rda")
load("/Users/garyzhou/Downloads/Dropbox/hippo_tau7/SubjectNames_hippocampus.rda")

### Masking threshold value derivation method (IQR-based)
# Define number of thresholds
num_thres <- 10

# Extract vectorized list of non-zero intensity values for all slices for each subject
nonzero_vals <- sapply(1:159, function(x) roi.list[[x]][roi.list[[x]] != 0])

# Summary of number of non-zero voxels per patient
summary(sort(sapply(1:159, function(x) length(nonzero_vals[[x]]))))

# Summary statistics for the non-zero intensity values across all slices for each subject
sapply(1:159, function(x) summary(nonzero_vals[[x]]))

# Get IQR for (non-zero) intensity vals for all patients 
nonzero_1st_quart <- summary(sort(sapply(1:159, function(x) summary(nonzero_vals[[x]]))))[2]
nonzero_3rd_quart <- summary(sort(sapply(1:159, function(x) summary(nonzero_vals[[x]]))))[5]
nonzero_iqr <- nonzero_3rd_quart - nonzero_1st_quart

# Evenly divide IQR by number of thresholds to get threshold values
threshold_val <- unname(sapply(1:(num_thres-1), function(x) nonzero_iqr*x/num_thres + nonzero_1st_quart))

### Visualization of unmasked images
for(subj in 1:length(subj.names)){
  subj.image <- roi.list[[subj]]
  setwd("/Users/garyzhou/Downloads/Dropbox/hippo_tau7/unmasked/images")
  dir.create(subj.names[subj])
  setwd(subj.names[subj])
  for(i in 1:(dim(mask)[2])){
    img <- mask[,i,]
    fname <- paste(i, "image.png", sep = "")
    png(fname)
    par(mar=c(0,0,0,0))
    image(img, axes = FALSE, col=gray.colors(2))
    dev.off()
  }
  print(subj)
}

### Apply mask and process
process_images <- function(q_val) {
  q <- round(q_val, 2)
  dir_base <- paste0("/Users/garyzhou/Downloads/Dropbox/hippo_tau7/Q", q, "/images")
  
  #Apply mask
  for(subj in 1:length(subj.names)){
    subj.image <- roi.list[[subj]]
    mask <- subj.image
    mask[subj.image > q] <- 1
    mask[subj.image <= q] <- 0
    
    #Create directory if path doesn't exist
    dir_subj <- file.path(dir_base, subj.names[subj])
    if(!dir.exists(dir_subj)) {
      dir.create(dir_subj)
    }
    
    #Create png images for each subject
    for(i in 1:(dim(mask)[2])){
      img <- mask[,i,]
      img1 <- matrix(0, dim(img)[1] + 15, dim(img)[2] + 10)
      img1[10:(9+dim(img)[1]), 5:(4+dim(img)[2])] <- img
      fname <- file.path(dir_subj, paste(i, "image.png", sep = ""))
      png(fname)
      par(mar=c(0,0,0,0))
      image(img1, axes = FALSE, col=gray.colors(2))
      dev.off()
    }
    print(subj)
  }
}

### Apply the function to each value in threshold_val (each maskiing threshold)
lapply(threshold_val, process_images)