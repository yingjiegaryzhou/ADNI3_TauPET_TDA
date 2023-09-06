# Package names
packages <- c("fslr", "oro.nifti", "neurobase", "RNifti", "abind", "stringr", "devtools")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load standard packages
library(fslr)
library(oro.nifti)
library(neurobase)
library(RNifti)
library(abind)
library(stringr)
library(devtools)

# ANTsR installation guide
# https://johnmuschelli.com/neuroc/installing_ANTsR/index.html
library(ANTsR)

# Extrantsr installation
# devtools::install_github("muschellij2/extrantsr")
library(extrantsr)

options(fsl.path ='/usr/local/fsl')
options(fsl.outputtype = 'NIFTI_GZ')

#####################################################################################
### Steps to register TAU PET image
# MRI used for masking in this method. Any value > 0 is a masked area of brain.

#1. Have template (t0), MRI (m0) and TAU PET image (tau0) ready
#2. Register tau0 to m0, obtain tau1
#2. Skull strip t0, obtain image t1
#3. Skull strip m0, obtain image m1
#4. Register m1 (ss-mri) to t1 (ss-template), obtain image m2, and mask (binary components of m2>0)
#5. Apply the mask in step 4 to tau1, obtain tau2
#6. Apply the transformation in step 4 (m1 registered to t1) on tau1, obtain the pre-processed TAU PET

# Set the path to project location (Windows)
setwd('/mnt/d')

# Get the template
aal.atlas <- readNIfTI('MNI152_T1_2mm.nii')

# Skull strip and remove neck
aal.template <- extrantsr::fslbet_robust(aal.atlas, remover = 'double_remove_neck')

### Get TAU PET file names of images
cur_path <- getwd()
image_path <- paste0(cur_path, '/all_tau/ADNI')
file_names <- dir(image_path)
image_names <- list()
for(i in file_names){
  cpath <- paste0(image_path, '/', i)
  cpath <- paste0(cpath, '/', dir(cpath))
  sub_files <- dir(cpath)
  sub_names <- c()
  for(j in 1:length(sub_files)){
    spath <- paste0(cpath, '/', sub_files[j])
    spath <- paste0(spath, '/', dir(spath))
    spath <- paste0(spath, '/', dir(spath))
    sub_names <- c(sub_names, spath)
  }
  image_names[[i]] <- sub_names
}

### Get MRI file names
mri_path <- paste0(cur_path, '/MRI_TAU_USE/ADNI')
file_names <- dir(mri_path)
mri_names <- list()
for(i in file_names){
  cpath <- paste0(mri_path, '/', i)
  cpath <- paste0(cpath, '/', dir(cpath))
  #get the first MRI of that subject
  cpath <- paste0(cpath, '/', dir(cpath)[1])
  nn <- length(dir(cpath))
  cpath <- paste0(cpath, '/', dir(cpath)[nn])
  cpath <- paste0(cpath, '/', dir(cpath))
  mri_names[[i]] <- cpath
}

### Check for differences between the datasets (names should match up)
all_mri <- names(mri_names)
all_tau <- names(image_names)
setdiff(all_mri, all_tau)
setdiff(all_tau, all_mri)

## Register the first one for each subject
## For each MRI, skull strip, apply mask to TAU PET,
## Register mask to template, apply transformation to TAU
for(name in names(mri_names)){
  new_folder <- paste0('MRI_TAU7/', name)
  if(file.exists(new_folder)){next}
  dir.create(new_folder)
  mri_file <- mri_names[[name]][1]
  tau_file <- image_names[[name]][1]
  mri.img <- readNIfTI(mri_file)
  tau.img <- readNIfTI(tau_file)
  out.file <- paste0(new_folder, '/aal', '.nii')
  
  # Brain extraction to MRI
  mri.mask <- extrantsr::fslbet_robust(mri.img, remover = 'double_remove_neck')
  
  # Register original tau to original mri
  tau.mri.orig <- flirt(infile = tau.img, reffile = mri.img)
  
  # Apply stripped MRI to ss-template w/ 7 dof
  mri.tmp <- flirt(infile = mri.mask, reffile = aal.template, omat = 'mri.mat', dof = 7)
  
  # Apply resulting mask to TAU(binary mask)
  mask <- mri.mask > 0
  tau.mask <- fslmask(tau.mri.orig, mask = mask)
  
  # Apply transformation from MRI registration to processed TAU
  aal.reg <- flirt_apply(infile = tau.mask, reffile = aal.template, initmat = 'mri.mat',
                         outfile = out.file)
}

### Visualize the pipeline each step of the way to doublecheck 
orthographic(mri.img)
orthographic(tau.img)
orthographic(aal.atlas)

orthographic(mri.mask)
orthographic(mask)

orthographic(tau.mri.orig)
orthographic(mri.tmp)

orthographic(tau.mask)
orthographic(aal.reg)
