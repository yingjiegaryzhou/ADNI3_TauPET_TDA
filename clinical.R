### Load the ADNI dataset that contains clinical and demographic information 
setwd("/Users/garyzhou/Downloads/Dropbox")
adni = read.csv("ADNIMERGE.csv")

### Obtain the diagnosis from the DXSUM database
### Description from DATADIC: Which best describes the participant's current diagnosis?
### 1=CN; 2=MCI; 3=Dementia
### The variable DXDDUE: 3b.  Suspected cause of dementia
### 1=Dementia due to Alzheimer's Disease; 2=Dementia due to other etiology
### The variable DXDSEV: 3a.  Dementia Severity - Clinician's Impression
### 1=Mild; 2=Moderate; 3=Severe

### Extract data for ADNI3
adni = adni[adni$ORIGPROT == "ADNI3",]
adni.dx = read.csv("DXSUM_PDXCONV_ADNIALL.csv")

### Load subject names
load("/Users/garyzhou/Downloads/Dropbox/hippocampus_tau7/SubjectNames_hippocampus.rda")

#### Dataset to save the clinical information
clinical = c()
dx = c()
for(subj in 1:length(subj.names)){
  
  ind = which(adni$PTID == subj.names[subj])[1]
  ind1 = which(adni.dx$PTID == subj.names[subj])
  clinical = rbind(clinical, adni[ind, c("PTID", "VISCODE", "AGE", "MMSE")])
  dx = rbind(dx, adni.dx[ind1, c("PTID", "VISCODE2", "DIAGNOSIS")])
  
}
names(dx)[2] = "VISCODE"
clinical = merge(clinical, dx, by = c("PTID", "VISCODE"))
#493 participants, 491 baseline and 12-month scans

