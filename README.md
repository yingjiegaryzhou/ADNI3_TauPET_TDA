# ADNI3 Texture Analysis of Tau PET Scans through SECT
Smooth Euler Characteristic Transform and subsequent Structured Functional Principal Component Analysis of baseline ADNI3 Tau PET scans. The files in this repo are samples for one ROI (hippocampus).

Order of Files
1. TAU_MRI.R: preprocessing
2. ExtractROI.R: further preprocessing
3. SavePNGImages.R: save as image files
4. EC3D.R: SECT algorithm
5. SECT_Hippo.R: SECT implementation
6. clinical.R: clean clinical and diagnosis files
7. plotECs_hippo.R: generate SECT matrices and plots
8. SFPCA_Code_Example.R: SFPCA implementation