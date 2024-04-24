# DW-SPECIAL

Material associated to the following publication: https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.29805 

Please cite the above should you use any material from this repository.

Jessie Mosso, EPFL, 23/04/2024

## Third party softwares:
- FID-A, Jamie Near: https://github.com/CIC-methods/FID-A
- LCModel

## Files description: 
> A_Create_study_bruker_JM_07122021_VF.m

is used to read ser/fid files and methods from Bruker experiments folders 

>B_Save_rawdatajob0_bruker_JM_28102023.m

is used to read rawdatajob0 from Bruker experiment folders (is optional, is used only when Bruker's online corrections (e.g ECC) have to be undone) 

> C_Preprocessing_FIDA_DWSPECIAL_13032024.m

is used to process the data out of **A_Create_study_bruker_JM_07122021_VF.m** and **B_Save_rawdatajob0_bruker_JM_28102023.m** and used the functions: 

> convertNicoStudyFidAformat_isisoff.m
> convertNicoStudyFidAformat_isison.m

and the functions in the folder: **custom FID-A functions** 

> Explist_dMRS_DWSPECIAL_JessieMosso_18122022.xlsx

contains the list with experiment numbers for each rat used in this paper
## Folders description:
> /custom_FID-A_functions
> 
contains custom FID-A functions used in **C_Preprocessing_FIDA_DWSPECIAL_13032024.m**

> /LCM_basisset_controlfile
> 
contains the control file and basis set for LCModel quantification of the output of **C_Preprocessing_FIDA_DWSPECIAL_13032024.m**, that can be found in *processed/sum/*

> /raw
> 
contains the raw data read after steps **A_Create_study_bruker_JM_07122021_VF.m** and **B_Save_rawdatajob0_bruker_JM_28102023.m** 

> /processed
> 
contains the processed data after **C_Preprocessing_FIDA_DWSPECIAL_13032024.m**
