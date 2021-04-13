# TwoPhase_postGWAS
R code for selected results from "Two-phase sample selection strategies for design and analysis in post-genome wide association fine-mapping studies"

This repository contains one file and two folders:
* twoPhaseGAS_1.07.tar.gz
* Sim2_PracticalScenario
* App_NFBC66

### twoPhaseGAS_1.07.tar.gz
This file contains the source code for a preliminary version of an R package implementing the proposed methods. This package is necessary for the codes above. Code that describes how to install the package is contained in file Sim2_PracticalScenario/Realistic_Simulation.R

### Sim2_PracticalScenario
Contains the necessary codes and information to reproduce the results in Section 4 of the manuscript. The main code is "Realistic_Simulation.R" while files "Step_[1-4].R" contain the auxiliary functions with the details of each step. These scripts rely on file "data_Realistic_R=1K_N=5K.RData", which containts the simulation replicates, to fully reproduce the results from the manuscript. This file is available as Supplementary Data along with the manuscript. However, the file can be generated using script "Step_0_Data_generation.R", although the results will not be reproducible. Files "genetic_map_GRCh37_chr16.txt.gz" and "EUR.chr16.HERPUD1-CETP.genotypes.vcf.gz" can be used to generate such new set of replicates. 

### App_NFBC66
This folder contains detailed bash scripts and R codes that describe the analysis steps for the illustration on the North Finland
Birth Cohort of 1966 (NFBC66) data. Note that no data from the study were made available. Researchers can request access to these data through dbGAP on [this link](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000276.v2.p1), which containts more information of the study.
