#!/bin/sh

PL=path/to/plink/1.09/
LO=path/to/liftOver

cd path/to/genotyping/data

## update the coordinates in the genotyping data from hg18 to hg19. Note that the chain is assumed to be in the same location as the lifover executable
## (http://genome.sph.umich.edu/wiki/LiftOver#Lift_Merlin.2FPLINK_format)
$LO/liftOver <(awk '{pos0=$4-1; print "chr"$1,pos0,$4,$2}' NFBC_dbGaP_20091127.bim) $LO/hg18ToHg19.over.chain.gz liftOver_output.BED liftOver_unlifted.BED

$PL/plink --noweb --allow-no-sex --bfile NFBC_dbGaP_20091127 --update-map <(awk '{gsub("chr","",$1); print $4,$1}' liftOver_output.BED) --update-chr --make-bed --out NFBC_dbGaP_20091127_hg19
$PL/plink --noweb --allow-no-sex --bfile NFBC_dbGaP_20091127_hg19 --update-map <(awk '{print $4,$3}' liftOver_output.BED) --make-bed --out NFBC_dbGaP_20091127_hg19

## run genome scans across the 9 traits
$PL/plink --bfile NFBC_dbGaP_20091127_hg19 --chr 1-22 \
--pheno NFBC66_phenotypes_raw.txt --all-pheno \
--covar NFBC66_covar_SexOCPG.txt \
--no-const-covar --allow-no-covars \
--linear hide-covar --out NFBC66_GWAS_adj_hg19

## The command above will generate files NFBC66_GWAS_adj_hg19.XX.assoc.linear, where XX is the corresponding QT (TG, HDL, LDL, CRP, GLU, INS, BMI, SBP, DBP)

## compute principal components
$PL/plink --bfile NFBC_dbGaP_20091127_hg19 --indep 50 5 2 --out NFBC_dbGaP_20091127_hg19_indep

$PL/plink --bfile NFBC_dbGaP_20091127_hg19 --extract NFBC_dbGaP_20091127_hg19_indep.prune.in --geno 0.01 --make-bed --out NFBC_dbGaP_20091127_hg19_pruned

$PL/plink --bfile NFBC_dbGaP_20091127_hg19_pruned --pca 40 --out NFBC_dbGaP_20091127_hg19_pca


