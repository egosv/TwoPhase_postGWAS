#!/bin/sh

export MDIR=/out/path/
export GCDIR=/path/to/gotcloud
export BAMLOC=/path/to/bam/files


#### Variant calling following ( https://genome.sph.umich.edu/wiki/GotCloud:_Variant_Calling_Pipeline )


####
#### Step 1. Make sure I have all the imput data:
####
## Aligned/Processed/Recalibrated BAM files
## BAM list file containing Sample IDs & BAM file names
## Reference files

### Aligned/Processed/Recalibrated BAM files

#*** This pipeline assumes that bam files have been already aligned (which is the case since I used sam-dump from sra-tools to do so). If not, fastq files are needed and an align call to gotcloud is also needed beforehand ***#

#*** Note, however, that the bam files provided have NOT been indexed, so, I'm using a sub-pipeline from GotCloud to do so. I'm going to do it per gene as well ***#
## Since the pipeline was not available, I'm just going to simply index using samtools

### Note that this used samtools v 0.1-19 for the sorting (which is presumably the one used in gotcloud)
cd $BAMLOC
 LPL APOA1 APOA5
for gene in GCKR; do
	awk -v gene=$gene '{print $0"_"gene"_aligned"}' path/to/SRA_unique.txt | xargs -t -n1 -P30 -I % sh -c 'samtools sort -f %.bam %_sorted.bam; samtools index %_sorted.bam'
done

### for a specific SRR
#srr=SRS422338
#srr=SRS427159
for gene in GCKR LPL APOA1 APOA5; do
	samtools sort -f ${srr}_${gene}_aligned.bam ${srr}_${gene}_aligned_sorted.bam; 
	samtools index ${srr}_${gene}_aligned_sorted.bam
done
# ls -lh ${srr}_*_aligned*.bam

### BAM list file
## Create it using file $MDIR/scripts/create_GotCloud_BAMlist.R

## Files are located in $MDIR/scripts/BAM_list_*.txt (one per gene)

### Use the latest reference files as per defaults
# cd ~/software/gotcloud
# wget ftp://anonymous@share.sph.umich.edu/gotcloud/ref/hs37d5-db142-v1.tgz
# tar xzf hs37d5-db142-v1.tgz

####
#### Step 2. Call gotcloud per gene:
####

#gene=GCKR
ulimit -n 6000
for gene in GCKR LPL APOA1 APOA5; do
	echo $gene
	mkdir -p $MDIR/processed_data/GotCloud/${gene}
	cd $MDIR/processed_data/GotCloud/${gene}
	## Define region per gene. Note that all these coordinates are in GRCh37.p13, i.e. hg19

	if [ "$gene" == "GCKR" ]; then chr=2; lo=$((27719470-5000)); up=$((27746556+5000)); 
	elif [ "$gene" == "LPL" ]; then chr=8; lo=$((19796582-5000)); up=$((19824770+5000)); 
	elif [ "$gene" == "APOA1" ]; then chr=11; lo=$((116706467-5000)); up=$((116708338+5000)); 
	elif [ "$gene" == "APOA5" ]; then chr=11; lo=$((116660086-5000)); up=$((116663136+5000)); 
	else echo "Not a gene we selected"; 
	fi
	
	## Create BED file per gene
	printf ${chr}'\t'${lo}'\t'${up}'\n' > $MDIR/scripts/${gene}.BED
	
	## Now, create a configuration file per gene
	CG_conf_file=$MDIR/scripts/GC_file_${gene}.conf

	echo "UNIFORM_TARGET_BED = $MDIR/scripts/${gene}.BED" > ${CG_conf_file}
	echo "WGS_SVM = TRUE" >> ${CG_conf_file}
	
	# Add these options for single sample processing as suggested in the documentation
	echo "UNIT_CHUNK = 20000000" >> ${CG_conf_file}
	echo "MODEL_GLFSINGLE = TRUE" >> ${CG_conf_file}
	echo "MODEL_SKIP_DISCOVER = FALSE" >> ${CG_conf_file}
	echo "MODEL_AF_PRIOR = TRUE" >> ${CG_conf_file}
	echo "VCF_EXTRACT = \$(REF_DIR)/snpOnly.vcf.gz" >> ${CG_conf_file}
	echo "EXT = \$(REF_DIR)/ALL.chrCHR.phase3.combined.sites.unfiltered.vcf.gz \$(REF_DIR)/chrCHR.filtered.sites.vcf.gz" >> ${CG_conf_file}
	
	### run all with defaults
	#gotcloud snpcall --outdir $MDIR/processed_data/GotCloud/${gene} --bam_list $MDIR/scripts/BAM_list_${gene}.txt --bamprefix $BAMLOC --numjobs 35 --maxlocaljobs 35 --verbose --chrs $chr --region $chr:$lo-$up
	
	 ### run with the config above
	gotcloud snpcall --conf ${CG_conf_file} --outdir $MDIR/processed_data/GotCloud/${gene} --bam_list $MDIR/scripts/BAM_list_${gene}.txt --bamprefix $BAMLOC --numjobs 35 --maxlocaljobs 35 --verbose --chrs $chr --region $chr:$lo-$up
done

### Delete everything in the output
# for gene in GCKR LPL APOA1 APOA5; do
# 	rm -rf $MDIR/processed_data/GotCloud/${gene}
# done

#### Not running this part for now
### run the genotype refinement (defaults seem fine for this step)
#gotcloud ldrefine --outdir $MDIR/processed_data/GotCloud/${gene} --bam_list $MDIR/scripts/BAM_list_${gene}.txt --bamprefix $BAMLOC --numjobs 30 --maxlocaljobs 30 --verbose --chrs $chr --region $chr:$lo-$up


####
#### Step 3. Annotate vcf file and convert to PLINK
####


PLINK=/path/to/plink/1.09
FDIR=$MDIR/processed_data

# gene=LPL; chr=8; 
for gene in GCKR LPL APOA1 APOA5; do
	
	echo $gene
	## Define region per gene. Note that all these coordinates are in GRCh37.p13, i.e. hg19
	if [ "$gene" == "GCKR" ]; then chr=2; lo=$((27719470-5000)); up=$((27746556+5000)); 
	elif [ "$gene" == "LPL" ]; then chr=8; lo=$((19796582-5000)); up=$((19824770+5000)); 
	elif [ "$gene" == "APOA1" ]; then chr=11; lo=$((116706467-5000)); up=$((116708338+5000)); 
	elif [ "$gene" == "APOA5" ]; then chr=11; lo=$((116660086-5000)); up=$((116663136+5000)); 
	else echo "Not a gene we selected"; 
	fi
	
	## Other filtered vcf (without PASS)
	#cd $MDIR/processed_data/GotCloud/${gene}/vcfs/chr${chr}/
	#$PLINK/plink --vcf chr${chr}.filtered.vcf.gz --geno 0.01 --mac 3 --make-bed --out $FDIR/${gene}_GotCloud_filtered
	
	### I believe files *filtered.PASS.vcf.gz only keep variants with INFO equats PASS (SNPs only, no INDELS nor anything else)
	cd $MDIR/processed_data/GotCloud/${gene}/split/chr${chr}/
	
	### annotate using bcftools with the dbSNPb37 from gotcloud reference file
	### but first fix some issues with the headers PL with the wrong number type, for instance
	gunzip -c chr${chr}.filtered.PASS.vcf.gz | head -1000 | grep -e "^#" > fixed.hdr ## advantage vcf doesn't need to be indexed
	sed -i 's/ID=PL,Number=./ID=PL,Number=G/g' fixed.hdr
	# In addition, some INFO fields are needed in tha header so, I'm putting them in a new file
	echo '##INFO=<ID=AZ,Number=1,Type=String,Description="Description for AZ info">' > add.hdr
	echo '##INFO=<ID=FIC,Number=1,Type=String,Description="Description for FIC info">' >> add.hdr
	echo '##INFO=<ID=SLRT,Number=1,Type=String,Description="Description for SLRT info">' >> add.hdr
	echo '##INFO=<ID=LBS,Number=1,Type=String,Description="Description for LBS info">' >> add.hdr
	echo '##INFO=<ID=OBS,Number=1,Type=String,Description="Description for OBS info">' >> add.hdr
	echo '##INFO=<ID=LQR,Number=1,Type=String,Description="Description for LQR info">' >> add.hdr
	echo '##INFO=<ID=SVM,Number=1,Type=String,Description="Description for SVM info">' >> add.hdr
	
	bcftools reheader -h fixed.hdr -o chr${chr}.filtered.PASS_PLfixed.vcf.gz chr${chr}.filtered.PASS.vcf.gz
	
	#index since required by annotate
	tabix -f chr${chr}.filtered.PASS_PLfixed.vcf.gz
	
	#this step performs the actual annotation 
	bcftools annotate -a $GCDIR/gotcloud.ref/dbsnp_142.b37.vcf.gz -c ID -x INFO,FILTER,FORMAT -h add.hdr -Oz -o chr${chr}.filtered.PASS_dbSNPb37anotated.vcf.gz -r $chr --threads 30 chr${chr}.filtered.PASS_PLfixed.vcf.gz
	
	#rm chr${chr}.filtered.PASS_PLfixed.vcf.gz chr${chr}.filtered.PASS_PLfixed.vcf.gz.tbi
	
	#### convert to plink
	$PLINK/plink --vcf chr${chr}.filtered.PASS_dbSNPb37anotated.vcf.gz --keep-allele-order --chr $chr --from-bp $lo --to-bp $up --set-missing-var-ids '@:#:$1:$2' --geno 0.01 --mac 3 --make-bed --out $FDIR/${gene}_GotCloud_filteredPASS
done

### change the ids to match the chip data as well as remove the 12 subjects without genotype
cd $FDIR
for gene in GCKR LPL APOA1 APOA5; do
	echo $gene
	
	$PLINK/plink --bfile $FDIR/${gene}_GotCloud_filteredPASS --remove IDs_to_remove_seqdatan=12.txt --update-ids Convert_IDs_seq_to_geno_n=4511.txt --make-bed --out ${gene}_GotCloud_filteredPASS_flt
	
	$PLINK/plink --bfile ${gene}_GotCloud_filteredPASS_flt --freqx --out ${gene}_GotCloud_filteredPASS_flt1
	$PLINK/plink --bfile ${gene}_GotCloud_filteredPASS_flt --freq --out ${gene}_GotCloud_filteredPASS_flt1
done

	
	#### Run using slurm
	#gotcloud snpcall --conf ${CG_conf_file} --outdir $MDIR/processed_data/GotCloud/${gene} --bam_list $MDIR/scripts/BAM_list_${gene}_test3.txt --bamprefix $BAMLOC --numjobs 100 --batchtype slurm --batchopts "-p biostat -n 100 --mem-per-cpu=1000M -t 7:00:00" --verbose --chrs $chr --region $chr:$lo-$up
	
	#gotcloud ldrefine --outdir $MDIR/processed_data/GotCloud/${gene} --bam_list $MDIR/scripts/BAM_list_${gene}_test3.txt --bamprefix $BAMLOC --numjobs 100 --batchtype slurm --batchopts "-p biostat -n 100 --mem-per-cpu=1000M -t 1:00:00" --verbose --chrs $chr --region $chr:$lo-$up


