#!/bin/env bash 

## DESCRIPTION
	# Author(s): Sara Bandres-Ciga, Mary B. Makarious, Monica Diez-Fairen
	# Project: LIN28A Letter (IPDGC)
	# PI(s): Mike Nalls, Andrew Singleton
	# Collaborators: Cornelis Blauwendraat 
	# Date Started: 13.04.2020
		# Date Last Updated: 21.04.2020
		# Update: Running Wald and score tests  

## PROPOSED WORKFLOW 
	# 0. Getting Started 
		# Create the covariate file, sex info file, pheno info file
		# Update the AMP files
	# 1. Extract the LIN28A region 
		# Generate --mac 3 for statistics, --mac 1 for burdens 
	# 2. Running with --model and --assoc in PLINK 
	# 3. Running with --logistic and --assoc in PLINK
	# 4. Generate VCFs
	# 5. Annotate via ANNOVAR
	# 6. Burden Analysis via RVTests 
	# 7. Single variant Wald and score tests via RVTests

#########################################################################################################################
##### 0. GETTING STARTED ################################################################################################
#########################################################################################################################

# Initializing some variables 
## REMOVED paths to files 
COV_NAMES="SEX,AGE,FAMILY_HISTORY,EDUCATION,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,STUDY_BioFIND,STUDY_PDBP,STUDY_PPMI"

# Load the necessary modules 
module load plink #v1.9
module load samtools #v1.9
module load annovar #v2018-04-16
module load rvtests #v2.1.0


# AMP-PD Info
	# BioFind, PDBP, and PPMI cohorts 
		# --geno 0.05 
		# European ancestry 
	# 28195229 variants
	# 2927 people total - 
		# 1647 cases, 1050 controls (230 missing phenotype)

#########################################################################################################################
##### 1. EXTRACT THE LIN28A REGION ######################################################################################
#########################################################################################################################

# Create a working directory for the LIN28A paper 
# Change into directory
cd $WORKING_DIR

# Extract the LIN28A region
	# Cornelis recommends --mac 3 for statistics and --mac 1 for burdens 
# --mac 1
plink --bfile $AMP_PLINK \
--chr 1 --from-bp 26410817 --to-bp 26429728 \
--mac 1 \
--make-bed --out $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1
	# Done 21.04.2020
	# PLINK output:
		# Total genotyping rate is 0.999946.
		# 108 variants and 2927 people pass filters and QC.
		# Among remaining phenotypes, 1647 are cases and 1050 are controls. (230 phenotypes are missing.)

# --mac 3 
plink --bfile $AMP_PLINK \
--chr 1 --from-bp 26410817 --to-bp 26429728 \
--mac 3 \
--make-bed --out $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3 

# wc -l LIN28A_WGS_AMP_PD_mac3.bim
	# 83 LIN28A_WGS_AMP_PD_mac3.bim
	# 83 variants total pass --mac 3 

#########################################################################################################################
##### 2. RUNNING WITH --MODEL AND --ASSOC IN PLINK ######################################################################
#########################################################################################################################

# We run --model and --assoc to get the frequencies of the SNPs and distribution of alleles of the dataset

# Running with --model
plink --bfile $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3 --model \
--out LIN28A_WGS_AMP_PD_mac3_noCov
	# File generated @ LIN28A_WGS_AMP_PD_mac3.model
		# --model does not take covariates 
		# Shows distribution of hom/het/hom 
	# Done 21.04.2020

# Running --assoc with no covariates 
plink --bfile $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3 --assoc \
--pheno $PHENO \
--ci 0.95 \
--out LIN28A_WGS_AMP_PD_mac3_noCov
	# File generated @ LIN28A_WGS_AMP_PD_mac3_noCov.assoc
		# --assoc is a quick way to get frequencies 
		# Does not take covariates 
	# Done 21.04.2020


#########################################################################################################################
##### 3. RUNNING WITH --LOGISTIC + --FISHER IN PLINK ####################################################################
#########################################################################################################################

# Running --logistic with covariates 
plink --bfile $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3 --logistic \
--pheno $PHENO \
--covar $COV_FILE \
--covar-name SEX,AGE,FAMILY_HISTORY,EDUCATION,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--ci 0.95 \
--hide-covar \
--out LIN28A_WGS_AMP_PD_mac3_wCovs
	# File generated @ LIN28A_WGS_AMP_PD_mac3_wCovs.assoc.logistic
		# --logistic does take covariates but wasn't run here 
		# --hide-covar flag to remove each and every individual covariate test from the final output
		# Adding study covariates with PLINK did not work... so were removed
			# But did work for RVTests, see section 7 
	# Done 21.04.2020

# Running --fisher with covariates 
plink --bfile $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3 --assoc fisher \
--pheno $PHENO \
--ci 0.95 \
--out LIN28A_WGS_AMP_PD_mac3 
	# File generated @ LIN28A_WGS_AMP_PD_mac3.assoc.fisher
		# --assoc fisher does not take covariates 
	# Done 21.04.2020


#########################################################################################################################
##### 4. GENERATE VCFS ##################################################################################################
#########################################################################################################################

# Generating a VCF file for --mac 1
plink --bfile $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1 \
--mac 1 \
--recode 'vcf-fid' --out $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1
	# Done 21.04.2020

# Generating a VCF file for --mac 3
plink --bfile $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3 \
--mac 3 \
--recode 'vcf-fid' --out $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3
	# Done 21.04.2020

# Zipping and tabix-ing the VCF for --mac 1
bgzip $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1.vcf
tabix -f -p vcf $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1.vcf.gz
	# Done 21.04.2020

# Zipping and tabix-ing the VCF for --mac 3
bgzip $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3.vcf
tabix -f -p vcf $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3.vcf.gz
	# Done 21.04.2020

#########################################################################################################################
##### 5. ANNOTATE VIA ANNOVAR ###########################################################################################
#########################################################################################################################

# Annotate with ANNOVAR
	# For burdens, only need --mac 1 
table_annovar.pl $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1.vcf.gz \
$ANNOVAR_DATA/hg38 --thread 20 -buildver hg38 \
-out $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1.annovar \
-arg '-splicing 15',,, \
-remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
-operation g,f,f,f -nastring . -vcfinput -polish
	# Done 21.04.2020
# Information: 
	# wc -l LIN28A_WGS_AMP_PD.annovar.hg38_multianno.txt
	# 135 LIN28A_WGS_AMP_PD.annovar.hg38_multianno.txt
	# Matches variant numbers above 

# Remove the VCF 
rm $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1.annovar.hg38_multianno.vcf

# Trim the annotation file
head -1 LIN28A_WGS_AMP_PD_mac1.annovar.hg38_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct LIN28A_WGS_AMP_PD_mac1.annovar.hg38_multianno.txt > LIN28A.trimmed.annotation.txt

#########################################################################################################################
##### 6. BURDEN ANALYSIS VIA RVTESTS ####################################################################################
#########################################################################################################################

# Generate CODING.txt - the range of exonic regions
	# For burdens, only need --mac 1 
	# synonymous + non-synonymous 
	# From trimmed file 
grep "exonic" $WORKING_DIR/LIN28A.trimmed.annotation.txt > $WORKING_DIR/CODING_regions.txt
cut -f 1,2,3,7 $WORKING_DIR/CODING_regions.txt > CODING.txt 
	# Checked to see if it matches Sara's file -- it does! 
# Information: 
	# wc -l CODING.txt
	# 4 CODING.txt
	# head CODING.txt
		# 1	26411441	26411441	LIN28A
		# 1	26425455	26425455	LIN28A
		# 1	26426394	26426394	LIN28A
		# 1	26426443	26426443	LIN28A

# Prep the files and generate PLINK with only coding variants 
plink --bfile $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1 --extract range CODING.txt --make-bed --out $WORKING_DIR/LIN28A_exonic_AMP_PD 
	# Files generated @ LIN28A_exonic_AMP_PD.*
	# Done 21.04.2020

# Create the VCF with only coding variants 
plink --bfile $WORKING_DIR/LIN28A_exonic_AMP_PD --recode 'vcf-fid' --out $WORKING_DIR/LIN28A_exonic_AMP_PD
	# Done 21.04.2020

# bgzip and tabix
bgzip $WORKING_DIR/LIN28A_exonic_AMP_PD.vcf
tabix -f -p vcf $WORKING_DIR/LIN28A_exonic_AMP_PD.vcf.gz
	# Done 21.04.2020

## Run RVTests
# NO MAF THRESHOLD

## ALL VARIANTS ##
rvtest --noweb --hide-covar --inVcf $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1.vcf.gz \
--pheno $COV_FILE \
--pheno-name PHENO \
--covar $COV_FILE \
--covar-name $COV_NAMES \
--kernel skat,skato --geneFile /data/LNG/makariousmb/refFlat_hg38.txt \
--out AMP_PD_BURDEN.LIN28A.NOFREQTHRESHOLD
	# Files generated @ AMP_PD_BURDEN.LIN28A.NOFREQTHRESHOLD.*.assoc
	# Done 21.04.2020

## CODING VARIANTS ##
rvtest --noweb --hide-covar --inVcf $WORKING_DIR/LIN28A_exonic_AMP_PD.vcf.gz \
--pheno $COV_FILE --covar $COV_FILE \
--pheno-name PHENO \
--covar-name $COV_NAMES \
--kernel skat,skato --geneFile /data/LNG/makariousmb/refFlat_hg38.txt \
--out AMP_PD_BURDEN.LIN28A.NOFREQTHRESHOLD_CODING
	# Files generated @ AMP_PD_BURDEN.LIN28A.NOFREQTHRESHOLD_CODING.*.assoc
	# Done 21.04.2020

# MAF < 0.01

## ALL VARIANTS ##
rvtest --noweb --hide-covar --inVcf $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1.vcf.gz \
--pheno $COV_FILE \
--pheno-name PHENO \
--covar $COV_FILE \
--covar-name $COV_NAMES \
--kernel skat,skato --geneFile /data/LNG/makariousmb/refFlat_hg38.txt \
--freqUpper 0.01 --out AMP_PD_BURDEN.LIN28A.maf001
	# Files generated @ AMP_PD_BURDEN.LIN28A.maf001.*.assoc
	# Done 21.04.2020

## CODING VARIANTS ##
rvtest --noweb --hide-covar --inVcf $WORKING_DIR/LIN28A_exonic_AMP_PD.vcf.gz \
--pheno $COV_FILE \
--pheno-name PHENO \
--covar $COV_FILE \
--covar-name $COV_NAMES \
--kernel skat,skato --geneFile /data/LNG/makariousmb/refFlat_hg38.txt \
--freqUpper 0.01 --out AMP_PD_BURDEN.LIN28A.maf001_CODING
	# Files generated @ AMP_PD_BURDEN.LIN28A.maf001_CODING.*.assoc
	# Done 21.04.2020

# MAF < 0.03

## ALL VARIANTS ##
rvtest --noweb --hide-covar --inVcf $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1.vcf.gz \
--pheno $COV_FILE \
--pheno-name PHENO \
--covar $COV_FILE \
--covar-name $COV_NAMES \
--kernel skat,skato --geneFile /data/LNG/makariousmb/refFlat_hg38.txt \
--freqUpper 0.03 --out AMP_PD_BURDEN.LIN28A.maf003
	# Files generated @ AMP_PD_BURDEN.LIN28A.maf003.*.assoc
	# Done 21.04.2020

## CODING VARIANTS ##
rvtest --noweb --inVcf $WORKING_DIR/LIN28A_exonic_AMP_PD.vcf.gz \
--pheno $COV_FILE --covar $COV_FILE \
--pheno-name PHENO \
--covar-name $COV_NAMES \
--kernel skat,skato \
--geneFile /data/LNG/makariousmb/refFlat_hg38.txt \
--freqUpper 0.03 --out AMP_PD_BURDEN.LIN28A.maf003_CODING
	# Files generated @ AMP_PD_BURDEN.LIN28A.maf003_CODING.*.assoc
	# Done 21.04.2020

#########################################################################################################################
##### 7. SINGLE VARIANT WALD AND SCORE TESTS VIA RVTESTS ################################################################
#########################################################################################################################

# For statistics, need --mac 3
	# Reason being if there is one allele, then not enough information to do stats, it is assigned a 0 and results in an NA

# RVTests
	# Tests done on single variants done on ALL variants 
		# Wald = same as logistic in PLINK, better for common variants 
		# score = better for rare variants
rvtest --noweb --hide-covar --inVcf $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3.vcf.gz \
--pheno $COV_FILE \
--pheno-name PHENO \
--covar $COV_FILE \
--covar-name $COV_NAMES \
--single wald,score --geneFile /data/LNG/makariousmb/refFlat_hg38.txt \
--out AMP_PD_mac3.LIN28A.LOGISTIC.allPCs
	# Files generated at AMP_PD_mac3.LIN28A.LOGISTIC.allPCs.*.assoc
		# --hide-covar flag to remove each and every individual covariate test from the final output
		# The variant line reported includes all PCs 
	# Done 21.04.2020

