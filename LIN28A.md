# LIN28A Rare Variant Analysis

## Information 
#### DESCRIPTION
-  **Author(s):** Sara Bandres-Ciga, Mary B. Makarious, Monica Diez-Fairen
-  **Project:** LIN28A Letter (IPDGC)
- **PI(s):** Mike Nalls, Andrew Singleton
- **Collaborators:** Cornelis Blauwendraat 
- **Date Last Updated:** 21.04.2020

#### INSPIRATION
- **Working Title:** Assessment of LIN28A variants in Parkinson's disease
- **Aim:** To scrutinize whether LIN28A (LOF) mutations are causative for Parkinson's disease
- **Brief Description:** A loss-of-function LIN28A mutation has been reported to cause early onset Parkinson's disease. However, replication in a large cohort is needed. So we were inspired to look into the IPDGC + AMP PD cohorts
- **Original Publication:** [https://www.ncbi.nlm.nih.gov/pubmed/31750563](https://www.ncbi.nlm.nih.gov/pubmed/31750563)

#### PROPOSED WORKFLOW 
#### [0. Getting Started](#0)
#### [1. Extract the LIN28A region](#1)
- Generate `--mac 3` for statistics, `--mac 1` for burdens 
#### [2. Running with `--model` and `--assoc` in PLINK](#2)
#### [3. Running with `--logistic` and `--assoc` in PLINK](#3)
#### [4. Generate VCFs](#4)
#### [5. Annotate via ANNOVAR](#5)
#### [6. Burden Analysis via RVTests](#6)
#### [7. Single variant Wald and score tests via RVTests](#7)

<a id="0"></a>
## 0. GETTING STARTED 
```bash
# Initializing some variables 
## REMOVED paths to files 
COV_NAMES="SEX,AGE,FAMILY_HISTORY,EDUCATION,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,STUDY_BioFIND,STUDY_PDBP,STUDY_PPMI"
```

#### Load the necessary modules (done on Biowulf)
```bash
module load plink #v1.9
module load samtools #v1.9
module load annovar #v2018-04-16
module load rvtests #v2.1.0
```
### Information on AMP-PD WGS Used:
- BioFind, PDBP, and PPMI cohorts 
	- `--geno 0.05` 
	- European ancestry 
- 28195229 variants
	-  2927 people total - 
	- 1647 cases, 1050 controls (230 missing phenotype)

<a id="1"></a>

## 1. EXTRACT THE LIN28A REGION 

```bash
# Create a working directory for the LIN28A paper 
# Change into directory
cd $WORKING_DIR

# Extract the LIN28A region
	# These are hg38 co-ordinates 
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
```
<a id="2"></a>

## 2. RUNNING WITH --MODEL AND --ASSOC IN PLINK 

We run --model and --assoc to get the frequencies of the SNPs and distribution of alleles of the dataset

#### Running with `--model`
```bash
plink --bfile $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3 --model \
--out LIN28A_WGS_AMP_PD_mac3_noCov
	# File generated @ LIN28A_WGS_AMP_PD_mac3.model
		# --model does not take covariates 
		# Shows distribution of hom/het/hom 
	# Done 21.04.2020
```

#### Running --assoc with no covariates 
```bash
plink --bfile $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3 --assoc \
--pheno $PHENO \
--ci 0.95 \
--out LIN28A_WGS_AMP_PD_mac3_noCov
	# File generated @ LIN28A_WGS_AMP_PD_mac3_noCov.assoc
		# --assoc is a quick way to get frequencies 
		# Does not take covariates 
	# Done 21.04.2020
```
<a id="3"></a>
## 3. RUNNING WITH `--LOGISTIC` + `--FISHER` IN PLINK 

#### Running `--logistic` with covariates 
```bash
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
```

#### Running `--fisher` with covariates 
```bash
plink --bfile $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3 --assoc fisher \
--pheno $PHENO \
--ci 0.95 \
--out LIN28A_WGS_AMP_PD_mac3 
	# File generated @ LIN28A_WGS_AMP_PD_mac3.assoc.fisher
		# --assoc fisher does not take covariates 
	# Done 21.04.2020
```

<a id="4"></a>

## 4. GENERATE VCFS 

#### Generating a VCF file for `--mac 1`
```bash
plink --bfile $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1 \
--mac 1 \
--recode 'vcf-fid' --out $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1
	# Done 21.04.2020
```

#### Generating a VCF file for `--mac 3`
```bash
plink --bfile $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3 \
--mac 3 \
--recode 'vcf-fid' --out $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3
	# Done 21.04.2020
```

#### Zipping and tabix-ing the VCF for `--mac 1`
```bash
bgzip $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1.vcf
tabix -f -p vcf $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1.vcf.gz
	# Done 21.04.2020
```

#### Zipping and tabix-ing the VCF for `--mac 3`
```bash
bgzip $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3.vcf
tabix -f -p vcf $WORKING_DIR/LIN28A_WGS_AMP_PD_mac3.vcf.gz
	# Done 21.04.2020
```
<a id="5"></a>

## 5. ANNOTATE VIA ANNOVAR 

#### Annotate with ANNOVAR
```bash
# For burdens, only need --mac 1 

table_annovar.pl $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1.vcf.gz \
$ANNOVAR_DATA/hg38 --thread 20 -buildver hg38 \
-out $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1.annovar \
-arg '-splicing 15',,, \
-remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
-operation g,f,f,f -nastring . -vcfinput -polish

# Remove the VCF 
rm $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1.annovar.hg38_multianno.vcf

	# Done 21.04.2020
```

##### Information: 
- `wc -l LIN28A_WGS_AMP_PD.annovar.hg38_multianno.txt`
	- 135 LIN28A_WGS_AMP_PD.annovar.hg38_multianno.txt
	-  Matches variant numbers above 



#### Trim the Annotation File
```bash
head -1 LIN28A_WGS_AMP_PD_mac1.annovar.hg38_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct LIN28A_WGS_AMP_PD_mac1.annovar.hg38_multianno.txt > LIN28A.trimmed.annotation.txt
```

<a id="6"></a>

## 6. BURDEN ANALYSIS VIA RVTESTS 

#### Generate CODING.txt - the range of exonic regions
 For burdens, only need `--mac 1` 
 - synonymous + non-synonymous 
 -  From trimmed file 
```bash
grep "exonic" $WORKING_DIR/LIN28A.trimmed.annotation.txt > $WORKING_DIR/CODING_regions.txt
cut -f 1,2,3,7 $WORKING_DIR/CODING_regions.txt > CODING.txt 
# Information: 
	# wc -l CODING.txt
	# 4 CODING.txt
	# head CODING.txt
		# 1	26411441	26411441	LIN28A
		# 1	26425455	26425455	LIN28A
		# 1	26426394	26426394	LIN28A
		# 1	26426443	26426443	LIN28A
```

#### Prep the files and generate PLINK with only coding variants 
```bash
plink --bfile $WORKING_DIR/LIN28A_WGS_AMP_PD_mac1 --extract range CODING.txt --make-bed --out $WORKING_DIR/LIN28A_exonic_AMP_PD 
	# Files generated @ LIN28A_exonic_AMP_PD.*
	# Done 21.04.2020
```

#### Create the VCF with only coding variants 
```bash
plink --bfile $WORKING_DIR/LIN28A_exonic_AMP_PD --recode 'vcf-fid' --out $WORKING_DIR/LIN28A_exonic_AMP_PD
	# Done 21.04.2020

# bgzip and tabix
bgzip $WORKING_DIR/LIN28A_exonic_AMP_PD.vcf
tabix -f -p vcf $WORKING_DIR/LIN28A_exonic_AMP_PD.vcf.gz
	# Done 21.04.2020
```

### Run RVTests

#### NO MAF THRESHOLD on ALL and CODING ONLY variants 
```bash
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
```
#### MAF < 0.01 on ALL and CODING ONLY variants 
```bash
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
```
#### MAF < 0.03 on ALL and CODING ONLY variants 
```bash
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
```

<a id="7"></a>
## 7. SINGLE VARIANT WALD AND SCORE TESTS VIA RVTESTS 

**For statistics, need --mac 3**: Reason being if there is one allele, then not enough information to do stats, it is assigned a 0 and results in an NA

#### Single Variant Wald + Score Tests RVTests
```bash
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
```

