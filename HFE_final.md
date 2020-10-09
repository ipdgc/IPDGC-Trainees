#!/bin/bash

####### GWAS #######

module load plink samtools rvtests R annovar

```
mkdir /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE
cd /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE
mkdir hardcallsNoNeuroX hardcallsNoNeuroX/bin hardcallsNoNeuroX/freq hardcallsNoNeuroX/score hardcallsNoNeuroX/burden hardcallsNoNeuroX/vcf hardcallsNoNeuroX/annotation hardcallsNoNeuroX/logistic
```
# ==================== 1. Subset PLINK Binaries + Convert to VCFs ====================

HFE positions on hg19

```
plink --bfile /data/LNG/saraB/HARDCALLS_PD_september_2018_no_cousins --remove-fam /data/LNG/saraB/NeuroX.fID.txt --geno 0.05 --chr 6 --from-bp 26087281 --to-bp 26098343 --make-bed --out hardcallsNoNeuroX/bin/HFE.GWAS
plink --bfile hardcallsNoNeuroX/bin/HFE.GWAS --recode vcf --out hardcallsNoNeuroX/vcf/HFE.GWAS
```
```
cd hardcallsNoNeuroX/vcf
bgzip HFE.GWAS.vcf
tabix -f -p vcf HFE.GWAS.vcf.gz
```
HFE_100KB positions on hg19

```
plink --bfile /data/LNG/saraB/HARDCALLS_PD_september_2018_no_cousins --remove-fam /data/LNG/saraB/NeuroX.fID.txt --geno 0.05 --chr 6 --from-bp 25987281 --to-bp 26198343 --make-bed --out hardcallsNoNeuroX/bin/HFE_100KB.GWAS
plink --bfile hardcallsNoNeuroX/bin/HFE_100KB.GWAS --recode vcf --out hardcallsNoNeuroX/vcf/HFE_100KB.GWAS
```

```
cd hardcallsNoNeuroX/vcf
bgzip HFE_100KB.GWAS.vcf
tabix -f -p vcf HFE_100KB.GWAS.vcf.gz
```

# ==================== 2. Fisher exact test and logistic regression ====================

```
cd /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE

plink --bfile hardcallsNoNeuroX/bin/HFE.GWAS --chr 6 --from-bp 26087281 --to-bp 26098343 --fisher --out hardcallsNoNeuroX/freq/HFE
plink --bfile hardcallsNoNeuroX/bin/HFE.GWAS --chr 6 --from-bp 26087281 --to-bp 26098343 --covar /data/LNG/saraB/IPDGC_all_samples_covariates.tab --covar-name sex,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --logistic --out hardcallsNoNeuroX/logistic/HFE_ci95 --ci 0.95
```
```
cd /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE

plink --bfile hardcallsNoNeuroX/bin/HFE_100KB.GWAS --chr 6 --from-bp 25987281 --to-bp 26198343 --fisher --out hardcallsNoNeuroX/freq/HFE_100KB
plink --bfile hardcallsNoNeuroX/bin/HFE_100KB.GWAS --chr 6 --from-bp 25987281 --to-bp 26198343 --covar /data/LNG/saraB/IPDGC_all_samples_covariates.tab --covar-name sex,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --logistic --out hardcallsNoNeuroX/logistic/HFE_100KB_ci95 --ci 0.95
```
# ==================== 3. Fisher exact test and logistic regression: Males versus Females ====================

# HFE

```
cd /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE

plink --bfile hardcallsNoNeuroX/bin/HFE.GWAS --chr 6 --from-bp 26087281 --to-bp 26098343 --filter-males --fisher --out hardcallsNoNeuroX/freq/HFE_males_ci95 --ci 0.95
plink --bfile hardcallsNoNeuroX/bin/HFE.GWAS --chr 6 --from-bp 26087281 --to-bp 26098343 --filter-males --covar /data/LNG/saraB/IPDGC_all_samples_covariates.tab --covar-name sex,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --logistic --out hardcallsNoNeuroX/logistic/HFE_males_ci95 --ci 0.95
plink --bfile hardcallsNoNeuroX/bin/HFE.GWAS --chr 6 --from-bp 26087281 --to-bp 26098343 --filter-females --fisher --out hardcallsNoNeuroX/freq/HFE_females_ci95 --ci 0.95
plink --bfile hardcallsNoNeuroX/bin/HFE.GWAS --chr 6 --from-bp 26087281 --to-bp 26098343 --filter-females --covar /data/LNG/saraB/IPDGC_all_samples_covariates.tab --covar-name sex,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --logistic --out hardcallsNoNeuroX/logistic/HFE_females_ci95 --ci 0.95
```
#HFE + 100kb

```
cd /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE

plink --bfile hardcallsNoNeuroX/bin/HFE_100KB.GWAS --chr 6 --from-bp 25987281 --to-bp 26198343 --filter-males --fisher --out hardcallsNoNeuroX/freq/HFE_100KB_males_ci95 --ci 0.95
plink --bfile hardcallsNoNeuroX/bin/HFE_100KB.GWAS --chr 6 --from-bp 25987281 --to-bp 26198343 --filter-males --covar /data/LNG/saraB/IPDGC_all_samples_covariates.tab --covar-name sex,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --logistic --out hardcallsNoNeuroX/logistic/HFE_100KB_males_ci95 --ci 0.95
plink --bfile hardcallsNoNeuroX/bin/HFE_100KB.GWAS --chr 6 --from-bp 25987281 --to-bp 26198343 --filter-females --fisher --out hardcallsNoNeuroX/freq/HFE_100KB_females_ci95 --ci 0.95
plink --bfile hardcallsNoNeuroX/bin/HFE_100KB.GWAS --chr 6 --from-bp 25987281 --to-bp 26198343 --filter-females --covar /data/LNG/saraB/IPDGC_all_samples_covariates.tab --covar-name sex,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --logistic --out hardcallsNoNeuroX/logistic/HFE_100KB_females_ci95 --ci 0.95
```

# ==================== 4. Annotation ====================

```
mkdir /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/annovar/

# Downloading databases 
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
# annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/ 
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 humandb/ 
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a humandb/
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_genome humandb/
```

```
module load annovar

cd /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/vcf

table_annovar.pl HFE.GWAS.vcf.gz /data/LNG/saraB/annovar/humandb/ -buildver hg19 \
--thread 16 \
-out HFE.annovar \
-remove -protocol avsnp147,refGene,ensGene,gnomad211_genome \
-operation f,g,g,f \
-nastring . \
-vcfinput

head -1 HFE.annovar.hg19_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct HFE.annovar.hg19_multianno.txt > HFE.GWAS.trimmed.annotation.txt
```

```
module load annovar

cd /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/vcf

table_annovar.pl HFE_100KB.GWAS.vcf.gz /data/LNG/saraB/annovar/humandb/ -buildver hg19 \
--thread 16 \
-out HFE_100KB.annovar \
-remove -protocol avsnp147,refGene,ensGene,gnomad211_genome \
-operation f,g,g,f \
-nastring . \
-vcfinput

head -1 HFE_100KB.annovar.hg19_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct HFE_100KB.annovar.hg19_multianno.txt > HFE_100KB.GWAS.trimmed.annotation.txt
```
# ==================== 5. Burden analysis ====================

#### Generating list of coding variants

```
awk '$7=="exonic" {print}' HFE.GWAS.trimmed.annotation.txt > HFE.GWAS.trimmed.annotation.coding.variants.txt
awk '{print $1" "$2" "$2" "$7}' HFE.GWAS.trimmed.annotation.coding.variants.txt > HFE.GWAS.trimmed.annotation.coding.variants.SNPs.txt

plink --bfile  /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/bin/HFE.GWAS --extract range HFE.GWAS.trimmed.annotation.coding.variants.SNPs.txt --recode vcf  --out HFE_CODING_GWAS

bgzip HFE_CODING_GWAS.vcf
tabix -f -p vcf HFE_CODING_GWAS.vcf.gz
```

```
awk '$7=="exonic" {print}' HFE_100KB.GWAS.trimmed.annotation.txt > HFE_100KB.GWAS.trimmed.annotation.coding.variants.txt
awk '{print $1" "$2" "$2" "$7}' HFE_100KB.GWAS.trimmed.annotation.coding.variants.txt > HFE_100KB.GWAS.trimmed.annotation.coding.variants.SNPs.txt

plink --bfile  /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/bin/HFE_100KB.GWAS --extract range HFE_100KB.GWAS.trimmed.annotation.coding.variants.SNPs.txt --recode vcf  --out HFE_100KB_CODING_GWAS

bgzip HFE_100KB_CODING_GWAS.vcf
tabix -f -p vcf HFE_100KB_CODING_GWAS.vcf.gz
```

#### MAF < 0.03

#### ALL VARIANTS 

```
rvtest --noweb --inVcf /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/vcf/HFE.GWAS.vcf.gz --pheno /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile /data/LNG/saraB/refFlat_hg19.txt --freqUpper 0.03 --out /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/burden/BURDEN.HFE.GWAS.maf003
rvtest --noweb --inVcf /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/vcf/HFE_100KB.GWAS.vcf.gz --pheno /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile /data/LNG/saraB/refFlat_hg19.txt --freqUpper 0.03 --out /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/burden/BURDEN.HFE_100KB.GWAS.maf003
```
#### CODING VARIANTS 

```
rvtest --noweb --inVcf /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/vcf/HFE_CODING_GWAS.vcf.gz --pheno /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile /data/LNG/saraB/refFlat_hg19.txt --freqUpper 0.03 --out /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/burden/BURDEN.HFE.GWAS.maf003_CODING
rvtest --noweb --inVcf /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/vcf/HFE_100KB_CODING_GWAS.vcf.gz --pheno /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile /data/LNG/saraB/refFlat_hg19.txt --freqUpper 0.03 --out /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/burden/BURDEN.HFE_100KB.GWAS.maf003_CODING
```
#### MAF < 0.01

#### ALL VARIANTS 

```
rvtest --noweb --inVcf /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/vcf/HFE.GWAS.vcf.gz --pheno /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile /data/LNG/saraB/refFlat_hg19.txt --freqUpper 0.01 --out /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/burden/BURDEN.HFE.GWAS.maf001
rvtest --noweb --inVcf /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/vcf/HFE_100KB.GWAS.vcf.gz --pheno /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile /data/LNG/saraB/refFlat_hg19.txt --freqUpper 0.01 --out /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/burden/BURDEN.HFE_100KB.GWAS.maf001
```
#### CODING VARIANTS 

```
rvtest --noweb --inVcf /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/vcf/HFE_CODING_GWAS.vcf.gz --pheno /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile /data/LNG/saraB/refFlat_hg19.txt --freqUpper 0.01 --out /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/burden/BURDEN.HFE.GWAS.maf001_CODING
rvtest --noweb --inVcf /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/vcf/HFE_100KB_CODING_GWAS.vcf.gz --pheno /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile /data/LNG/saraB/refFlat_hg19.txt --freqUpper 0.01 --out /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/hardcallsNoNeuroX/burden/BURDEN.HFE_100KB.GWAS.maf001_CODING
```

# ==================== ####### WGS ####### ==================== #
```
cd /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/HFE/WGS/
```

# ==================== 6. Fisher exact test and logistic regression ====================

HFE and HFE_100KB positions on hg38

```
plink --bfile HFE_WGS_AMP_PD --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --make-bed --out HFE_WGS_AMP_PD_pheno
plink --bfile HFE_WGS_AMP_PD_pheno --update-sex /data/LNG/saraB/AMP_PD/toupdatesex.txt --make-bed --out HFE_WGS_AMP_PD_pheno_sex
plink --bfile HFE_WGS_AMP_PD_pheno_sex --fisher --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --out HFE_WGS_AMP_PD 
plink --bfile HFE_WGS_AMP_PD_pheno_sex --logistic --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX, AGE, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 --out HFE_WGS_AMP_PD --ci 0.95
```

```
plink --bfile HFE_100kb_WGS_AMP_PD --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --make-bed --out HFE_100kb_WGS_AMP_PD_pheno
plink --bfile HFE_100kb_WGS_AMP_PD_pheno --update-sex /data/LNG/saraB/AMP_PD/toupdatesex.txt --make-bed --out HFE_100kb_WGS_AMP_PD_pheno_sex
plink --bfile HFE_100kb_WGS_AMP_PD_pheno_sex --fisher --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --out HFE_100kb_WGS_AMP_PD
plink --bfile HFE_100kb_WGS_AMP_PD_pheno_sex --logistic --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX, AGE, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 --out HFE_100kb_WGS_AMP_PD --ci 0.95
```

# ==================== 7. Males versus Females - Fisher exact test and logistic regression ====================

## HFE
```
plink --bfile HFE_WGS_AMP_PD_pheno_sex --fisher --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --out HFE_WGS_AMP_PD_females --filter-females
plink --bfile HFE_WGS_AMP_PD_pheno_sex --logistic --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX, AGE, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 --out HFE_WGS_AMP_PD_females --ci 0.95 --filter-females
plink --bfile HFE_WGS_AMP_PD_pheno_sex --fisher --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --out HFE_WGS_AMP_PD_males --filter-males
plink --bfile HFE_WGS_AMP_PD_pheno_sex --logistic --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX, AGE, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 --out HFE_WGS_AMP_PD_males --ci 0.95 --filter-males
```

## HFE + 100KB
```
plink --bfile HFE_100kb_WGS_AMP_PD_pheno_sex --fisher --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --out HFE_100kb_WGS_AMP_PD_females --filter-females
plink --bfile HFE_100kb_WGS_AMP_PD_pheno_sex --logistic --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX, AGE, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 --out HFE_100kb_WGS_AMP_PD_females --ci 0.95 --filter-females
plink --bfile HFE_100kb_WGS_AMP_PD_pheno_sex --fisher --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --out HFE_100kb_WGS_AMP_PD_males --filter-males
plink --bfile HFE_100kb_WGS_AMP_PD_pheno_sex --logistic --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX, AGE, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 --out HFE_100kb_WGS_AMP_PD_males --ci 0.95 --filter-males
```
# ==================== 8. Annotation ====================

```
plink --bfile HFE_WGS_AMP_PD_pheno_sex --recode 'vcf-fid' --out HFE_WGS_AMP_PD_pheno_sex
bgzip HFE_WGS_AMP_PD_pheno_sex.vcf
tabix -f -p vcf HFE_WGS_AMP_PD_pheno_sex.vcf.gz

table_annovar.pl HFE_WGS_AMP_PD_pheno_sex.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
--thread 16 \
-out HFE_WGS.annovar \
-remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
-operation g,f,f,f \
-nastring . \
-vcfinput

head -1 HFE_WGS.annovar.hg38_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct HFE_WGS.annovar.hg38_multianno.txt > HFE.WGS.trimmed.annotation.txt
```

```
plink --bfile HFE_100kb_WGS_AMP_PD_pheno_sex --recode 'vcf-fid' --out HFE_100kb_WGS_AMP_PD_pheno_sex
bgzip HFE_100kb_WGS_AMP_PD_pheno_sex.vcf
tabix -f -p vcf HFE_100kb_WGS_AMP_PD_pheno_sex.vcf.gz

table_annovar.pl HFE_100kb_WGS_AMP_PD_pheno_sex.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
--thread 16 \
-out HFE_100kb_WGS.annovar \
-remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
-operation g,f,f,f \
-nastring . \
-vcfinput

head -1 HFE_100kb_WGS.annovar.hg38_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct HFE_100kb_WGS.annovar.hg38_multianno.txt > HFE_100kb_WGS.trimmed.annotation.txt
```

# ==================== 9. Burden analysis ====================

#### Generating list of coding variants

```
awk '$6=="exonic" {print}' HFE.WGS.trimmed.annotation.txt > HFE.trimmed.annotation.coding.variants.WGS.txt
awk '{print $1" "$2" "$2" "$7}' HFE.trimmed.annotation.coding.variants.WGS.txt > HFE.trimmed.annotation.coding.variants.WGS.SNPs.txt

plink --bfile HFE_WGS_AMP_PD_pheno_sex --extract range HFE.trimmed.annotation.coding.variants.WGS.txt --recode 'vcf-fid' --out HFE_CODING_AMP

bgzip HFE_CODING_AMP.vcf
tabix -f -p vcf HFE_CODING_AMP.vcf.gz
```

```
awk '$6=="exonic" {print}' HFE_100kb_WGS.trimmed.annotation.txt > HFE_100kb.trimmed.annotation.coding.variants.WGS.txt
awk '{print $1" "$2" "$2" "$7}' HFE_100kb.trimmed.annotation.coding.variants.WGS.txt > HFE_100kb.trimmed.annotation.coding.variants.WGS.SNPs.txt

plink --bfile HFE_100kb_WGS_AMP_PD_pheno_sex --extract range HFE_100kb.trimmed.annotation.coding.variants.WGS.txt --recode 'vcf-fid' --out HFE_100kb_CODING_AMP

bgzip HFE_100kb_CODING_AMP.vcf
tabix -f -p vcf HFE_100kb_CODING_AMP.vcf.gz
```

#### MAF < 0.03

#### ALL VARIANTS 

```
rvtest --noweb --inVcf HFE_WGS_AMP_PD_pheno_sex.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 --out AMP_PD_BURDEN_WGS.maf003
rvtest --noweb --inVcf HFE_100kb_WGS_AMP_PD_pheno_sex.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 --out AMP_PD_BURDEN.HFE100kb_WGS.maf003
```
#### CODING VARIANTS 

```
rvtest --noweb --inVcf HFE_CODING_AMP.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 --out AMP_PD_BURDEN.HFE.WGS.maf003_CODING
rvtest --noweb --inVcf HFE_100kb_CODING_AMP.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 --out AMP_PD_BURDEN.HFE100kb_WGS.maf003_CODING
```

#### MAF < 0.01

#### ALL VARIANTS 

```
rvtest --noweb --inVcf HFE_WGS_AMP_PD_pheno_sex.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.01 --out AMP_PD_BURDEN_WGS.maf001
rvtest --noweb --inVcf HFE_100kb_WGS_AMP_PD_pheno_sex.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.01 --out AMP_PD_BURDEN.HFE100kb_WGS.maf001
```

#### CODING VARIANTS 

```
rvtest --noweb --inVcf HFE_CODING_AMP.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.01 --out AMP_PD_BURDEN.HFE.WGS.maf001_CODING
rvtest --noweb --inVcf HFE_100kb_CODING_AMP.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.01 --out AMP_PD_BURDEN.HFE100kb_WGS.maf001_CODING
```