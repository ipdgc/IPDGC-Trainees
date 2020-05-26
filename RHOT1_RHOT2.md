# RHOT1

Teresa Perinan (Sevilla), Sara Bandres-Ciga (NIH)

## WGS analysis 

### Fisher exact test 

RHOT1 positions on hg38

module load plink samtools rvtests R annovar

```
plink --bfile RHOT1_WGS_AMP_PD --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --make-bed --out RHOT1_WGS_AMP_PD_pheno
plink --bfile RHOT1_WGS_AMP_PD_pheno --update-sex /data/LNG/saraB/AMP_PD/toupdatesex.txt --make-bed --out RHOT1_WGS_AMP_PD_pheno_sex
plink --bfile RHOT1_WGS_AMP_PD_pheno_sex --fisher --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX, AGE, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 --out RHOT1_WGS_AMP_PD --ci 0.95
```

### Annotation (ANNOVAR)

```
plink --bfile RHOT1_WGS_AMP_PD_pheno_sex --recode 'vcf-fid' --out RHOT1_WGS_AMP_PD_pheno_sex
bgzip RHOT1_WGS_AMP_PD_pheno_sex.vcf
tabix -f -p vcf RHOT1_WGS_AMP_PD_pheno_sex.vcf.gz

table_annovar.pl RHOT1_WGS_AMP_PD_pheno_sex.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
--thread 16 \
-out RHOT1_WGS.annovar \
-remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
-operation g,f,f,f \
-nastring . \
-vcfinput

head -1 RHOT1_WGS.annovar.hg38_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct RHOT1_WGS.annovar.hg38_multianno.txt > RHOT1.trimmed.annotation.txt
```

### Burden analysis 

#### Generating list of coding variants

```
awk '$6=="exonic" {print}' RHOT1.trimmed.annotation.txt > RHOT1.trimmed.annotation.coding.variants.WGS.txt
awk '{print $1" "$2" "$2" "$7}' RHOT1.trimmed.annotation.coding.variants.WGS.txt > RHOT1.trimmed.annotation.coding.variants.WGS.SNPs.txt

plink --bfile RHOT1_WGS_AMP_PD_pheno_sex --extract range RHOT1.trimmed.annotation.coding.variants.WGS.txt --recode 'vcf-fid' --out RHOT1_CODING_AMP

bgzip RHOT1_CODING_AMP.vcf
tabix -f -p vcf RHOT1_CODING_AMP.vcf.gz
```

#### MAF < 0.03

#### ALL VARIANTS 

```
rvtest --noweb --inVcf RHOT1_WGS_AMP_PD_pheno_sex.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 --out AMP_PD_BURDEN.RHOT1.maf003
```

#### CODING VARIANTS 

```
rvtest --noweb --inVcf RHOT1_CODING_AMP.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 --out AMP_PD_BURDEN.RHOT1.maf003_CODING
```

### Extract your gene(s) from PD summary stats without 23andMe data 

```
* listofgenes.txt will have the following format: chr start stop gene

R

library(data.table)
library(dplyr)
library(tidyr)
library(plyr)
data <- fread("NallsEtAl2019.tab", header = T)
data$CHR <- ldply(strsplit(as.character(data$SNP), split = ":"))[[1]]
data$BP <- ldply(strsplit(as.character(data$SNP), split = ":"))[[2]]
genesTemp <- read.table("listofgenes.txt", sep = "", header = F)
colnames(genesTemp) <- c("CHR","START","STOP","GENE")
genes <- subset(genesTemp, CHR != "X")
for(i in 1:length(genes$GENE))
{
  thisGene <- genes$GENE[i]
  thisChr <- genes$CHR[i]
  lower <- genes$START[i]
  upper <- genes$STOP[i]
  output <- subset(data, CHR == thisChr & BP >= lower & BP <= upper)
  fwrite(output, file = paste(thisGene,"_variants.tab", sep = ""), na = "NA", quote = F, row.names = F, sep = "\t")
}
```

## GWAS analysis 

RHOT1 positions on hg19

### Fisher exact test 

```
plink --bfile /data/LNG/saraB/HARDCALLS_PD_september_2018_no_cousins --remove-fam /data/LNG/saraB/NeuroX.fID.txt --chr 5 --geno 0.15 --from-bp 159990127 --to-bp 160279221 --make-bed --out RHOT1.GWAS

plink --bfile RHOT1.GWAS --fisher --covar /data/LNG/saraB/IPDGC_all_samples_covariates.tab --covar-name sex,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --out RHOT1_GWAS --make-bed --ci 0.95
```

### Annotate VCF with ANNOVAR 

```
plink --bfile RHOT1.GWAS --recode 'vcf-fid' --out RHOT1.GWAS
bgzip RHOT1.GWAS.vcf
tabix -f -p vcf RHOT1.GWAS.vcf.gz

table_annovar.pl /data/LNG/saraB/RHOT1/hardcallsNoNeuroX/vcf/RHOT1.GWAS.vcf.gz /data/LNG/saraB/annovar/humandb/ -buildver hg19 \
--thread 16 \
-out /data/LNG/saraB/RHOT1/hardcallsNoNeuroX/annotation/RHOT1.GWAS.annovar \
-remove -protocol avsnp147,refGene,ensGene,gnomad211_genome \
-operation f,g,g,f \
-nastring . \
-vcfinput

cd /data/LNG/saraB/RHOT1/hardcallsNoNeuroX/annotation
head -1 RHOT1.GWAS.annovar.hg19_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct RHOT1.GWAS.annovar.hg19_multianno.txt > RHOT1.GWAS.trimmed.annotation.txt
```

### Burden analysis

#### Generating list of coding variants

```
awk '$6=="exonic" {print}' RHOT1.GWAS.trimmed.annotation.txt > RHOT1.GWAS.trimmed.annotation.coding.variants.txt
awk '{print $1" "$2" "$2" "$7}' RHOT1.GWAS.trimmed.annotation.coding.variants.txt > RHOT1.trimmed.annotation.coding.variants.SNPs.txt

plink --bfile RHOT1.GWAS  --recode 'vcf-fid' --extract range RHOT1.trimmed.annotation.coding.variants.SNPs.txt --out RHOT1.CODING.GWAS
bgzip RHOT1.CODING.GWAS.vcf
tabix -f -p vcf RHOT1.CODING.GWAS.vcf.gz
```

#### MAF < 0.03

#### ALL VARIANTS 

```
cd /data/LNG/saraB/RHOT1

rvtest --noweb --inVcf /data/LNG/saraB/RHOT1/hardcallsNoNeuroX/vcf/RHOT1.GWAS.vcf.gz --pheno /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile /data/LNG/saraB/refFlat_hg19.txt --freqUpper 0.03 --out hardcallsNoNeuroX/burden/BURDEN.RHOT1.maf03

```
#### CODING VARIANTS 

```
cd /data/LNG/saraB/RHOT1

rvtest --noweb --inVcf /data/LNG/saraB/RHOT1/hardcallsNoNeuroX/vcf/RHOT1.CODING.GWAS.vcf.gz --pheno /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile /data/LNG/saraB/refFlat_hg19.txt --freqUpper 0.03 --out hardcallsNoNeuroX/burden/BURDEN.RHOT1.CODING.maf03
```

------------------------------------------------------------------------------------------------------------------------------

# RHOT2

## WGS analysis 

### Fisher exact test 

RHOT2 positions on hg38

module load plink samtools rvtests R annovar

```
plink --bfile RHOT2_WGS_AMP_PD --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --make-bed --out RHOT2_WGS_AMP_PD_pheno
plink --bfile RHOT2_WGS_AMP_PD_pheno --update-sex /data/LNG/saraB/AMP_PD/toupdatesex.txt --make-bed --out RHOT2_WGS_AMP_PD_pheno_sex
plink --bfile RHOT2_WGS_AMP_PD_pheno_sex --fisher --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX, AGE, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 --out RHOT2_WGS_AMP_PD --ci 0.95
```

### Annotation (ANNOVAR)

```
plink --bfile RHOT2_WGS_AMP_PD_pheno_sex --recode 'vcf-fid' --out RHOT2_WGS_AMP_PD_pheno_sex
bgzip RHOT2_WGS_AMP_PD_pheno_sex.vcf
tabix -f -p vcf RHOT2_WGS_AMP_PD_pheno_sex.vcf.gz

table_annovar.pl RHOT2_WGS_AMP_PD_pheno_sex.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
--thread 16 \
-out RHOT2_WGS.annovar \
-remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
-operation g,f,f,f \
-nastring . \
-vcfinput

head -1 RHOT2_WGS.annovar.hg38_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct RHOT2_WGS.annovar.hg38_multianno.txt > RHOT2.trimmed.annotation.txt
```

### Burden analysis 

#### Generating list of coding variants

```
awk '$6=="exonic" {print}' RHOT2.trimmed.annotation.txt > RHOT2.trimmed.annotation.coding.variants.WGS.txt
awk '{print $1" "$2" "$2" "$7}' RHOT2.trimmed.annotation.coding.variants.WGS.txt > RHOT2.trimmed.annotation.coding.variants.WGS.SNPs.txt

plink --bfile RHOT2_WGS_AMP_PD_pheno_sex --extract range RHOT2.trimmed.annotation.coding.variants.WGS.txt --recode 'vcf-fid' --out RHOT2_CODING_AMP

bgzip RHOT2_CODING_AMP.vcf
tabix -f -p vcf RHOT2_CODING_AMP.vcf.gz
```

#### MAF < 0.03

#### ALL VARIANTS 

```
rvtest --noweb --inVcf RHOT2_WGS_AMP_PD_pheno_sex.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 --out AMP_PD_BURDEN.RHOT2.maf003
```

#### CODING VARIANTS 

```
rvtest --noweb --inVcf RHOT2_CODING_AMP.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 --out AMP_PD_BURDEN.RHOT2.maf003_CODING
```
## GWAS analysis 

RHOT2 positions on hg19

### Fisher exact test 

```
plink --bfile /data/LNG/saraB/HARDCALLS_PD_september_2018_no_cousins --remove-fam /data/LNG/saraB/NeuroX.fID.txt --chr 5 --geno 0.15 --from-bp 159990127 --to-bp 160279221 --make-bed --out RHOT2.GWAS

plink --bfile RHOT2.GWAS --fisher --covar /data/LNG/saraB/IPDGC_all_samples_covariates.tab --covar-name sex,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --out RHOT2_GWAS --make-bed --ci 0.95
```

### Annotate VCF with ANNOVAR 

```
plink --bfile RHOT2.GWAS --recode 'vcf-fid' --out RHOT2.GWAS
bgzip RHOT2.GWAS.vcf
tabix -f -p vcf RHOT2.GWAS.vcf.gz

table_annovar.pl /data/LNG/saraB/RHOT2/hardcallsNoNeuroX/vcf/RHOT1.GWAS.vcf.gz /data/LNG/saraB/annovar/humandb/ -buildver hg19 \
--thread 16 \
-out /data/LNG/saraB/RHOT2/hardcallsNoNeuroX/annotation/RHOT2.GWAS.annovar \
-remove -protocol avsnp147,refGene,ensGene,gnomad211_genome \
-operation f,g,g,f \
-nastring . \
-vcfinput

cd /data/LNG/saraB/RHOT2/hardcallsNoNeuroX/annotation
head -1 RHOT2.GWAS.annovar.hg19_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct RHOT2.GWAS.annovar.hg19_multianno.txt > RHOT2.GWAS.trimmed.annotation.txt
```

### Burden analysis

#### Generating list of coding variants

```
awk '$6=="exonic" {print}' RHOT2.GWAS.trimmed.annotation.txt > RHOT2.GWAS.trimmed.annotation.coding.variants.txt
awk '{print $1" "$2" "$2" "$7}' RHOT2.GWAS.trimmed.annotation.coding.variants.txt > RHOT2.trimmed.annotation.coding.variants.SNPs.txt


plink --bfile RHOT2.GWAS  --recode 'vcf-fid' --extract range RHOT2.trimmed.annotation.coding.variants.SNPs.txt --out RHOT2.CODING.GWAS
bgzip RHOT2.CODING.GWAS.vcf
tabix -f -p vcf RHOT2.CODING.GWAS.vcf.gz
```

#### MAF < 0.03

#### ALL VARIANTS 

```
cd /data/LNG/saraB/RHOT2

rvtest --noweb --inVcf /data/LNG/saraB/RHOT2/hardcallsNoNeuroX/vcf/RHOT2.GWAS.vcf.gz --pheno /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile /data/LNG/saraB/refFlat_hg19.txt --freqUpper 0.03 --out hardcallsNoNeuroX/burden/BURDEN.RHOT2.maf03

```
#### CODING VARIANTS 

```
cd /data/LNG/saraB/RHOT2

rvtest --noweb --inVcf /data/LNG/saraB/RHOT2/hardcallsNoNeuroX/vcf/RHOT2.CODING.GWAS.vcf.gz --pheno /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile /data/LNG/saraB/refFlat_hg19.txt --freqUpper 0.03 --out hardcallsNoNeuroX/burden/BURDEN.RHOT2.CODING.maf03
```
