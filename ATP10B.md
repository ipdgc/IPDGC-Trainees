# ATP10B

Raquel Real (UCL), Anni Moore (NIH), Sara Bandres-Ciga (NIH)

## WGS analysis 

### Fisher exact test 

ATP10B positions on hg38

module load plink samtools rvtests R annovar

```
plink --bfile ATP10B_WGS_AMP_PD --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --make-bed --out ATP10B_WGS_AMP_PD_pheno
plink --bfile ATP10B_WGS_AMP_PD_pheno --update-sex /data/LNG/saraB/AMP_PD/toupdatesex.txt --make-bed --out ATP10B_WGS_AMP_PD_pheno_sex
plink --bfile ATP10B_WGS_AMP_PD_pheno_sex --fisher --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX, AGE, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 --out ATP10B_WGS_AMP_PD --ci 0.95
```

### Annotation (ANNOVAR)

```
plink --bfile ATP10B_WGS_AMP_PD_pheno_sex --recode 'vcf-fid' --out ATP10B_WGS_AMP_PD_pheno_sex
bgzip ATP10B_WGS_AMP_PD_pheno_sex.vcf
tabix -f -p vcf ATP10B_WGS_AMP_PD_pheno_sex.vcf.gz

table_annovar.pl /data/LNG/saraB/AMP_PD/genesfortrainees/ATP10B/ATP10B_WGS_AMP_PD_pheno_sex.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
--thread 16 \
-out /data/LNG/saraB/AMP_PD/genesfortrainees/ATP10B/ATP10B_WGS.annovar \
-remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
-operation g,f,f,f \
-nastring . \
-vcfinput

head -1 ATP10B_WGS.annovar.hg38_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct ATP10B_WGS.annovar.hg38_multianno.txt > ATP10B.trimmed.annotation.txt
```

### Burden analysis 

#### Generating list of coding variants

```
awk '$6=="exonic" {print}' ATP10B.trimmed.annotation.txt > ATP10B.trimmed.annotation.coding.variants.WGS.txt
awk '{ print $1":"$2 }' ATP10B.trimmed.annotation.coding.variants.WGS.txt > ATP10B.trimmed.annotation.coding.variants..WGS.SNPs.txt

plink --bfile ATP10B_WGS_AMP_PD_pheno_sex --extract range ATP10B.trimmed.annotation.coding.variants.WGS.SNPs.txt  ---recode 'vcf-fid' --out ATP10B_CODING_AMP

bgzip ATP10B_CODING_AMP.vcf
tabix -f -p vcf ATP10B_CODING_AMP.vcf.gz
```

#### MAF < 0.03

#### ALL VARIANTS 

```
rvtest --noweb --inVcf ATP10B_WGS_AMP_PD_pheno_sex.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 --out AMP_PD_BURDEN.ATP10B.maf003
```

#### CODING VARIANTS 

```
rvtest --noweb --inVcf ATP10B_CODING_AMP.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 --out AMP_PD_BURDEN.ATP10B.maf003_CODING
```

### Extract your gene(s) from PD + AAO GWAS summary stats 

* listofgenes.txt will have the following format: chr start stop gene

```
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

data2 <- fread("IPDGC_AAO_GWAS_sumstats_april_2018.txt", header = T)
data2$CHR <- ldply(strsplit(as.character(data2$MarkerName), split = ":"))[[1]]
data2$BP <- ldply(strsplit(as.character(data2$MarkerName), split = ":"))[[2]]
#colnames(data2) <- c("chr","pos","no","bp","A1","A2")
genesTemp <- read.table("listofgenes.txt", sep = "", header = F)
colnames(genesTemp) <- c("CHR","START","STOP","GENE")
genes <- subset(genesTemp, CHR != "X")
for(i in 1:length(genes$GENE))
{
  thisGene <- genes$GENE[i]
  thisChr <- genes$CHR[i]
  lower <- genes$START[i]
  upper <- genes$STOP[i]
  output <- subset(data2, CHR == thisChr & BP >= lower & BP <= upper)
  fwrite(output, file = paste(thisGene,"_variants_AAO.tab", sep = ""), na = "NA", quote = F, row.names = F, sep = "\t")
}
```

## GWAS analysis 

ATP10B positions on hg19

### Fisher exact test 

```
plink --bfile /data/LNG/saraB/HARDCALLS_PD_september_2018_no_cousins --remove-fam /data/LNG/saraB/NeuroX.fID.txt --chr 5 --geno 0.15 --from-bp 159990127 --to-bp 160279221 --make-bed --out ATP10B.GWAS

plink --bfile ATP10B.GWAS --fisher --covar /data/LNG/saraB/IPDGC_all_samples_covariates.tab --covar-name sex,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --out ATP10B_GWAS --make-bed --ci 0.95
```

### Annotate VCF with ANNOVAR 

```
plink --bfile ATP10B.GWAS --recode 'vcf-fid' --out ATP10B.GWAS
bgzip ATP10B.GWAS.vcf
tabix -f -p vcf ATP10B.GWAS.vcf.gz

table_annovar.pl /data/LNG/saraB/ATP10B/hardcallsNoNeuroX/vcf/ATP10B.GWAS.vcf.gz /data/LNG/saraB/annovar/humandb/ -buildver hg19 \
--thread 16 \
-out /data/LNG/saraB/ATP10B/hardcallsNoNeuroX/annotation/ATP10B.GWAS.annovar \
-remove -protocol avsnp147,refGene,ensGene,gnomad211_genome \
-operation f,g,g,f \
-nastring . \
-vcfinput

cd /data/LNG/saraB/ATP10B/hardcallsNoNeuroX/annotation
head -1 ATP10B.GWAS.annovar.hg19_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct ATP10B.GWAS.annovar.hg19_multianno.txt > ATP10B.GWAS.trimmed.annotation.txt
```

### Burden analysis

#### Generating list of coding variants

```
awk '$6=="exonic" {print}' ATP10B.GWAS.trimmed.annotation.txt > ATP10B.GWAS.trimmed.annotation.coding.variants.txt
awk '{ print $1":"$2 }' ATP10B.GWAS.trimmed.annotation.coding.variants.txt > ATP10B.trimmed.annotation.coding.variants.SNPs.txt

plink --bfile ATP10B.GWAS  --recode 'vcf-fid' --extract range ATP10B.trimmed.annotation.coding.variants.SNPs.txt --out ATP10B.CODING.GWAS
bgzip ATP10B.CODING.GWAS.vcf
tabix -f -p vcf ATP10B.CODING.GWAS.vcf.gz
```

#### MAF < 0.03

#### ALL VARIANTS 

cd /data/LNG/saraB/ATP10B

```
rvtest --noweb --inVcf /data/LNG/saraB/ATP10B/hardcallsNoNeuroX/vcf/ATP10B.GWAS.vcf.gz --pheno /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile /data/LNG/saraB/refFlat_hg19.txt --freqUpper 0.03 --out hardcallsNoNeuroX/burden/BURDEN.ATP10B.maf03
```

#### CODING VARIANTS 

cd /data/LNG/saraB/ATP10B

```
rvtest --noweb --inVcf /data/LNG/saraB/ATP10B/hardcallsNoNeuroX/vcf/ATP10B.CODING.GWAS.vcf.gz --pheno /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar /data/LNG/saraB/IPDGC_all_samples_covariates.vcf.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile /data/LNG/saraB/refFlat_hg19.txt --freqUpper 0.03 --out hardcallsNoNeuroX/burden/BURDEN.ATP10B.CODING.maf03
```

### Look for enrichment of compound heterozygous in cases versus controls in WGS data

mkdir /data/LNG/anni/ATP10B_2/

#### Look at the frequencies of hom/compound hets in the WGS

(extract coding variants in homozygous state in cases and controls)

#### Generating list of coding variants

```
awk '$6=="exonic" {print}' ATP10B.trimmed.annotation.txt > ATP10B.trimmed.annotation.coding.variants.WGS.txt
awk '{ print $1":"$2 }' ATP10B.trimmed.annotation.coding.variants.WGS.txt > ATP10B.trimmed.annotation.coding.variants..WGS.SNPs.txt
```

#### Generating list of nonsynonymous variants

```
awk '$9=="nonsynonymous" {print}' ATP10B.trimmed.annotation.coding.variants.WGS.txt > ATP10B.trimmed.annotation.coding.variants.WGS.nonsynonymous.txt
awk '{ print $1":"$2 }' ATP10B.trimmed.annotation.coding.variants.WGS.nonsynonymous.txt > ATP10B.trimmed.annotation.coding.variants.WGS.nonsynonymous.SNPs.txt
```

#### Reformat 2nd column of plink file from rs ids to chr:pos to be able to extract

```
awk '{print $1, $1":"$4",$3,$4,$5}' ATP10B_WGS_AMP_PD_pheno_sex.bim > ATP10B_WGS_AMP_PD_pheno_sex_mod.bim
```

#### Change names to match new bim file

```
cp ATP10B_WGS_AMP_PD_pheno_sex.fam ATP10B_WGS_AMP_PD_pheno_sex_mod.fam
cp ATP10B_WGS_AMP_PD_pheno_sex.bed ATP10B_WGS_AMP_PD_pheno_sex_mod.bed
```

#### Filter and extract coding variants

```
plink --bfile /data/LNG/anni/ATP10B_2/ATP10B_WGS_AMP_PD_pheno_sex_mod --chr 5 --geno 0.15 --extract /data/LNG/anni/ATP10B_2/ATP10B.trimmed.annotation.coding.variants.SNPs.WGS.txt --make-bed --out /data/LNG/anni/ATP10B_2/ATP10B.WGS.CODING.VARIANTS.all
```

#### Recode to get genotypes of coding variants

```
plink --bfile /data/LNG/anni/ATP10B_2/ATP10B.WGS.CODING.VARIANTS.all --recode A --out /data/LNG/anni/ATP10B_2/ATP10B.WGS.CODING.VARIANTS.all.recode
plink --bfile /data/LNG/anni/ATP10B_2/ATP10B.WGS.CODING.VARIANTS.maxmaf0.03 --recode A --out /data/LNG/anni/ATP10B_2/ATP10B.WGS.CODING.VARIANTS.maxmaf0.03.recode
```

#### Get genotypes as 2 (homoz for alternate allele), 1 (heterozygous) and 0 homoz for reference allele

```
*extract individuals (cases and controls) with genotypes =2 and individuals with genotypes =1 in more than 1 variant
*done in python locally
*/Users/mooreank/Desktop/Sara/ATP10B_2/ATP10B_genotype_check.ipynb
```
