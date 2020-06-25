# *NUS1* data processing, variant extraction, selection and association with Parkinson's disease
Bernabe Bustos, Ph.D. (Northwestern University)


## **Burden of rare (freq < 0.01) nonsynonymous damaging variants**
### **IPDGC-exome data**
**Variant extraction, annotation and selection**
```
#Genomic coordinates in hg19

cat > NUS1_hg19.bed
chr6	117986617	118041886

#Extract variants from VCF file

module load bcftools

bcftools view -R NUS1_hg19.bed IPDGC_VarsQC_FullIndivQc.rawID_NoMaf_autosomes.NoMono.idiopathic.NoSuperDup.vcf.gz -o NUS1_var.vcf

#Annotation with annovar

/projects/b1049/genetics_programs/annovar_2017/annovar/table_annovar.pl NUS1_var.vcf /projects/b1049/genetics_programs/annovar_2017/annovar/humandb/ -buildver hg19 \
--thread 10 \
-out NUS1_var.annovar \
-remove -protocol refGene,dbnsfp30a,gnomad_exome,gnomad_genome,clinvar_20190305 \
-operation g,f,f,f,f \
-nastring . \
-vcfinput

#Selection of rare noncoding damaging variants

cut -f1-10,29,67,72 NUS1_var.annovar.hg19_multianno.txt | cat

Chr	Start	End	Ref	Alt	Func.refGene	Gene.refGene	GeneDetail.refGene	ExonicFunc.refGene	AAChange.refGene	CADD_phred	Freq ID
6	118014227	118014227	C	T	exonic	NUS1	.	synonymous SNV	NUS1:NM_138459:exon2:c.C438T:p.S146S	.	0.0003479	6_118014227_C_T
6	118014278	118014278	T	C	exonic	NUS1	.	synonymous SNV	NUS1:NM_138459:exon2:c.T489C:p.D163D	.	0.0006798	6_118014278_T_C
6	118014295	118014295	C	G	exonic	NUS1	.	nonsynonymous SNV	NUS1:NM_138459:exon2:c.C506G:p.P169R	19.26	0.000342	6_118014295_C_G
6	118015292	118015292	A	G	exonic	NUS1	.	nonsynonymous SNV	NUS1:NM_138459:exon3:c.A640G:p.K214E	18.02	0.0003376	6_118015292_A_G
6	118024775	118024775	T	C	exonic	NUS1	.	synonymous SNV	NUS1:NM_138459:exon4:c.T699C:p.N233N	.	0.0003353	6_118024775_T_C
6	118024832	118024832	C	T	exonic	NUS1	.	synonymous SNV	NUS1:NM_138459:exon4:c.C756T:p.G252G	.	0.0003349	6_118024832_C_T


cut -f1-10,29,67,72 NUS1_var.annovar.hg19_multianno.txt | \
awk '$6 == "exonic" && $9 == "nonsynonymous" && $12 > 12.37 && $13 < 0.01{print $NF}' > NUS1_funcVars.txt
cat NUS1_funcVars.txt

6_118014295_C_G
6_118015292_A_G

```
**Allelic counts and contingency table**
```
module load plink/1.9-beta-4.6

#vcf to plink

plink --vcf IPDGC_VarsQC_FullIndivQc.rawID_NoMaf_autosomes.NoMono.idiopathic.NoSuperDup.vcf.gz \
--double-id \
--biallelic-only strict list \
--vcf-half-call missing \
--make-bed \
--out IPDGC_VarsQC_FullIndivQc.rawID_NoMaf_autosomes.NoMono.idiopathic.NoSuperDup

#allelic counts in cases

plink --bfile IPDGC_VarsQC_FullIndivQc.rawID_NoMaf_autosomes.NoMono.idiopathic.NoSuperDup \
--pheno IPDGC_exomes_pheno.txt \
--extract NUS1_funcVars.txt \
--filter-cases \
--counts \
--out NUS1_funcVars.cases_counts.txt

#allelic counts in controls

plink --bfile IPDGC_VarsQC_FullIndivQc.rawID_NoMaf_autosomes.NoMono.idiopathic.NoSuperDup \
--pheno IPDGC_exomes_pheno.txt \
--extract NUS1_funcVars.txt \
--filter-controls \
--counts \
--out NUS1_funcVars.controls_counts.txt

#contingency table


           | Present  | Absent   |
| -------- | -------- | -------- |
| Cases    |    2     |   2082   |
| -------- | -------- | -------- |
| Controls |    0     |    904   |

#Fisher's exact test in -> https://www.graphpad.com/quickcalcs/contingency1/
P=1.00
```
### **AMP-PD WGS data**
**Variant extraction, annotation and selection**

```
#Genomic coordinates in hg38

cat > NUS1_hg38.bed
chr6    117665469    117720727

#Extract variants from Plink file

module load plink/1.9-beta-4.6 samtools

plink --bfile releases_2019_v1release_1015_wgs_plink_bfile_all_vcfs.IndivQC \
--extract range NUS1_hg38.bed \
--make-bed \
--out NUS1_WGS_AMP_PD_pheno_sex

#Annotation with annovar

plink --bfile NUS1_WGS_AMP_PD_pheno_sex --recode 'vcf-fid' --out NUS1_WGS_AMP_PD_pheno_sex
bgzip NUS1_WGS_AMP_PD_pheno_sex.vcf
tabix -f -p vcf NUS1_WGS_AMP_PD_pheno_sex.vcf.gz

table_annovar.pl NUS1_WGS_AMP_PD_pheno_sex.vcf.gz \ 
humandb/ -buildver hg38 \
--thread 10 \
-out NUS1_WGS_AMP_PD_pheno_sex.annovar \
-remove -protocol refGene,dbnsfp35a,gnomad30_genome,clinvar_20200316 \
-operation g,f,f,f \
-nastring . \
-vcfinput


cut -f1-10,29,72 NUS1_WGS_AMP_PD_pheno_sex.annovar.hg38_multianno.txt | \
head

Chr	Start	End	Ref	Alt	Func.refGene	Gene.refGene	GeneDetail.refGene	ExonicFunc.refGene	AAChange.refGene	CADD_phred	Otherinfo ID
6	117665484	117665484	C	G	intergenic	GOPC;LOC101927919	dist=62973;dist=7329	.	.	.	0.04765	rs117390884
6	117665536	117665536	G	C	intergenic	GOPC;LOC101927919	dist=63025;dist=7277	.	.	.	0.01224	rs180774613
6	117665553	117665553	G	C	intergenic	GOPC;LOC101927919	dist=63042;dist=7260	.	.	.	0.2718	rs9374675
6	117665672	117665672	C	T	intergenic	GOPC;LOC101927919	dist=63161;dist=7141	.	.	.	0	rs184306021
6	117666119	117666120	AT	-	intergenic	GOPC;LOC101927919	dist=63608;dist=6693	.	.	.	0.458	rs71723796
6	117666237	117666237	T	C	intergenic	GOPC;LOC101927919	dist=63726;dist=6576	.	.	.	0	rs189259197
6	117666334	117666334	T	A	intergenic	GOPC;LOC101927919	dist=63823;dist=6479	.	.	.	0.0005562	rs114693732
6	117666345	117666345	A	G	intergenic	GOPC;LOC101927919	dist=63834;dist=6468	.	.	.	0.04153	rs77576981
6	117666671	117666671	C	T	intergenic	GOPC;LOC101927919	dist=64160;dist=6142	.	.	.	0.1409	rs11153691


cut -f1-10,52,99,104 NUS1_WGS_AMP_PD_pheno_sex.annovar.hg38_multianno.txt | \
awk '$6 == "exonic" && $9 == "nonsynonymous" && $12 > 12.37 && $13 < 0.01{print $NF}' > NUS1_AMP_PD_funcVars.txt

cat NUS1_AMP_PD_funcVars.txt

rs150646335
rs146171115

```
**Allelic counts and contingency table**
```
module load plink/1.9-beta-4.6

#allelic counts in cases

plink --bfile NUS1_WGS_AMP_PD_pheno_sex \
--extract NUS1_AMP_PD_funcVars.txt \
--filter-cases \
--counts \
--out NUS1_AMP_PD_funcVars.cases_counts.txt

#allelic counts in controls

plink --bfile NUS1_WGS_AMP_PD_pheno_sex \
--extract NUS1_AMP_PD_funcVars.txt \
--filter-controls \
--counts \
--out NUS1_AMP_PD_funcVars.controls_counts.txt

#contingency table


           | Present  | Absent   |
| -------- | -------- | -------- |
| Cases    |    6     |   3288   |
| -------- | -------- | -------- |
| Controls |    6     |   2094   |

#Fisher's exact test in -> https://www.graphpad.com/quickcalcs/contingency1/
P=0.41
```
**Combined burden test**
```
           | Present  | Absent   |
| -------- | -------- | -------- |
| Cases    |    8     |   5370   |
| -------- | -------- | -------- |
| Controls |    6     |   2998   |

#Fisher's exact test in -> https://www.graphpad.com/quickcalcs/contingency1/
P=0.58
```

## **Burden of rare (freq < 0.01) coding variants**
### **IPDGC-exome data**

**Selection and extraction of variants**
```
cut -f1-10,29,67,72 NUS1_var.annovar.hg19_multianno.txt | \
awk '$6 == "exonic" && $13 < 0.01{print $NF}' > NUS1_SKATO_Vars.txt

module load vcftools/1.0.7 samtools 

vcftools --gzvcf NUS1_hg19.bed IPDGC_VarsQC_FullIndivQc.rawID_NoMaf_autosomes.NoMono.idiopathic.NoSuperDup.vcf.gz \
--snps NUS1_SKATO_Vars.txt \
--recode \
--recode-INFO-all \
--stdout | bgzip > IPDGC_exomes_NUS1_SKATO_Vars.vcf.gz

tabix -p vcf IPDGC_exomes_NUS1_SKATO_Vars.vcf.gz
```

**Run RVTests**
```
rvtest --inVcf IPDGC_exomes_NUS1_SKATO_Vars.vcf.gz \
--pheno ../PDGC_QC_CR90.rawID_wSEX_wPHENO_MIND05_NoHighHe_NoIBD_NoPCA_SexCheck_StevenQC_All_Vars.fam \
--covar IPDGC_allQC.covar \
--covar-name SEX,PC1,PC2,PC3,PC4,PC5 \
--out NUS1_coding_allVars_wCovs \
--geneFile refFlat_hg19.txt \
--kernel skato
```
### **AMP-PD WGS data**

**Selection and extraction of variants**
```
cut -f1-10,52,99,104 NUS1_WGS_AMP_PD_pheno_sex.annovar.hg38_multianno.txt | \
awk '$6 == "exonic" && $13 < 0.01{print $NF}' > NUS1_AMP_PD_SKATO_vars.txt

module load vcftools/1.0.7 samtools 

vcftools --gzvcf NUS1_WGS_AMP_PD_pheno_sex.vcf.gz \
--snps NUS1_AMP_PD_SKATO_vars.txt \
--recode \
--recode-INFO-all \
--stdout | bgzip > AMP_PD_WGS_NUS1_SKATO_Vars.vcf.gz

tabix -p vcf AMP_PD_WGS_NUS1_SKATO_Vars.vcf.gz
```

**Run RVTests**
```
rvtest --inVcf AMP_PD_WGS_NUS1_SKATO_Vars.vcf.gz \
--pheno covs_Mike.txt \
--covar covs_Mike.txt \
--covar-name SEX AAO PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 \
--out AMP_PD_NUS1_coding_allVars_wCovs \
--geneFile refFlat_hg38.txt \
--kernel skato
```
### **IPDGC-GWAS cohorts**

```
## No coding variants present ##
```

## **Single variant tests**

### **IPDGC-exome data**
**Variants with MAF < 1%**
```
module load plink/1.9-beta-4.6

# Extraction of all vars in NUS1

plink --bfile IPDGC_VarsQC_FullIndivQc.rawID_NoMaf_autosomes.NoMono.idiopathic.NoSuperDup --extract range NUS1_hg19.bed --make-bed --out NUS1_IPDGC_exomes_allVars

# Model association

plink --bfile NUS1_IPDGC_exomes_allVars \
--model fisher \
--max-maf 0.01 \
--out NUS1_IPDGC_exomes_RareVars
```
**Variants with MAF > 1%**
```
## No variants left for analysis ##
```

### **AMP-PD WGS data**
**Variants with MAF < 1%**
```
plink --bfile NUS1_WGS_AMP_PD_pheno_sex \
--model fisher \
--max-maf 0.01 \
--out NUS1_WGS_AMP_PD_pheno_sex_RareVars
```
**Variants with MAF > 1%**
```
plink --bfile NUS1_WGS_AMP_PD_pheno_sex \
--logistic hide-covar \
--maf 0.01 \
--covar covs_Mike.txt \
--covar-name SEX, AAO, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 \
--ci 0.95 \
--out NUS1_WGS_AMP_PD_pheno_sex_CommonVars

```
### **IPDGC-GWAS cohorts**
**Variants with MAF < 1%**
```
# Extraction of all vars in NUS1

plink --bfile HARDCALLS_PD_september_2018_no_cousins \
--remove-fam NeuroX.fID.txt \
--chr 6 \
--geno 0.15 \
--from-bp 117986617 \
--to-bp 118041886 \
--make-bed \
--out NUS1.GWAS

# Model association

plink --bfile NUS1.GWAS \
--model fisher \
--max-maf 0.01 \
--out NUS1_GWAS_model_RareVars
```
**Variants with MAF > 1%**
```
plink --bfile NUS1.GWAS \
--logistic hide-covar \
--covar IPDGC_all_samples_covariates.tab \
--covar-name sex AGE PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 Dataset \
--ci 0.95 \
--out NUS1_GWAS_CommonVars
```

## ***NUS1* summary stats from PD GWAS meta-analysis without 23andMe data**
```
module load R/3.6.3

cat listofgenes.txt
chr6	117986617	118041886    NUS1

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