# PTRHD1
Yuri Sosero, Lynne Krohn (McGill)

## GWAS analysis - IPDGC cohort (Nalls et al., 2019).

### 1. Subset PLINK Binaries + Convert to VCFs

PTRHD1 positions on hg19

module load plink samtools rvtests R annovar

```{bash}
plink --bfile /data/LNG/saraB/HARDCALLS_PD_september_2018_no_cousins --remove-fam /data/LNG/saraB/NeuroX.fID.txt --geno 0.05 --chr 2 --from-bp 24913136 --to-bp 25116251 --make-bed --out hardcallsNoNeuroX/bin/PTRHD1.GWAS

plink --bfile hardcallsNoNeuroX/bin/PTRHD1.GWAS --recode vcf --out hardcallsNoNeuroX/vcf/PTRHD1.GWAS

cd hardcallsNoNeuroX/vcf 
bgzip PTRHD1.GWAS.vcf 
tabix -f -p vcf PTRHD1.GWAS.vcf.gz
```

### 2. Fisher exact test and logistic regression 

```{bash}
cd /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/PTRHD1

plink --bfile hardcallsNoNeuroX/bin/PTRHD1.GWAS --chr 2 --from-bp 24913136 --to-bp 25116251 --fisher --out hardcallsNoNeuroX/freq/PTRHD1

plink --bfile hardcallsNoNeuroX/bin/PTRHD1.GWAS --chr 2 --from-bp 24913136 --to-bp 25116251 --covar /data/LNG/saraB/IPDGC_all_samples_covariates.tab --covar-name sex,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --logistic --out hardcallsNoNeuroX/logistic/PTRHD1
```


### 3. Annotation 

```{bash}
mkdir /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/PTRHD1/annovar/

# Downloading databases 
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
# annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/ 
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 humandb/ 
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a humandb/
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_genome humandb/

module load annovar

cd /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/PTRHD1/hardcallsNoNeuroX/vcf

table_annovar.pl PTRHD1.GWAS.vcf.gz /data/LNG/saraB/annovar/humandb/ -buildver hg19 \
--thread 16 \
-out PTRHD1.annovar \
-remove -protocol avsnp147,refGene,ensGene,gnomad211_genome \
-operation f,g,g,f \
-nastring . \
-vcfinput

head -1 PTRHD1.annovar.hg19_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct PTRHD1.annovar.hg19_multianno.txt > PTRHD1.GWAS.trimmed.annotation.txt
```

### 4. Selection of homozygous and compound heterozygous variants 

```{python}
import pandas as pd

import os
from google.colab import drive
drive.mount('/content/drive/')
os.chdir("/content/drive/My Drive/")

## RECAP of code ran on plink

## Generating list of coding variants without synonymous variants
##awk '$10=="nonsynonymous" {print}' PTRHD1.GWAS.trimmed.annotation.coding.variants.txt > PPTRHD1.GWAS.trimmed.annotation.coding.variants_non_syn.txt
##awk '{print $1" "$2" "$2" "$7}' PPTRHD1.GWAS.trimmed.annotation.coding.variants_non_syn.txt > PPTRHD1.GWAS.trimmed.annotation.coding.variants_non_syn.SNPs.txt

## Filter and extract coding variants without synonymous variants
##plink --bfile PTRHD1.GWAS --chr 2 --extract range PPTRHD1.GWAS.trimmed.annotation.coding.variants_non_syn.SNPs.txt --make-bed --out PTRHD1.GWAS.CODING.not_syn

## Recode to get genotypes of coding variants
##plink --bfile PTRHD1.GWAS.CODING.not_syn --recode A --out PTRHD1.GWAS.CODING.not_syn

"""### Load in genotype file and get genotypes as 2 (homoz for alternate allele), 1 (heterozygous) and 0 homoz for reference allele"""

df = pd.read_csv('PTRHD1.GWAS.CODING.not_syn.raw', sep = ' ')
df.head()

df.shape

"""#### making subset table of sample info"""

info = df[['FID','IID','PAT','MAT','SEX','PHENOTYPE']]

df.columns

"""#### making subset table without info"""

df_sub = df[['FID', 'IID', '2:24974958_T',
       '2:24980955_G', '2:25022598_G']]

df_sub.head()

df_sub.shape

"""### Identifying patients homozygous for alt allele"""

homo = df_sub[df_sub.values  == 2]
homo = homo.drop_duplicates()
print(homo.shape)
homo.head()

homo.columns

homo_info = pd.merge(homo, info,  how='left', left_on=['FID','IID'], right_on = ['FID','IID'])
homo_info = homo_info[['FID','IID','PAT','MAT','SEX','PHENOTYPE','2:24974958_T',
       '2:24980955_G', '2:25022598_G']]


homo_info = homo_info.dropna()
homo_info.head()

homo_info.shape

"""##### removing sample info"""

del homo_info['FID']
del homo_info['IID']
homo_info.head()

homo_info_control = homo_info.loc[homo_info['PHENOTYPE'] == 1]
print(homo_info_control.shape)
homo_info_case = homo_info.loc[homo_info['PHENOTYPE'] == 2]
print(homo_info_case.shape)

homo_info.to_csv('PTRHD1.GWAS.CODING.VARIANTS.recode.homo.txt', sep = '\t', index=False)

homo_info_control.to_csv('CONTROLS_PTRHD1.GWAS.CODING.VARIANTS.recode.homo.txt', sep = '\t', index=False)

homo_info_case.to_csv('CASES_PTRHD1.GWAS.CODING.VARIANTS.recode.homo.txt', sep = '\t', index=False)

"""### Identifying patients heterozygous for an alt allele"""

het = df_sub[df_sub.values == 1]
het = het.drop_duplicates()
print(het.shape)
het.head()

het.columns

het['count'] = het[['2:24974958_T',
       '2:24980955_G', '2:25022598_G']].sum(axis=1)
het.head()

"""#### finding samples with more than one alt allele"""

het_comp = het.loc[het['count'] > 1]
het_comp = het_comp.drop_duplicates()
print(het_comp.shape)
het_comp.head()

new_df = pd.merge(het_comp, info,  how='left', left_on=['FID','IID'], right_on = ['FID','IID'])
new_df = new_df[['FID','IID','PAT','MAT','SEX','PHENOTYPE','2:24974958_T',
       '2:24980955_G', '2:25022598_G', 'count']]

new_df = new_df.dropna()
new_df.head()

print(new_df.shape)

new_df.head()

"""#### removing sample info"""

del new_df['FID']
del new_df['IID']

controls = new_df.loc[new_df['PHENOTYPE'] == 1]
cases = new_df.loc[new_df['PHENOTYPE'] == 2]

print(controls.shape)
print(cases.shape)

new_df.to_csv('PTRHD1.GWAS.CODING.VARIANTS.recode.compound.hets.txt', sep = '\t', index=True)
new_df.head()

cases.to_csv('CASES_PTRHD1.GWAS.CODING.VARIANTS.recode.compound.hets.txt', sep = '\t', index=True)

controls.to_csv('CONTROLS_PTRHD1.GWAS.CODING.VARIANTS.recode.compound.hets.txt', sep = '\t', index=True)
```


## WGS analysis - AMP-PD cohort (www.amp-pd.org)


### 1. Fisher exact test and logistic regression

PTRHD1 positions on hg38

module load plink samtools rvtests R annovar

```{bash}
plink --bfile PTRHD1_WGS_AMP_PD --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --make-bed --out PTRHD1_WGS_AMP_PD_pheno
plink --bfile PTRHD1_WGS_AMP_PD_pheno --update-sex /data/LNG/saraB/AMP_PD/toupdatesex.txt --make-bed --out PTRHD1_WGS_AMP_PD_pheno_sex

plink --bfile PTRHD1_WGS_AMP_PD_pheno_sex --fisher --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --out PTRHD1_WGS_AMP_PD
plink --bfile PTRHD1_WGS_AMP_PD_pheno_sex --logistic --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX, AGE, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 --out PTRHD1_WGS_AMP_PD --ci 0.95
```



### 2. Annotation

```{bash}
plink --bfile PTRHD1_WGS_AMP_PD_pheno_sex --recode 'vcf-fid' --out PTRHD1_WGS_AMP_PD_pheno_sex
bgzip PTRHD1_WGS_AMP_PD_pheno_sex.vcf
tabix -f -p vcf PTRHD1_WGS_AMP_PD_pheno_sex.vcf.gz

table_annovar.pl PTRHD1_WGS_AMP_PD_pheno_sex.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
--thread 16 \
-out PTRHD1_WGS.annovar \
-remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
-operation g,f,f,f \
-nastring . \
-vcfinput

head -1 PTRHD1_WGS.annovar.hg38_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct PTRHD1_WGS.annovar.hg38_multianno.txt > PTRHD1.WGS.trimmed.annotation.txt
```


### 3. Burden analysis

```{bash}
#### Generating list of coding variants

awk '$6=="exonic" {print}' PTRHD1.WGS.trimmed.annotation.txt > PTRHD1.trimmed.annotation.coding.variants.WGS.txt
awk '{print $1" "$2" "$2" "$7}' PTRHD1.trimmed.annotation.coding.variants.WGS.txt > PTRHD1.trimmed.annotation.coding.variants.WGS.SNPs.txt

plink --bfile PTRHD1_WGS_AMP_PD_pheno_sex --extract range PTRHD1.trimmed.annotation.coding.variants.WGS.txt --recode 'vcf-fid' --out PTRHD1_CODING_AMP

bgzip PTRHD1_CODING_AMP.vcf
tabix -f -p vcf PTRHD1_CODING_AMP.vcf.gz

#### MAF < 0.01

#### ALL VARIANTS 

rvtest --noweb --inVcf PTRHD1_WGS_AMP_PD_pheno_sex.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.01 --out AMP_PD_BURDEN.PTRHD1.WGS.maf001

#### CODING VARIANTS 

rvtest --noweb --inVcf PTRHD1_CODING_AMP.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.01 --out AMP_PD_BURDEN.PTRHD1.WGS.maf001_CODING
```



### 4. Selection of homozygous and compound heterozygous variants

```{python}
import pandas as pd

## RECAP of code run on plink

## Generating list of coding variants
##awk '$6=="exonic" {print}' PTRHD1.WGS.trimmed.annotation.txt > PTRHD1.WGS.trimmed.annotation.coding.variants.WGS.txt
##awk '{print $1" "$2" "$2" "$7}' PTRHD1.WGS.trimmed.annotation.coding.variants.WGS.txt > PTRHD1.WGS.trimmed.annotation.coding.variants.WGS.SNPs.txt

## Generating list of coding variants without synonymous variants
##awk '$9=="nonsynonymous" {print}' PTRHD1.WGS.trimmed.annotation.txt > PTRHD1.WGS.trimmed.annotation.coding.variants.WGS_nonsnynonymous.txt
##awk '{print $1" "$2" "$2" "$7}' PTRHD1.WGS.trimmed.annotation.coding.variants.WGS_nonsnynonymous.txt > PTRHD1.WGS.trimmed.annotation.coding.variants.WGS_nonsnynonymous.SNPs.txt

## Filter and extract coding variants without synonymous variants
##plink --bfile PTRHD1_WGS_AMP_PD_pheno_sex --chr 2 --extract range PTRHD1.WGS.trimmed.annotation.coding.variants.WGS_nonsnynonymous.SNPs.txt --make-bed --out PTRHD1.WGS.CODING.not_syn

## Recode to get genotypes of coding variants
##plink --bfile PTRHD1.WGS.CODING.not_syn --recode A --out PTRHD1.WGS.CODING.not_syn

"""### Load in genotype file and get genotypes as 2 (homoz for alternate allele), 1 (heterozygous) and 0 homoz for reference allele"""

df = pd.read_csv('PTRHD1.WGS.CODING.not_syn.raw', sep = ' ')
df.head()

df.shape

"""#### making subset table of sample info"""

info = df[['FID','IID','PAT','MAT','SEX','PHENOTYPE']]

df.columns

"""#### making subset table without info"""

df_sub = df[['FID', 'IID', 'rs149553409_C',
       'rs202008308_G', 'rs148155916_C', 'rs142018995_T', 'rs149214507_C',
       'rs146794228_A', 'rs201790363_T', 'rs145705009_A', 'rs138104918_C',
       'rs149490066_G', 'rs370093383_G', 'rs56099330_A', 'rs1804645_T',
       'rs150066931_G', 'rs148719445_C', 'rs201252444_C', 'rs138426814_C',
       'rs201704717_C', 'rs200722485_A', 'rs143782986_T', 'rs1550116_G',
       'rs201920035_A', 'rs76741040_G', 'rs147672058_G', 'rs149340371_T',
       'rs148113144_C', 'rs146239018_T', 'rs139463930_T', 'rs143399833_A',
       'rs146553503_T', 'rs147179390_T', 'rs199524620_A', 'rs146165057_T',
       'rs185180064_T', 'rs190125602_A', 'rs78678013_T', 'rs115329263_A',
       'rs200642280_A', 'rs199731805_T', 'rs139407103_A', 'rs114530292_C',
       'rs372743845_T', 'rs201606553_T', 'rs149726253_A']]

df_sub.head()

df_sub.shape

"""### Identifying patients homozygous for alt allele"""

homo = df_sub[df_sub.values  == 2]
homo = homo.drop_duplicates()
print(homo.shape)
homo.head()

homo.columns

homo_info = pd.merge(homo, info,  how='left', left_on=['FID','IID'], right_on = ['FID','IID'])
homo_info = homo_info[['FID','IID','PAT','MAT','SEX','PHENOTYPE','rs149553409_C', 'rs202008308_G', 'rs148155916_C',
       'rs142018995_T', 'rs149214507_C', 'rs146794228_A', 'rs201790363_T',
       'rs145705009_A', 'rs138104918_C', 'rs149490066_G', 'rs370093383_G',
       'rs56099330_A', 'rs1804645_T', 'rs150066931_G', 'rs148719445_C',
       'rs201252444_C', 'rs138426814_C', 'rs201704717_C', 'rs200722485_A',
       'rs143782986_T', 'rs1550116_G', 'rs201920035_A', 'rs76741040_G',
       'rs147672058_G', 'rs149340371_T', 'rs148113144_C', 'rs146239018_T',
       'rs139463930_T', 'rs143399833_A', 'rs146553503_T', 'rs147179390_T',
       'rs199524620_A', 'rs146165057_T', 'rs185180064_T', 'rs190125602_A',
       'rs78678013_T', 'rs115329263_A', 'rs200642280_A', 'rs199731805_T',
       'rs139407103_A', 'rs114530292_C', 'rs372743845_T', 'rs201606553_T',
       'rs149726253_A']]


homo_info = homo_info.dropna()
homo_info.head()

homo_info.shape

"""##### removing sample info"""

del homo_info['FID']
del homo_info['IID']
homo_info.head()

homo_info_control = homo_info.loc[homo_info['PHENOTYPE'] == 1]
print(homo_info_control.shape)
homo_info_case = homo_info.loc[homo_info['PHENOTYPE'] == 2]
print(homo_info_case.shape)

homo_info.to_csv('PTRHD1.WGS.CODING.VARIANTS.recode.homo.txt', sep = '\t', index=False)

homo_info_control.to_csv('CONTROLS_PTRHD1.WGS.CODING.VARIANTS.recode.homo.txt', sep = '\t', index=False)

homo_info_case.to_csv('CASES_PTRHD1.WGS.CODING.VARIANTS.recode.homo.txt', sep = '\t', index=False)

"""### Identifying patients heterozygous for an alt allele"""

het = df_sub[df_sub.values == 1]
het = het.drop_duplicates()
print(het.shape)
het.head()

het.columns

het['count'] = het[['rs149553409_C', 'rs202008308_G', 'rs148155916_C',
       'rs142018995_T', 'rs149214507_C', 'rs146794228_A', 'rs201790363_T',
       'rs145705009_A', 'rs138104918_C', 'rs149490066_G', 'rs370093383_G',
       'rs56099330_A', 'rs1804645_T', 'rs150066931_G', 'rs148719445_C',
       'rs201252444_C', 'rs138426814_C', 'rs201704717_C', 'rs200722485_A',
       'rs143782986_T', 'rs1550116_G', 'rs201920035_A', 'rs76741040_G',
       'rs147672058_G', 'rs149340371_T', 'rs148113144_C', 'rs146239018_T',
       'rs139463930_T', 'rs143399833_A', 'rs146553503_T', 'rs147179390_T',
       'rs199524620_A', 'rs146165057_T', 'rs185180064_T', 'rs190125602_A',
       'rs78678013_T', 'rs115329263_A', 'rs200642280_A', 'rs199731805_T',
       'rs139407103_A', 'rs114530292_C', 'rs372743845_T', 'rs201606553_T',
       'rs149726253_A']].sum(axis=1)
het.head()

"""#### finding samples with more than one alt allele"""

het_comp = het.loc[het['count'] > 1]
het_comp = het_comp.drop_duplicates()
print(het_comp.shape)
het_comp.head()

new_df = pd.merge(het_comp, info,  how='left', left_on=['FID','IID'], right_on = ['FID','IID'])
new_df = new_df[['FID','IID','PAT','MAT','SEX','PHENOTYPE','rs149553409_C', 'rs202008308_G', 'rs148155916_C',
       'rs142018995_T', 'rs149214507_C', 'rs146794228_A', 'rs201790363_T',
       'rs145705009_A', 'rs138104918_C', 'rs149490066_G', 'rs370093383_G',
       'rs56099330_A', 'rs1804645_T', 'rs150066931_G', 'rs148719445_C',
       'rs201252444_C', 'rs138426814_C', 'rs201704717_C', 'rs200722485_A',
       'rs143782986_T', 'rs1550116_G', 'rs201920035_A', 'rs76741040_G',
       'rs147672058_G', 'rs149340371_T', 'rs148113144_C', 'rs146239018_T',
       'rs139463930_T', 'rs143399833_A', 'rs146553503_T', 'rs147179390_T',
       'rs199524620_A', 'rs146165057_T', 'rs185180064_T', 'rs190125602_A',
       'rs78678013_T', 'rs115329263_A', 'rs200642280_A', 'rs199731805_T',
       'rs139407103_A', 'rs114530292_C', 'rs372743845_T', 'rs201606553_T',
       'rs149726253_A', 'count']]

new_df = new_df.dropna()
new_df.head()

print(new_df.shape)

new_df.head()

"""#### removing sample info"""

del new_df['FID']
del new_df['IID']

controls = new_df.loc[new_df['PHENOTYPE'] == 1]
cases = new_df.loc[new_df['PHENOTYPE'] == 2]

print(controls.shape)
print(cases.shape)

new_df.to_csv('PTRHD1.WGS.CODING.VARIANTS.recode.compound.hets.txt', sep = '\t', index=True)
new_df.head()

cases.to_csv('CASES_PTRHD1.WGS.CODING.VARIANTS.recode.compound.hets.txt', sep = '\t', index=True)

controls.to_csv('CONTROLS_PTRHD1.WGS.CODING.VARIANTS.recode.compound.hets.txt', sep = '\t', index=True)
```
