# ANG Variant Analysis

Francis Grenn (NIH), Cornelis Blauwendraat (NIH)
---
#### positions:
```
chrm = 14
#hg19
start = 21152336 
end = 21167130
#hg38
start_hg38=20684177
end_hg38=20698971
```
## 1) Analysis in IPDGC Cohort
### A) Subset plink binaries
```
plink --bfile $PATH/HARDCALLS_PD_september_2018_no_cousins --remove-fam $PATH/NeuroX.fID.txt --chr 14  --from-bp 21152336 --geno 0.15 --to-bp 21167130 --make-bed --out $PATH/hardcallsNoNeuroX/ANG.GWAS
```
### B) Convert to CVF
```
plink --bfile $PATH/hardcallsNoNeuroX/ANG.GWAS --recode vcf-fid --out $PATH/hardcallsNoNeuroX/vcf/ANG.GWAS
```
### C) Output Frequency
```
plink --bfile $PATH/HARDCALLS_PD_september_2018_no_cousins --chr 14  --from-bp 21152336 --to-bp 21167130  --remove-fam $PATH/NeuroX.fID.txt --freq --geno 0.15 --out $PATH/freq/ANG
plink --bfile $PATH/HARDCALLS_PD_september_2018_no_cousins --chr 14  --from-bp 21152336 --to-bp 21167130  --remove-fam $PATH/NeuroX.fID.txt --freq case-control --geno 0.15 --out $PATH/freq/ANG
```
### D) Output Logistic Regression
```
plink --bfile $PATH/HARDCALLS_PD_september_2018_no_cousins --remove-fam $PATH/NeuroX.fID.txt --chr 14  --from-bp 21152336 --to-bp 21167130 --geno 0.15 --covar $PATH/IPDGC_all_samples_covariates.tab --covar-name sex,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --assoc --out $PATH/logistic/ANG
```
### E) Annotate VCFs with ANNOVAR
```
module load annovar
table_annovar.pl $PATH/vcf/ANG.GWAS.vcf.gz $PATH/annovar/2019-10-24/hg19/ -buildver hg19 --thread 16 -out $PATH/annotation/ANG.annovar -remove -protocol avsnp147,refGene,ensGene,gnomad211_genome -operation f,g,g,f -nastring . -vcfinput

head -1 $PATH/hardcallsNoNeuroX/annotation/ANG.annovar.hg19_multianno.txt > $PATH/hardcallsNoNeuroX/annotation/header.txt
colct="$(wc -w $PATH/hardcallsNoNeuroX/annotation/header.txt| cut -f1 -d' ')"
cut -f1-$colct $PATH/hardcallsNoNeuroX/annotation/ANG.annovar.hg19_multianno.txt > $PATH/hardcallsNoNeuroX/annotation/ANG.trimmed.annotation.txt
```
### F) Burden Analysis
```
#no upper freq bound
rvtest --noweb --inVcf $PATH/hardcallsNoNeuroX/vcf/ANG.GWAS.vcf.gz --pheno $PATH/IPDGC_all_samples_covariates.tab --covar $PATH/IPDGC_all_samples_covariates.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile $PATH/hardcallsNoNeuroX/burden/refFlat_hg19_nochr.txt --out $PATH/hardcallsNoNeuroX/burden/BURDEN.ANG.mafnobound

#0.05 freq bound
rvtest --noweb --inVcf $PATH/hardcallsNoNeuroX/vcf/ANG.GWAS.vcf.gz --pheno $PATH/IPDGC_all_samples_covariates.tab --covar $PATH/IPDGC_all_samples_covariates.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile $PATH/hardcallsNoNeuroX/burden/refFlat_hg19_nochr.txt --freqUpper 0.05 --out $PATH/hardcallsNoNeuroX/burden/BURDEN.ANG.maf05

#0.03 freq bound
rvtest --noweb --inVcf $PATH/hardcallsNoNeuroX/vcf/ANG.GWAS.vcf.gz --pheno $PATH/IPDGC_all_samples_covariates.tab --covar $PATH/IPDGC_all_samples_covariates.tab --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE --burden cmc,zeggini,mb,fp,cmcWald --kernel skat,skato --geneFile $PATH/hardcallsNoNeuroX/burden/refFlat_hg19_nochr.txt --freqUpper 0.03 --out $PATH/hardcallsNoNeuroX/burden/BURDEN.ANG.maf03
```

## 2) Analysis in AMP-PD Data
### A) Fisher exact test
```
plink --bfile $PATH/ANG_WGS_AMP_PD --pheno $PATH/pheno_Mike.txt --make-bed --out $PATH/ANG_WGS_AMP_PD_pheno
plink --bfile $PATH/ANG_WGS_AMP_PD_pheno --update-sex $PATH/toupdatesex.txt --make-bed --out $PATH/ANG_WGS_AMP_PD_pheno_sex
plink --bfile $PATH/ANG_WGS_AMP_PD_pheno_sex --fisher --pheno $PATH/pheno_Mike.txt --covar $PATH/covs_Mike.txt --covar-name SEX, AGE, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 --out $PATH/ANG_WGS_AMP_PD --ci 0.95
```

check for significance
```
#PYTHON
fisher = pd.read_table("$PATH/ANG_WGS_AMP_PD.assoc.fisher", delim_whitespace=True)
print(fisher.shape)

#check lowest p value
print(fisher.sort_values(by='P', ascending=True).head())

#multiple test correction
cutoff = 0.05/ 168
print(cutoff)

multi_test_correct = fisher[fisher['P']< cutoff]
print(multi_test_correct.shape)
```
### B) Annotate
```
plink --bfile $PATH/ANG_WGS_AMP_PD_pheno_sex --recode 'vcf-fid' --out $PATH/ANG_WGS_AMP_PD_pheno_sex
bgzip $PATH/ANG_WGS_AMP_PD_pheno_sex.vcf
tabix -f -p vcf $PATH/ANG_WGS_AMP_PD_pheno_sex.vcf.gz


table_annovar.pl $PATH/ANG_WGS_AMP_PD_pheno_sex.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 -out $PATH/ANG_WGS.annovar -remove -protocol refGene,ljb26_all,gnomad30_genome,clinvar_20190305 -operation g,f,f,f -nastring . -vcfinput

head -1 ANG_WGS.annovar.hg38_multianno.txt > header.txt  
colct="$(wc -w header.txt| cut -f1 -d' ')"  
cut -f1-$colct ANG_WGS.annovar.hg38_multianno.txt > ANG.trimmed.annotation.txt
```

### C) Burden Analysis
make list of coding variants
```
awk '$6=="exonic" {print}' ANG.trimmed.annotation.txt > ANG.trimmed.annotation.coding.variants.WGS.txt
awk '{print $1" "$2" "$2" "$7}' ANG.trimmed.annotation.coding.variants.WGS.txt > ANG.trimmed.annotation.coding.variants.WGS.SNPs.txt

plink --bfile ANG_WGS_AMP_PD_pheno_sex --extract range ANG.trimmed.annotation.coding.variants.WGS.SNPs.txt --recode 'vcf-fid' --out ANG_CODING_AMP

bgzip ANG_CODING_AMP.vcf
tabix -f -p vcf ANG_CODING_AMP.vcf.gz
```

MAF <0.03  
All variants
```
rvtest --noweb --inVcf $PATH/ANG_WGS_AMP_PD_pheno_sex.vcf.gz --pheno $PATH/covs_Mike.txt --covar $PATH/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile $PATH/refFlat_hg38.txt --freqUpper 0.03 --out $PATH/AMP_PD_BURDEN.ANG.maf003
```

MAF < 0.03
coding variants
```
rvtest --noweb --inVcf $PATH/ANG_CODING_AMP.vcf.gz --pheno $PATH/covs_Mike.txt --covar $PATH/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile $PATH/refFlat_hg38.txt --freqUpper 0.03 --out $PATH/AMP_PD_BURDEN.ANG.maf003_CODING
```

### D) Get Frequencies
```
plink --bfile $PATH/ANG_WGS_AMP_PD_pheno_sex --chr 14 --from-bp 20684177 --to-bp 20698971 --freq --geno 0.15 --out $PATH/ANG_WGS_AMP_PD_pheno_sex
plink --bfile $PATH/ANG_WGS_AMP_PD_pheno_sex --chr 14 --from-bp 20684177 --to-bp 20698971 --freq case-control --geno 0.15 --out $PATH/ANG_WGS_AMP_PD_pheno_sex
```

## 3) Locus Zoom
```
module load locuszoom
locuszoom --metal $PATH/nallsEtAl2019_excluding23andMe_allVariants.tab --pvalcol p --markercol SNP --refgene ANG --pop EUR --build hg19 --source 1000G_March2012 --plotonly
locuszoom --metal $PATH/sorted_AAO_april3_18_final_discovery.txt --pvalcol P-value --markercol MarkerName --refgene ANG --pop EUR --build hg19 --source 1000G_March2012 --plotonly
```

## 4) Look at PD Risk and AOO Summary Statistics for ANG

### A) PD Risk
```
#R
library(data.table)
library(dplyr)
library(tidyr)
library(plyr)
data <- fread("$PATH/nallsEtAl2019_excluding23andMe_allVariants.tab", header = T)
data$CHR <- ldply(strsplit(as.character(data$SNP), split = ":"))[[1]]
data$BP <- ldply(strsplit(as.character(data$SNP), split = ":"))[[2]]
genesTemp <- read.table("/data/LNG/Frank/ANG/listofgenes.txt", sep = "", header = F)
colnames(genesTemp) <- c("CHR","START","STOP","GENE")
genes <- subset(genesTemp, CHR != "X")
for(i in 1:length(genes$GENE))
{
  thisGene <- genes$GENE[i]
  thisChr <- genes$CHR[i]
  lower <- genes$START[i]
  upper <- genes$STOP[i]
  output <- subset(data, CHR == thisChr & BP >= lower & BP <= upper)
  fwrite(output, file = paste(thisGene,"_risk_variants.tab", sep = ""), na = "NA", quote = F, row.names = F, sep = "\t")
}

```
### B) AOO
```
#R
library(data.table)
library(dplyr)
library(tidyr)
library(plyr)
data <- fread("$PATH/sorted_AAO_april3_18_final_discovery.txt", header = T)
data$CHR <- ldply(strsplit(as.character(data$MarkerName), split = ":"))[[1]]
data$BP <- ldply(strsplit(as.character(data$MarkerName), split = ":"))[[2]]
genesTemp <- read.table("/data/LNG/Frank/ANG/listofgenes.txt", sep = "", header = F)
colnames(genesTemp) <- c("CHR","START","STOP","GENE")
genes <- subset(genesTemp, CHR != "X")
for(i in 1:length(genes$GENE))
{
  thisGene <- genes$GENE[i]
  thisChr <- genes$CHR[i]
  lower <- genes$START[i]
  upper <- genes$STOP[i]
  output <- subset(data, CHR == thisChr & BP >= lower & BP <= upper)
  fwrite(output, file = paste(thisGene,"_aoo_variants.tab", sep = ""), na = "NA", quote = F, row.names = F, sep = "\t")
}


```

## 5) Compare Protein Sequence
protein sequence obtained from https://uswest.ensembl.org/Homo_sapiens/Transcript/Sequence_Protein?db=core;g=ENSG00000214274;r=14:20684177-20698971;t=ENST00000336811

compare protein sequence to the coding variants AA change we found and to the coding variant AA changes from Van Es et al.

Comparing AA changes because Van Es et al. only identify variants by AA change, not chr:bp or rsid.

Also Van Es et al. AA changes are 24 or 25 positions off from our AA changes, so need to compare to find out why and make sure we are comparing the correct variants in the results. 

```
#PYTHON
protein="MVMGLGVLLLVFVLGLGLTPPTLAQDNSRYTHFLTQHYDAKPQGRDDRYCESIMRRRGLTSPCKDINTFIHGNKRSIKAICENKNGNPHRENLRISKSSFQVTTCKLHGGSPWPPCQYRATAGFRNVVVACENGLPVHLDQSIFRRP"

#make dataframe with AA and position 
pos_df = pd.DataFrame()
i=1
for c in protein:
    pos_df = pos_df.append(pd.DataFrame({'sequence_AA':[c],'sequence_pos':[i]}))
    i = i+1
print(pos_df.shape)
print(pos_df.head())

#read our variants
our_ang_coding_vars = pd.read_csv("$PATH/ANG.trimmed.annotation.coding.variants.WGS.txt",delim_whitespace=True,header=None)
print(our_ang_coding_vars.shape)
print(our_ang_coding_vars.head())

#use regex to extract AA change
our_ang_coding_vars_df = our_ang_coding_vars.iloc[:,10].str.extract(r'p\.(\w)(\d+)(\w)')
our_ang_coding_vars_df.columns = ["our_AA_orig","our_position","our_AA_new"]
our_ang_coding_vars_df["our_position"] = pd.to_numeric(our_ang_coding_vars_df["our_position"])
our_ang_coding_vars_df

#merge
merged_left = pd.merge(left = pos_df, right = our_ang_coding_vars_df, left_on = "sequence_pos", right_on = "our_position", how = "left")
merged_right = pd.merge(left = pos_df, right = our_ang_coding_vars_df, left_on = "sequence_pos", right_on = "our_position", how = "right")

#read previous paper's data
other_data = pd.read_csv("$PATH/VanEs_Data.csv")
print(other_data.shape)
print(other_data.head())
#use regex to extract AA change from their data
other_data[['van_es_AA_orig','van_es_position','van_es_AA_new']] = other_data['Variant'].str.extract(r'(\w)\(?(\-?\d+)\)?(\w)')
#add 25 because it looks like positions are off by 25 or 24
other_data['van_es_position_plus_25'] = pd.to_numeric(other_data['van_es_position']) + 25
other_data_to_merge = other_data[['Variant','van_es_AA_orig','van_es_position','van_es_AA_new','van_es_position_plus_25']]
other_data_to_merge.columns = ['van_es_variant','van_es_AA_orig','van_es_position','van_es_AA_new','van_es_position_plus_25']

#merge all
merged_all = pd.merge(left = merged_left, right = other_data_to_merge, left_on = 'sequence_pos', right_on = 'van_es_position_plus_25', how = 'left')
merged_all.shape

merged_all.to_csv("$PATH/compare_positions.csv",index=None)
```

compare positions of the AA changes in the compare_positions.csv file to identify the correct offset so we know which of our variants match variants from other papers.

seems that all variants from Van Es et al. that have a negative AA change position are offset by 25. And variants from Van Es et al. that have a positive AA change position are offset by 24. So take this into account when generating tables.
