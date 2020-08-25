# UQCRC1

Konstantin Senkevich, Lynne Krohn (McGill)

## GWAS analysis. IPDGC cohort (Nalls et al., 2019).

###  1. Subset PLINK Binaries + Convert to VCFs 

```
UQCRC1 positions on hg19
plink --bfile /data/LNG/saraB/HARDCALLS_PD_september_2018_no_cousins --remove-fam /data/LNG/saraB/NeuroX.fID.txt --geno 0.05 --chr 3 --from-bp 48536432 --to-bp 48747098 --make-bed --out hardcallsNoNeuroX/bin/UQCRC1.GWAS
plink --bfile hardcallsNoNeuroX/bin/UQCRC1.GWAS --recode vcf --out hardcallsNoNeuroX/vcf/UQCRC1.GWAS
    
cd hardcallsNoNeuroX/vcf
bgzip UQCRC1.GWAS.vcf
tabix -f -p vcf UQCRC1.GWAS.vcf.gz
```

###  2. Fisher exact test (to calculate allele frequency) and logistic regression

module load plink samtools rvtests R annovar

```        
plink --bfile hardcallsNoNeuroX/bin/UQCRC1.GWAS --chr 3 --from-bp 48536432 --to-bp 48747098 --fisher --out hardcallsNoNeuroX/freq/UQCRC1   
plink --bfile hardcallsNoNeuroX/bin/UQCRC1.GWAS --chr 3 --from-bp 48536432 --to-bp 48747098 --covar /data/LNG/saraB/IPDGC_all_samples_covariates.tab --covar-name sex,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --logistic --out hardcallsNoNeuroX/logistic/UQCRC1
```
#### 3. Annotation 

Downloading databases 

```
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/ 
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 humandb/ 
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_genome humandb/
```  

module load annovar

```    
cd /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/UQCRC1/hardcallsNoNeuroX/vcf
    
table_annovar.pl UQCRC1.GWAS.vcf.gz /data/LNG/saraB/annovar/humandb/ -buildver hg19 \
--thread 16 \
-out UQCRC1.annovar \
-remove -protocol avsnp147,refGene,ensGene,gnomad211_genome \
-operation f,g,g,f \
-nastring . \
-vcfinput
 
head -1 UQCRC1.annovar.hg19_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct UQCRC1.annovar.hg19_multianno.txt > UQCRC1.GWAS.trimmed.annotation.txt
``` 
  
## WGS analysis. AMP-PD cohort  ([www.amp-pd.org](http://www.amp-pd.org/))

###  1. Fisher exact test to caclucalte allele frequency and logistic regression 

UQCRC1 positions on hg38
    
module load plink samtools rvtests R annovar
  
```   
plink --bfile UQCRC1_WGS_AMP_PD --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --make-bed --out UQCRC1_WGS_AMP_PD_pheno
plink --bfile UQCRC1_WGS_AMP_PD_pheno --update-sex /data/LNG/saraB/AMP_PD/toupdatesex.txt --make-bed --out UQCRC1_WGS_AMP_PD_pheno_sex
plink --bfile UQCRC1_WGS_AMP_PD_pheno_sex --fisher --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --out UQCRC1_WGS_AMP_PD
plink --bfile UQCRC1_WGS_AMP_PD_pheno_sex --logistic --pheno /data/LNG/saraB/AMP_PD/pheno_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX, AGE, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 --out UQCRC1_WGS_AMP_PD --ci 0.95
``` 

###  2. Annotation 

```
plink --bfile UQCRC1_WGS_AMP_PD_pheno_sex --recode 'vcf-fid' --out UQCRC1_WGS_AMP_PD_pheno_sex
bgzip UQCRC1_WGS_AMP_PD_pheno_sex.vcf
tabix -f -p vcf UQCRC1_WGS_AMP_PD_pheno_sex.vcf.gz
    
table_annovar.pl UQCRC1_WGS_AMP_PD_pheno_sex.vcf.gz $ANNOVAR_DATA/hg38 -buildver hg38 \
--thread 16 \
-out UQCRC1_WGS.annovar \
-remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
-operation g,f,f,f \
-nastring . \
-vcfinput
    
head -1 UQCRC1_WGS.annovar.hg38_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct UQCRC1_WGS.annovar.hg38_multianno.txt > UQCRC1.WGS.trimmed.annotation.txt
``` 

### 3. Burden analysis

#### Generating list of coding variants

```
awk '$6=="exonic" {print}' UQCRC1.trimmed.annotation.txt > UQCRC1.trimmed.annotation.coding.variants.WGS.txt
awk '{print $1" "$2" "$2" "$7}' UQCRC1.trimmed.annotation.coding.variants.WGS.txt > UQCRC1.trimmed.annotation.coding.variants.WGS.SNPs.txt
    
plink --bfile UQCRC1_WGS_AMP_PD_pheno_sex --extract range UQCRC1.trimmed.annotation.coding.variants.WGS.txt --recode 'vcf-fid' --out UQCRC1_CODING_AMP
bgzip UQCRC1_CODING_AMP.vcf
tabix -f -p vcf UQCRC1_CODING_AMP.vcf.gz
```

#### Generating list of nonsynonymous variants

```
awk '$9=="nonsynonymous" {print}' UQCRC1.trimmed.annotation.txt > UQCRC1.trimmed.annotation.nonsynonymous.variants.WGS.txt
awk '{print $1" "$2" "$2" "$7}' UQCRC1.trimmed.annotation.nonsynonymous.variants.WGS.txt > UQCRC1.trimmed.annotation.nonsynonymous.variants.WGS.SNPs.txt

plink --bfile UQCRC1_WGS_AMP_PD_pheno_sex --extract range UQCRC1.trimmed.annotation.nonsynonymous.variants.WGS.txt --recode 'vcf-fid' --out UQCRC1_NONSYN_AMP
bgzip UQCRC1_NONSYN_AMP.vcf
tabix -f -p vcf UQCRC1_NONSYN_AMP.vcf.gz
```  

#### MAF < 0.03

#### ALL VARIANTS 

```
rvtest --noweb --inVcf UQCRC1_WGS_AMP_PD_pheno_sex.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 --out AMP_PD_BURDEN.UQCRC1.WGS.maf003
```

#### CODING VARIANTS 

```
rvtest --noweb --inVcf UQCRC1_CODING_AMP.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 --out AMP_PD_BURDEN.UQCRC1.WGS.maf003_CODING
```

#### NONSYNONYMOUS VARIANTS 

```
rvtest --noweb --inVcf UQCRC1_NONSYN_AMP.vcf.gz --pheno /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar /data/LNG/saraB/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile /data/LNG/makariousmb/refFlat_hg38.txt --freqUpper 0.03 --out AMP_PD_BURDEN.UQCRC1.WGS.maf003_NONSYN
```