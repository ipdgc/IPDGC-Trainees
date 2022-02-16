# Private-variants data processing
Jing Hu (Northwestern University & Zhejiang University )

## Get private Lost of function variants for sequencing data

### 1. Get just variants in chr1 and chr6 form vcf file (the chormosomes of *PINK1, PARK7 and PRKN*)
```
module load bcftools

bcftools view -r 1,6 -Oz -o chr1_6.file.vcf.gz file.vcf.gz
```

### 2. Annotation with annovar (change coordinates accordingly)

```
module load bcftools/1.4-6
module load perl/5.26

perl convert2annovar.pl -format vcf4 chr1_6.file.vcf.gz -outfile chr1_6.file.avinput --allsample -withfreq -includeinfo

awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $11 }' chr1_6.file.avinput  > short.chr1_6.file.avinput
perl table_annovar.pl short.chr1_6.file.avinput  humandb/ -buildver hg38 -out chr1_6.file -remove -protocol refGene,genomicSuperDups,gnomad211_exome,gnomad30_genome  -operation g,r,f,f -nastring . -polish
```

### 3. Selection of rare missense, stopgain, and nonsynonymous variants

```
# Replace all the space in multianno.txt with ~
cat chr1_6.file.hg38_multianno.txt | sed -e 's/ /~/g' > nospace.chr1_6.file.hg38_multianno.txt

# Add header
echo "chr    start    end    ref    alt    frequency-in-sample    id    " > header.txt
cat header.txt short.chr1_6.file.avinput > addheader.short.chr1_6.file.avinput
paste addheader.short.chr1_6.file.avinput nospace.chr1_6.file.hg38_multianno.txt > withfre.nospace.chr1_6.file.hg38_multianno.txt

# Get variants in three genes
cat withfre.nospace.chr1_6.file.hg38_multianno.txt | awk '$14 == "PRKN" || $14 == "PINK1" || $14 == "PARK7" {print}'> 3gene.annotation.txt

# Get variants with maf < 0.01 in both genome and exome and pass genomicSuperDups
cat 3gene.annotation.txt | awk '$19 < 0.01 || $19 == "." {print}' | awk '$36 < 0.01 || $36 == "." {print}' | awk '$18 == "." {print}' > MAF0.01.3gene.annotation.txt

# Get coding splicing variants with MAF<0.01
cat MAF0.01.3gene.annotation.txt | awk '$13 == "exonic;splicing" || $13 == "splicing" {print}' > splicing.MAF0.01.3gene.annotation.txt

# Get missence variants with MAF<0.01
cat MAF0.01.3gene.annotation.txt | awk '$16 == "nonsynonymous~SNV" {print}' > nonsynonymous.MAF0.01.3gene.annotation.txt

# Get stopgain variants with MAF<0.01 
cat MAF0.01.3gene.annotation.txt | awk '$16 == "stopgain" {print}' > stopgain.MAF0.01.3gene.annotation.txt

# Merge nonsynonymous, stopgain and coding splicing variants and get their ID
cat nonsynonymous.MAF0.01.3gene.annotation.txt stopgain.MAF0.01.3gene.annotation.txt splicing.MAF0.01.3gene.annotation.txt | sort| uniq > nonsy_stopgain_splicing.MAF0.01.3gene.annotation.txt
awk '{print $7}' nonsy_stopgain_splicing.MAF0.01.3gene.annotation.txt > nonsy_stopgain_splicing.MAF0.01.3gene.id.txt

# Get nonsy_stopgain_splicing.MAF0.01.3gene variants form original chr1_6 vcf file

module load bcftools
bcftools view -i 'ID=@nonsy_stopgain_splicing.MAF0.01.3gene.id.txt' -Oz -o  nonsy_stopgain_splicing.MAF0.01.3gene.chr1_6.file.vcf.gz chr1_6.file.vcf.gz
```

### 4. Get private variants (only singletons) 

```
module load vcftools 

vcftools --gzvcf nonsy_stopgain_splicing.MAF0.01.3gene.chr1_6.file.vcf.gz --singletons --out nonsy_stopgain_splicing.MAF0.01.3gene.chr1_6

# Generate information to grep full variants ID in nonsy_stopgain_splicing.MAF0.01.3gene.chr1_6.file.vcf.gz
awk '{print "|""chr"$1"_"$2}' nonsy_stopgain_splicing.MAF0.01.3gene.chr1_6.file.singletons > nonsy_stopgain_splicing.MAF0.01.3gene.chr1_6.id

# Generate the grep command 
awk BEGIN{RS=EOF}'{gsub(/\n/," ");print}' nonsy_stopgain_splicing.MAF0.01.3gene.chr1_6.id > new.nonsy_stopgain_splicing.MAF0.01.3gene.chr1_6.id
sed 's/ //g' new.nonsy_stopgain_splicing.MAF0.01.3gene.singletons.id > command.txt
zcat nonsy_stopgain_splicing.MAF0.01.3gene.chr1_6.file.vcf.gz | grep -E 'command.txt' | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' > check_private.nonsy_stopgain_splicing.MAF0.01.3gene.txt

### some variants share the same position, carefully check for the correct variants and remove the wrong one. Also remove doubletons.

# Get only private variants in vcf file and convert to plink bfile

module load bcftools

bcftools view -i 'ID=@final.private.nonsy_stopgain_splicing.MAF0.01.3gene.id.txt' -Oz -o  private.nonsy_stopgain_splicing.MAF0.01.3gene.chr1_6.file.vcf.gz nonsy_stopgain_splicing.MAF0.01.3gene.chr1_6.file.vcf.gz

/projects/b1049/genetics_programs/plink_Jun_2019/plink --vcf  private.nonsy_stopgain_splicing.MAF0.01.3gene.chr1_6.file.vcf.gz --allow-no-sex --double-id --vcf-half-call missing  --make-bed --out  private.nonsy_stopgain_splicing.MAF0.01.3gene
```

### 5. CADD annotation (upload to CADD website)

```
# Get short VCF file for CADD website annotation (CADD website can annotate indels)
zcat  private.nonsy_stopgain_splicing.MAF0.01.3gene.chr1_6.file.vcf.gz | cut -f 1,2,3,4,5 > short. private.nonsy_stopgain_splicing.MAF0.01.3gene.chr1_6.file.vcf

# Upload to CADD website and get CADD > 20 variants
awk '$6 > 20 {print "chr"$1"_"$2"_"$3"_"$4 }' *.tsv  | grep -v "#" | sort | uniq >  uniq.CADD20.private.nonsy_stopgain_splicing.MAF0.01.3gene.txt

# Get only private variants with CADD > 20
/projects/b1049/genetics_programs/plink_Jun_2019/plink --bfile private.nonsy_stopgain_splicing.MAF0.01.3gene --extract uniq.CADD20.private.nonsy_stopgain_splicing.MAF0.01.3gene.id.txt --make-bed --allow-no-sex --out private.CADD20.nonsy_stopgain_splicing.MAF0.01.3gene

# check the row numbers of CADD20.private.nonsy_stopgain_splicing.MAF0.01.3gene.txt and private.CADD20.nonsy_stopgain_splicing.MAF0.01.3gene.bim file to see if we get all the variants (if different, some variants have reversed ref and alt)
wc -l private.CADD20.nonsy_stopgain_splicing.MAF0.01.3gene.bim
wc -l uniq.CADD20.private.nonsy_stopgain_splicing.MAF0.01.3gene.txt

# If the numbers are not the same, do:
awk '{print $2}' private.CADD20.nonsy_stopgain_splicing.MAF0.01.3gene.bim > rightseq.private.CADD20.nonsy_stopgain_splicing.MAF0.01.3gene.id

sort uniq.CADD20.private.nonsy_stopgain_splicing.MAF0.01.3gene.id.txt rightseq.private.CADD20.nonsy_stopgain_splicing.MAF0.01.3gene.id rightseq.private.CADD20.nonsy_stopgain_splicing.MAF0.01.3gene.id | uniq -u > reverseseq.rightseq.private.CADD20.nonsy_stopgain_splicing.MAF0.01.3gene.id

sed 's/_/ /g' reverseseq.rightseq.private.CADD20.nonsy_stopgain_splicing.MAF0.01.3gene.id | awk '{print $1"_"$2"_"$4"_"$3}' > change.reverseseq.rightseq.private.CADD20.nonsy_stopgain_splicing.MAF0.01.3gene.id

cat rightseq.private.CADD20.nonsy_stopgain_splicing.MAF0.01.3gene.id  change.reverseseq.rightseq.private.CADD20.nonsy_stopgain_splicing.MAF0.01.3gene.id > final.private.CADD20.nonsy_stopgain_splicing.MAF0.01.3gene.id

# Get final only private variants with CADD >20
/projects/b1049/genetics_programs/plink_Jun_2019/plink --bfile private.nonsy_stopgain_splicing.MAF0.01.3gene --extract final.private.CADD20.nonsy_stopgain_splicing.MAF0.01.3gene.id --make-bed --allow-no-sex --out private.CADD20.nonsy_stopgain_splicing.MAF0.01.3gene

# check again if we grep all the CADD > 20 variants, should have the same row number now
wc -l private.CADD20.nonsy_stopgain_splicing.MAF0.01.3gene.bim
wc -l uniq.CADD20.private.nonsy_stopgain_splicing.MAF0.01.3gene.txt
```

#### For the summary counts data of private variants extracted from Parkinson's Disease variant browser, we also filtered them by functional annotation and CADD annotation similar to the above steps.

## SKAT-O script for individual cohorts

```
library(SKAT)
#### AMP-PD-WGS whole cohort
File.Bed <- "private.nonsy_stopgain_splicing.MAF0.01.3gene.bed"
File.Bim <- "private.nonsy_stopgain_splicing.MAF0.01.3gene.bim"
File.Fam <- "private.nonsy_stopgain_splicing.MAF0.01.3gene.fam"
File.Cov <- "new.AMP-PD_April2021.Covars_Pheno.txt"
File.SetID <- "AMP-PD.private.nonsy_stopgain_splicing.MAF0.01.3gene.SetID"
File.SSD <- "AMP-PD.private.nonsy_stopgain_splicing.MAF0.01.3gene.SSD"
File.Info <- "AMP-PD.private.nonsy_stopgain_splicing.MAF0.01.3gene.Info"

Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)
FAM_Cov <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary = TRUE, cov_header = TRUE)

#open SSD file
SSD.INFO <- Open_SSD(File.SSD, File.Info)

#Run with covariates
obj <- SKAT_Null_Model(PHENO ~ SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, out_type = "D", data = FAM_Cov, Adjustment = FALSE)
out.skato <- SKATBinary.SSD.All(SSD.INFO, obj, method = "SKATO")
out.df = out.skato$results
write.table(out.df, file = "AMP-PD.private.nonsy_stopgain_splicing.MAF0.01.3gene.skato.txt", col.names = TRUE, row.names = FALSE)
quit()
n

#### IPDGC-WES whole cohort
library(SKAT)
File.Bed <- "private.maf0.01.3genes.nonsy_stopgain_splicing.new.IPDGC.chr1-22.bed"
File.Bim <- "private.maf0.01.3genes.nonsy_stopgain_splicing.new.IPDGC.chr1-22.bim"
File.Fam <- "private.maf0.01.3genes.nonsy_stopgain_splicing.new.IPDGC.chr1-22.fam"
File.Cov <- "RVAS.covariate.txt"
File.SetID <- "private.maf0.01.3genes.nonsy_stopgain_splicing.new.IPDGC.chr1-22.SetID"
File.SSD <- "private.maf0.01.3genes.nonsy_stopgain_splicing.new.IPDGC.chr1-22.SSD"
File.Info <- "private.maf0.01.3genes.nonsy_stopgain_splicing.new.IPDGC.chr1-22.Info"

Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)
FAM_Cov <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary = TRUE, cov_header = TRUE)

#open SSD file
SSD.INFO <- Open_SSD(File.SSD, File.Info)

#Run with covariates
obj <- SKAT_Null_Model(Phenotype.x ~ Sex.x + covarage + pc1 + pc2 +pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, out_type = "D", data = FAM_Cov, Adjustment = FALSE)
out.skato <- SKATBinary.SSD.All(SSD.INFO, obj, method = "SKATO")
out.df = out.skato$results
write.table(out.df, file = "private.maf0.01.3genes.nonsy_stopgain_splicing.new.IPDGC.chr1-22.skato.txt", col.names = TRUE, row.names = FALSE)
quit()
n

#### McGill ISR whole cohort
library(SKAT)
File.Bed <- "private.nonsy_stopgain_splicing.3genes.ISR.DP50.bed"
File.Bim <- "private.nonsy_stopgain_splicing.3genes.ISR.DP50.bim"
File.Fam <- "private.nonsy_stopgain_splicing.3genes.ISR.DP50.fam"
File.Cov <- "recode.noNA.PINK1.ISR.FullSamples.reduced.txt"
File.SetID <- "private.nonsy_stopgain_splicing.3genes.ISR.DP50.SetID"
File.SSD <- "private.nonsy_stopgain_splicing.3genes.ISR.DP50.SSD"
File.Info <- "private.nonsy_stopgain_splicing.3genes.ISR.DP50.Info"

Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)
FAM_Cov <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary = TRUE, cov_header = TRUE)

#open SSD file
SSD.INFO <- Open_SSD(File.SSD, File.Info)

#Run with covariates
obj <- SKAT_Null_Model(PHENO ~ SEX + Age + Ethn + GBA_status, out_type = "D", data = FAM_Cov, Adjustment = FALSE)
out.skato <- SKATBinary.SSD.All(SSD.INFO, obj, method = "SKATO")
out.df = out.skato$results
write.table(out.df, file = "private.nonsy_stopgain_splicing.3genes.ISR.DP50.skato.txt", col.names = TRUE, row.names = FALSE)
quit()
n

#### McGill FC whole cohort
library(SKAT)
File.Bed <- "private.nonsy_stopgain_splicing.3genes.FC.bed"
File.Bim <- "private.nonsy_stopgain_splicing.3genes.FC.bim"
File.Fam <- "private.nonsy_stopgain_splicing.3genes.FC.fam"
File.Cov <- "recode.noNA.PINK1.FC.FullSamples.reduced.txt"
File.SetID <- "private.nonsy_stopgain_splicing.3genes.FC.DP50.SetID"
File.SSD <- "private.nonsy_stopgain_splicing.3genes.FC.DP50.SSD"
File.Info <- "private.nonsy_stopgain_splicing.3genes.FC.DP50.Info"

Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)
FAM_Cov <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary = TRUE, cov_header = TRUE)

#open SSD file
SSD.INFO <- Open_SSD(File.SSD, File.Info)

#Run with covariates
obj <- SKAT_Null_Model(PHENO ~ SEX + Age + Ethn + GBA_status, out_type = "D", data = FAM_Cov, Adjustment = FALSE)
out.skato <- SKATBinary.SSD.All(SSD.INFO, obj, method = "SKATO")
out.df = out.skato$results
write.table(out.df, file = "private.nonsy_stopgain_splicing.3genes.FC.DP50.skato.txt", col.names = TRUE, row.names = FALSE)
quit()
n

#### McGill NY whole cohort
library(SKAT)
File.Bed <- "private.nonsy_stopgain_splicing.3genes.NY.bed"
File.Bim <- "private.nonsy_stopgain_splicing.3genes.NY.bim"
File.Fam <- "private.nonsy_stopgain_splicing.3genes.NY.fam"
File.Cov <- "recode.noNA.PINK1.NY.FullSamples.reduced.txt"
File.SetID <- "private.nonsy_stopgain_splicing.3genes.NY.DP50.SetID"
File.SSD <- "private.nonsy_stopgain_splicing.3genes.NY.DP50.SSD"
File.Info <- "private.nonsy_stopgain_splicing.3genes.NY.DP50.Info"

Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)
FAM_Cov <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary = TRUE, cov_header = TRUE)

#open SSD file
SSD.INFO <- Open_SSD(File.SSD, File.Info)

#Run with covariates
obj <- SKAT_Null_Model(PHENO ~ SEX + Age + Ethn + GBA_status, out_type = "D", data = FAM_Cov, Adjustment = FALSE)
out.skato <- SKATBinary.SSD.All(SSD.INFO, obj, method = "SKATO")
out.df = out.skato$results
write.table(out.df, file = "private.nonsy_stopgain_splicing.3genes.NY.DP50.skato.txt", col.names = TRUE, row.names = FALSE)
quit()
n

```
## META-SKAT script

```
library(MetaSKAT)
library(SKAT)

### ######################
### cohort 1
### ######################
File1.Bed <- "private.nonsy_stopgain_splicing.MAF0.01.3gene.bed"
File1.Bim <- "private.nonsy_stopgain_splicing.MAF0.01.3gene.bim"
File1.Fam <- "private.nonsy_stopgain_splicing.MAF0.01.3gene.fam"
File1.Cov <- "new.AMP-PD_April2021.Covars_Pheno.txt"
File1.SetID <- "AMP-PD.private.nonsy_stopgain_splicing.MAF0.01.3gene.SetID"
File1.MSSD <- "AMP-PD.private.nonsy_stopgain_splicing.MAF0.01.3gene.MSSD"
File1.MInfo <- "AMP-PD.private.nonsy_stopgain_splicing.MAF0.01.3gene.MInfo"

FAM_Cov1 <- Read_Plink_FAM_Cov(File1.Fam, File1.Cov, Is.binary = TRUE, cov_header = TRUE)
obj1 <- SKAT_Null_Model(PHENO ~ SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, out_type = "D", data = FAM_Cov1, Adjustment = FALSE)
#### Generate summary statistics files
N.Sample<- 6053
Generate_Meta_Files(obj1, File1.Bed, File1.Bim, File1.SetID, File1.MSSD, File1.MInfo, N.Sample, File.Permu = NULL, data=NULL, impute.method="fixed")

### ######################
### cohort 2
### ######################
File2.Bed <- "private.maf0.01.3genes.nonsy_stopgain_splicing.new.IPDGC.chr1-22.bed"
File2.Bim <- "private.maf0.01.3genes.nonsy_stopgain_splicing.new.IPDGC.chr1-22.bim"
File2.Fam <- "private.maf0.01.3genes.nonsy_stopgain_splicing.new.IPDGC.chr1-22.fam"
File2.Cov <- "RVAS.covariate.txt"
File2.SetID <- "private.maf0.01.3genes.nonsy_stopgain_splicing.new.IPDGC.chr1-22.SetID"
File2.MSSD <- "private.maf0.01.3genes.nonsy_stopgain_splicing.new.IPDGC.chr1-22.MSSD"
File2.MInfo <- "private.maf0.01.3genes.nonsy_stopgain_splicing.new.IPDGC.chr1-22.MInfo"

FAM_Cov2 <- Read_Plink_FAM_Cov(File2.Fam, File2.Cov, Is.binary = TRUE, cov_header = TRUE)
obj2 <- SKAT_Null_Model(Phenotype.x ~ Sex.x + covarage + pc1 + pc2 +pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, out_type = "D", data = FAM_Cov2, Adjustment = FALSE)

#### Generate summary statistics files
N.Sample<- 1564
Generate_Meta_Files(obj2, File2.Bed, File2.Bim, File2.SetID, File2.MSSD, File2.MInfo, N.Sample, File.Permu = NULL, data=NULL, impute.method="fixed")

### ######################
### cohort 3
### ######################
File3.Bed <- "private.nonsy_stopgain_splicing.3genes.ISR.DP50.bed"
File3.Bim <- "private.nonsy_stopgain_splicing.3genes.ISR.DP50.bim"
File3.Fam <- "private.nonsy_stopgain_splicing.3genes.ISR.DP50.fam"
File3.Cov <- "recode.noNA.PINK1.ISR.FullSamples.reduced.txt"
File3.SetID <- "private.nonsy_stopgain_splicing.3genes.ISR.DP50.SetID"
File3.MSSD <- "private.nonsy_stopgain_splicing.3genes.ISR.DP50.MSSD"
File3.MInfo <- "private.nonsy_stopgain_splicing.3genes.ISR.DP50.MInfo"

FAM_Cov3 <- Read_Plink_FAM_Cov(File3.Fam, File3.Cov, Is.binary = TRUE, cov_header = TRUE)
obj3 <- SKAT_Null_Model(PHENO ~ SEX + Age + Ethn + GBA_status, out_type = "D", data = FAM_Cov3, Adjustment = FALSE)

#### Generate summary statistics files
N.Sample<- 1151
Generate_Meta_Files(obj3, File3.Bed, File3.Bim, File3.SetID, File3.MSSD, File3.MInfo, N.Sample, File.Permu = NULL, data=NULL, impute.method="fixed")

### ######################
### cohort 4
### ######################
File4.Bed <- "private.nonsy_stopgain_splicing.3genes.FC.bed"
File4.Bim <- "private.nonsy_stopgain_splicing.3genes.FC.bim"
File4.Fam <- "private.nonsy_stopgain_splicing.3genes.FC.fam"
File4.Cov <- "recode.noNA.PINK1.FC.FullSamples.reduced.txt"
File4.SetID <- "private.nonsy_stopgain_splicing.3genes.FC.DP50.SetID"
File4.MSSD <- "private.nonsy_stopgain_splicing.3genes.FC.DP50.MSSD"
File4.MInfo <- "private.nonsy_stopgain_splicing.3genes.FC.DP50.MInfo"

FAM_Cov4 <- Read_Plink_FAM_Cov(File4.Fam, File4.Cov, Is.binary = TRUE, cov_header = TRUE)
obj4 <- SKAT_Null_Model(PHENO ~ SEX + Age + Ethn + GBA_status, out_type = "D", data = FAM_Cov4, Adjustment = FALSE)

#### Generate summary statistics files
N.Sample<- 3001
Generate_Meta_Files(obj4, File4.Bed, File4.Bim, File4.SetID, File4.MSSD, File4.MInfo, N.Sample, File.Permu = NULL, data=NULL, impute.method="fixed")

### ######################
### cohort 5
### ######################
File5.Bed <- "private.nonsy_stopgain_splicing.3genes.NY.bed"
File5.Bim <- "private.nonsy_stopgain_splicing.3genes.NY.bim"
File5.Fam <- "private.nonsy_stopgain_splicing.3genes.NY.fam"
File5.Cov <- "recode.noNA.PINK1.NY.FullSamples.reduced.txt"
File5.SetID <- "private.nonsy_stopgain_splicing.3genes.NY.DP50.SetID"
File5.MSSD <- "private.nonsy_stopgain_splicing.3genes.NY.DP50.MSSD"
File5.MInfo <- "private.nonsy_stopgain_splicing.3genes.NY.DP50.MInfo"

FAM_Cov5 <- Read_Plink_FAM_Cov(File5.Fam, File5.Cov, Is.binary = TRUE, cov_header = TRUE)
obj5 <- SKAT_Null_Model(PHENO ~ SEX + Age + Ethn + GBA_status, out_type = "D", data = FAM_Cov5, Adjustment = FALSE)

#### Generate summary statistics files
N.Sample<- 1323
Generate_Meta_Files(obj5, File5.Bed, File5.Bim, File5.SetID, File5.MSSD, File5.MInfo, N.Sample, File.Permu = NULL, data=NULL, impute.method="fixed")

#### ####################
#### For meta
#### ####################
File.MSSD.vec<-rep("",5)
File.MInfo.vec<-rep("",5)
File.MSSD.vec[1]<-"AMP-PD.private.nonsy_stopgain_splicing.MAF0.01.3gene.MSSD"
File.MInfo.vec[1]<-"AMP-PD.private.nonsy_stopgain_splicing.MAF0.01.3gene.MInfo"
File.MSSD.vec[2]<-"private.maf0.01.3genes.nonsy_stopgain_splicing.new.IPDGC.chr1-22.MSSD"
File.MInfo.vec[2]<-"private.maf0.01.3genes.nonsy_stopgain_splicing.new.IPDGC.chr1-22.MInfo"
File.MSSD.vec[3]<-"private.nonsy_stopgain_splicing.3genes.ISR.DP50.MSSD"
File.MInfo.vec[3]<-"private.nonsy_stopgain_splicing.3genes.ISR.DP50.MInfo"
File.MSSD.vec[4]<-"private.nonsy_stopgain_splicing.3genes.FC.DP50.MSSD"
File.MInfo.vec[4]<-"private.nonsy_stopgain_splicing.3genes.FC.DP50.MInfo"
File.MSSD.vec[5]<-"private.nonsy_stopgain_splicing.3genes.NY.DP50.MSSD"
File.MInfo.vec[5]<-"private.nonsy_stopgain_splicing.3genes.NY.DP50.MInfo"

### Open_MSSD_File_2Read Read Meta SSD and Info files
Cohort.Info<- Open_MSSD_File_2Read(File.MSSD.vec, File.MInfo.vec)

#### RUN meta analysis
out<- MetaSKAT_MSSD_ALL(Cohort.Info, combined.weight=TRUE, weights.beta=c(1,25), method="optimal", r.corr=0, is.separate = FALSE, Group_Idx=NULL, MAF.cutoff=1, missing_cutoff=0.15)
```

## Fisher's exact test 

```
data <- read.table("all.forfisher.txt", header = FALSE)
for (i in 1:nrow(data)) {
row <- data[i,]
dataM <- matrix(c(row$V2, row$V3, row$V4, row$V5),nrow = 2)
Pvalues<-fisher.test(dataM, alternative="greater", conf.level = 0.95)$p.value
OR <- fisher.test(dataM)$estimate
CI95_L <- "["(tail((fisher.test(dataM)$conf.int),2), 1)

CI95_U <- "["(tail((fisher.test(dataM)$conf.int),2), 2)

write.table(Pvalues, file="greaterall.fisher_Pvalues.txt", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(OR, file="greaterall.fisher_ORvalues.txt", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(CI95_L, file="greaterall.fisher_CI95_L.txt", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(CI95_U, file="greaterall.fisher_CI95_U.txt", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
}
quit()
n

echo "Gene  CarriersInCases  CarriersInControls  NonCarriersInCases  NonCarriersInControls  p_value  OR  CI95_L  CI95_U" > fisher_result.header.txt 

paste all.forfisher.txt greaterall.fisher_Pvalues.txt greaterall.fisher_ORvalues.txt greaterall.fisher_CI95_L.txt greaterall.fisher_CI95_U.txt > greaterall.fisher_wPvalues.txt

cat fisher_result.header.txt greaterall.fisher_wPvalues.txt > withheader.greaterall.fisher_wPvalues.txt
```

## Generate forest plot using Metafor

```
library(metafor)
### copy meta-analysis data into 'dat'
dat <- read.table("newall.PARK7.txt", header = TRUE)
 
### calculate log risk ratios and corresponding sampling variances (and use the 'slab' argument to store study labels as part of the data frame)
dat <- escalc(measure="OR", ai=CarriersInCases, bi=NonCarriersInCases, ci=CarriersInControls, di=NonCarriersInControls, data=dat, slab=paste(Gene, sep=", "))

### fit fixed-effects model analysis
res <- rma(yi, vi, data=dat, method="FE")
print(res, digits=3)

### forest plot with extra annotations
pdf(file="/projects/b1049/jing/final.all.PARK7.pdf")
forest(res, atransf=exp, at=log(c(.05, .25, 1, 4)), xlim=c(-16,6),
       ilab=cbind(dat.bcg$CarriersInCases, dat.bcg$NonCarriersInCases, dat.bcg$CarriersInControls, dat.bcg$NonCarriersInControls),
       ilab.xpos=c(-9.5,-8,-6,-4.5), cex=.75, header="Gene and Cohort",
       mlab="")
 
### add text with Q-value, dfs, p-value, and I^2 statistic
text(-16, -1, pos=4, cex=0.75, bquote(paste("Test for Heterogeneity (Q = ",
     .(formatC(res$QE, digits=2, format="f")), ", df = ", .(res$k - res$p),
     ", p = ", .(formatC(res$QEp, digits=2, format="f")), "; ", I^2, " = ",
     .(formatC(res$I2, digits=1, format="f")), "%)")))	

text(-16, -1.5, pos=4, cex=0.75, bquote(paste("FE Model (Z = ",
     .(formatC(res$zval, digits=2, format="f")), 
	 ", p = ", .(formatC(res$pval, digits=2, format="f")), ")")))	
     
quit()
n

```