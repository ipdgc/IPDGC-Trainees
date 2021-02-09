
# KTN1 GWAS and whole genome analysis
## 1. GWAS (IPDGC genotype data)
#### 1. Subset PLINK Binaries + Convert to VCFs
#### 2. Output Frequency with PLINK
#### 3. Logistic Regression with PLINK
#### 4. Annotate VCFs with Annovar 

    #!/bin/bash

    #hg19
    
    ##set up environment
    module load plink samtools rvtests R annovar

    mkdir $KTN1_DIR/
    cd $KTN1_DIR/

    mkdir hardcallsNoNeuroX hardcallsNoNeuroX/bin 
    mkdir hardcallsNoNeuroX/freq hardcallsNoNeuroX/score 
    mkdir hardcallsNoNeuroX/burden hardcallsNoNeuroX/vcf 
    mkdir hardcallsNoNeuroX/annotation hardcallsNoNeuroX/burden 
    mkdir hardcallsNoNeuroX/logistic hardcallsNoNeuroX/freq

    # ==================== 1. Subset PLINK Binaries to cis-KTN1 region ====================
    ## Remove NeuroX + keep males + genotype quality of 15% or less missingness
    
    ##KTN1 region 56025790 - 56168244 +/- MB
    ##hg19 coordinates

    plink --bfile $FILE_DIR/HARDCALLS_PD_september_2018_no_cousins \
    --remove-fam $FILE_DIR/NeuroX.fID.txt \
    --chr 14 --geno 0.15 --from-bp 55025790 --to-bp 57168244 \
    --make-bed --out hardcallsNoNeuroX/bin/KTN1.GWAS.MB
        # Done 11.01.2020
        # PLINK output:
            # Total genotyping rate in remaining samples is 0.461474.
            # 15399 variants removed due to missing genotype data (--geno).
            # 6188 variants and 32338 people pass filters and QC.
            # Among remaining phenotypes, 14671 are cases and 17667 are controls.


    plink --bfile hardcallsNoNeuroX/bin/KTN1.GWAS.MB --recode vcf --out hardcallsNoNeuroX/vcf/KTN1.GWAS.MB

    cd hardcallsNoNeuroX/vcf
    bgzip KTN1.GWAS.MB.vcf
    tabix -f -p vcf KTN1.GWAS.MB.vcf.gz

    # ==================== 2. Output Frequency via PLINK ==================== 
    cd $KTN1_DIR/
    # Overall
    plink --bfile $FILE_DIR/HARDCALLS_PD_september_2018_no_cousins \
    --chr 14 --from-bp 55025790 --to-bp 57168244 \
    --remove-fam $FILE_DIR/NeuroX.fID.txt \
    --geno 0.15 --freq \
    --out hardcallsNoNeuroX/freq/KTN1.MB
        # Done 11.01.2021


    # Case-control
    plink --bfile $FILE_DIR/HARDCALLS_PD_september_2018_no_cousins \
    --chr 14 --from-bp 55025790 --to-bp 57168244 \
    --remove-fam $FILE_DIR/NeuroX.fID.txt \
    --geno 0.15 --freq case-control  \
    --out hardcallsNoNeuroX/freq/KTN1.MB
        # Done 11.01.2021

    # ==================== 3. Output Logistic Regression via PLINK ====================
    plink --bfile $FILE_DIR/HARDCALLS_PD_september_2018_no_cousins \
    --remove-fam $FILE_DIR/NeuroX.fID.txt \
    --chr 14 --from-bp 55025790 --to-bp 57168244 \
    --geno \
    --covar $FILE_DIR/IPDGC_all_samples_covariates.tab \
    --covar-name sex,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE \
    --assoc \
    --out hardcallsNoNeuroX/logistic/KTN1
        # Done 12.01.2020

    # ==================== 4. Annotate VCFs via Annovar ====================

    # Downloading databases 
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
    # annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/ 
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 humandb/ 
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a humandb/
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_genome humandb/

    table_annovar.pl $KTN1_DIR/hardcallsNoNeuroX/vcf/KTN1.GWAS.MB.vcf.gz $FILE_DIR/annovar/humandb/ -buildver hg19 \
    --thread 16 \
    -out $KTN1_DIR/hardcallsNoNeuroX/annotation/KTN1.MB.annovar \
    -remove -protocol avsnp147,refGene,ensGene,gnomad211_genome \
    -operation f,g,g,f \
    -nastring . \
    -vcfinput

    cd $KTN1_DIR/hardcallsNoNeuroX/annotation
    head -1 KTN1.MB.annovar.hg19_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct KTN1.MB.annovar.hg19_multianno.txt > KTN1.MB.trimmed.annotation.txt
        # Done 12.01.2020


## 2. Whole Genome Analysis (AMP-PD)
#### 1. Extract KTN1 Region with PLINK
#### 2. Annotate with Annovar
#### 3. Burden analysis with RVTESTS
#### 4. Single Variant Wald and Score Tests with RVTESTS 
 
 
    #!/bin/bash
    
    ##set up environment
    module load plink samtools rvtests annovar

    ## setting up directories
    cd $KTN1_DIR/

    mkdir AMPPD

    # Initializing some variables 
    COV_NAMES="SEX,AGE,FAMILY_HISTORY,EDUCATION,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,STUDY_BioFIND,STUDY_PDBP,STUDY_PPMI"

    ## ======================= 0. Prep file  ========================= 
    ##whole data exported from gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/plink/bfile/all_vcfs*

    ##add pheno data to data
    plink --bfile $KTN1_DIR/amppd/all_vcfs --pheno $FILE_DIR/AMP_PD/pheno_Mike.txt --make-bed --out $KTN1_DIR/amppd/plink/WGS_AMP_PD_pheno
        # Done 07.09.2020
        # PLINK output:
            # 28791969 variants and 3074 people pass filters and QC.
            # Among remaining phenotypes, 1735 are cases and 1100 are controls.  (239
            # phenotypes are missing.)

    ## add sex data to data
    plink --bfile $KTN1_DIR/amppd/plink/WGS_AMP_PD_pheno --update-sex $FILE_DIR/AMP_PD/toupdatesex.txt --make-bed --out $KTN1_DIR/amppd/plink/WGS_AMP_PD_pheno_sex
        # Done 07.09.2020


    # ==================== 1. Extract cis-KTN1 Region via PLINK  ====================

    ##extract KTN1 region
    # Extract the KTN1 region
        # These are hg38 co-ordinates 
        ##extract KTN1 region : 55559072 - 56801526 +/-1Mb
    # no filter
    plink --bfile $KTN1_DIR/amppd/plink/WGS_AMP_PD_pheno_sex \
    --chr 14 --from-bp 55975694 --to-bp 55975694 \
    --make-bed --out $KTN1_DIR1/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_MB
        # Done 11.01.2020
        # PLINK output:
            # Total genotyping rate is 0.995527.
            # 22803 variants and 3074 people pass filters and QC.
            # Among remaining phenotypes, 1735 are cases and 1100 are controls.

    # --mac 1
    plink --bfile $KTN1_DIR/amppd/plink/WGS_AMP_PD_pheno_sex \
    --chr 14 --from-bp 54559072 --to-bp 56701526 \
    --mac 1 \
    --make-bed --out $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_MB_mac1
        # Done 14.01.2020
        # PLINK output:
            # 22789 variants and 3074 people pass filters and QC.
            # Among remaining phenotypes, 1735 are cases and 1100 are controls.


    # --mac 3 
    plink --bfile $KTN1_DIR/amppd/plink/WGS_AMP_PD_pheno_sex \
    --chr 14 --from-bp 54559072 --to-bp 56701526 \
    --mac 3 \
    --make-bed --out $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_MB_mac3 
            # Done 11.01.20
            # PLINK output:
                # 17094 variants and 3074 people pass filters and QC.
                # Among remaining phenotypes, 1735 are cases and 1100 are controls.


    # ==================== 2. Annotation via Annovar ====================

    ## no filters
    plink --bfile KTN1_WGS_AMP_PD_pheno_sex_MB --recode 'vcf-fid' --out KTN1_WGS_AMP_PD_pheno_sex_MB
    bgzip KTN1_WGS_AMP_PD_pheno_sex_MB.vcf
    tabix -f -p vcf KTN1_WGS_AMP_PD_pheno_sex_MB.vcf.gz

    #Now using ANNOVAR

    table_annovar.pl $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_MB.vcf.gz $KTN1_DIR/amppd/humandb -buildver hg38 \
    --thread 16 \
    -out $KTN1_DIR/amppd/annotation/KTN1_WGS_MB.annovar \
    -remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput

    head -1 $KTN1_DIR/amppd/annotation/KTN1_WGS_MB.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct KTN1_WGS_MB.annovar.hg38_multianno.txt > $KTN1_DIR/amppd/annotation/KTN1.MB.trimmed.annotation.txt

    ##annotate mac1 vcf
    
    plink --bfile KTN1_WGS_AMP_PD_pheno_sex_MB_mac1 --recode 'vcf-fid' --out KTN1_WGS_AMP_PD_pheno_sex_MB_mac1
    bgzip KTN1_WGS_AMP_PD_pheno_sex_MB_mac1.vcf
    tabix -f -p vcf KTN1_WGS_AMP_PD_pheno_sex_MB_mac1.vcf.gz

    table_annovar.pl $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_MB_mac1.vcf.gz $KTN1_DIR/amppd/humandb -buildver hg38 \
    --thread 16 \
    -out $KTN1_DIR/amppd/annotation/KTN1_WGS_MB_mac1.annovar \
    -remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput

    head -1 $KTN1_DIR/amppd/annotation/KTN1_WGS_MB_mac1.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct KTN1_WGS_MB_mac1.annovar.hg38_multianno.txt > $KTN1_DIR/amppd/annotation/KTN1.MB.mac1.trimmed.annotation.txt 
    
    ##annotate mac3 vcf

    plink --bfile KTN1_WGS_AMP_PD_pheno_sex_MB_mac3 --recode 'vcf-fid' --out KTN1_WGS_AMP_PD_pheno_sex_MB_mac3
    bgzip KTN1_WGS_AMP_PD_pheno_sex_MB_mac3.vcf
    tabix -f -p vcf KTN1_WGS_AMP_PD_pheno_sex_MB_mac3.vcf.gz

    table_annovar.pl $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_MB_mac3.vcf.gz $KTN1_DIR/amppd/humandb -buildver hg38 \
    --thread 16 \
    -out $KTN1_DIR/amppd/annotation/KTN1_WGS_MB_mac3.annovar \
    -remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput

    head -1 $KTN1_DIR/amppd/annotation/KTN1_WGS_MB_mac3.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct KTN1_WGS_MB_mac3.annovar.hg38_multianno.txt > $KTN1_DIR/amppd/annotation/KTN1.MB.mac3.trimmed.annotation.txt

    

    # ==================== 4. Burden analysis on all variants via RVTESTS ====================
 

    # mac 1
    ## run burden with no freq, 1%, 3%
    
    ##freq none
    rvtest --noweb --inVcf $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_MB_mac1.vcf.gz \
    --pheno $FILE_DIR/AMP_PD/covs_Mike.txt \
    --covar $FILE_DIR/AMP_PD/covs_Mike.txt \
    --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald \
    --geneFile $FILE_DIR/refFlat_hg38.txt \
    --out AMP_PD_BURDEN.KTN1.MB.mac1


    ##freq = 1%
    rvtest --noweb --inVcf $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_MB_mac1.vcf.gz \
    --pheno $FILE_DIR/AMP_PD/covs_Mike.txt \
    --covar $FILE_DIR/AMP_PD/covs_Mike.txt \
    --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald \
    --geneFile $FILE_DIR/refFlat_hg38.txt \
    --freqUpper 0.01 \
    --out AMP_PD_BURDEN.KTN1.MB.mac1.freq001


    ##freq = 3%
    rvtest --noweb --inVcf $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_MB_mac1.vcf.gz \
    --pheno $FILE_DIR/AMP_PD/covs_Mike.txt \
    --covar $FILE_DIR/AMP_PD/covs_Mike.txt \
    --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald \
    --geneFile $FILE_DIR/refFlat_hg38.txt \
    --freqUpper 0.03 \
    --out AMP_PD_BURDEN.KTN1.MB.mac1.freq003


    # ==================== 5. Single Variant Wald and Score Tests via RVTESTS ====================
    
    # For statistics, need --mac 3: Reason being if there is one allele, then not enough information to do stats, 
    # it is assigned a 0 and results in an NA


    # Tests done on single variants done on ALL variants 
            # Wald = same as logistic in PLINK, better for common variants 
            # score = better for rare variants
    rvtest --noweb --hide-covar --inVcf $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_MB_mac3.vcf.gz \
    --pheno $FILE_DIR/AMP_PD/covs_Mike.txt \
    --pheno-name PHENO \
    --covar $FILE_DIR/AMP_PD/covs_Mike.txt \
    --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --single wald,score --geneFile $FILE_DIR/refFlat_hg38.txt \
    --out AMP_PD_mac3.KTN1.MB.LOGISTIC.allPCs

# Summary Statistics 
### 1. GWAS Summary Statistics (Nalls 2019)
    # ==================== 1. Extract gene from GWAS sum stats + Zoom plots ====================

    # KTN1 positions on hg19 (14: 56025790 - 56168244) 
    # cis-KTN1 (14: 55025790 - 57168244) --> KTN1.MB.txt

    ##making from $KTN1_DIR/annotation/KTN1.trimmed.annotation/txt
    #awk '{ print $1,$2,$3,$8 }' $KTN1_DIR/hardcallsNoNeuroX/annotation/KTN1.MB.trimmed.annotation.txt > $KTN1_DIR/sum_stats/KTN1.MB.txt

    R
    ##run in biowulf
    ## Nalls 2019
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(plyr)
    #data <- fread("Nalls2019_no23andMe.tbl", header = T)
    data <- fread("$FILE_DIR/nallsEtAl2019_excluding23andMe_allVariants.tab", header = T)
    data$CHR <- ldply(strsplit(as.character(data$SNP), split = ":"))[[1]]
    data$CHR <- sub("chr", "", data$CHR)
    data$BP <- ldply(strsplit(as.character(data$SNP), split = ":"))[[2]]
    genesTemp <- read.table("KTN1.MB.txt", sep = " ", header = T)
    colnames(genesTemp) <- c("CHR","START","STOP","GENE")
    genes <- subset(genesTemp, CHR != "X")
    for(i in 1:length(genes$GENE))
    {
      thisGene <- genes$GENE[i]
      thisChr <- genes$CHR[i]
      lower <- genes$START[i]
      upper <- genes$STOP[i]
      output <- subset(data, CHR == thisChr & BP >= lower & BP <= upper)
      file_name <- paste(thisGene,"_risk_variants_mb.tab", sep = "")
      fwrite(output, file = file_name, na = "NA", quote = F, row.names = F, sep = "\t")
    }
    
    ######## Zoom plot ##########
    module load locuszoom
    locuszoom --metal $KTN1_DIR/sum_stats/cisKTN1_variants.Nalls.tab --pvalcol p --markercol SNP --chr 14 --start 55025790 --end 57168244  --pop EUR --build hg19 --source 1000G_March2012 --plotonly &


### 2. Age of Onset GWAS Summary Statistics (Blauwendraat 2019)

    ## Blauwendraat 2019
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(plyr)
    #data <- fread("Blauwendraat2019_no23andMe.tbl", header = T)
    data <- fread("IPDGC_AAO_GWAS_sumstats_april_2018.txt", header = T)
    data$CHR <- ldply(strsplit(as.character(data$MarkerName), split = ":"))[[1]]
    data$BP <- ldply(strsplit(as.character(data$MarkerName), split = ":"))[[2]]
    chr14 <-subset(data, CHR == 'chr14')

    genesTemp <- read.table("KTN1.MB.txt", sep = "", header = F)
    colnames(genesTemp) <- c("CHR","START","STOP","GENE")
    genes <- subset(genesTemp, CHR != "X")

    for(i in 1:length(genes$GENE))
    {
      thisGene <- genes$GENE[i]
      thisChr <- genes$CHR[i]
      lower <- genes$START[i]
      upper <- genes$STOP[i]
      output <- subset(chr14, CHR == thisChr & BP >= lower & BP <= upper)
      file_name <- paste(thisGene,"_aao_variants_mb.tab", sep = "")
      fwrite(output, file = file_name, na = "NA", quote = F, row.names = F, sep = "\t")
    }
    
    ######### Zoom plot ###########
   
    locuszoom --metal $KTN1_DIR/sum_stats/cisKTN1_variants.Blauwendraat.tab --pvalcol P-value --markercol SNP --chr 14 --start 55025790 --end 57168244 --flank 1000kb --pop EUR --build hg19 --source 1000G_March2012 --plotonly &

# KTN1 expression analysis
## 1. eQTL (GTEx v8)
## 2. Expression covariate correction + T-test (AMP-PD)
    all_ktn1_expression_notebook_github.ipynb
