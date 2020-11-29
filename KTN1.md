
# KTN1 GWAS and whole genome analysis
## 1. GWAS (IPDGC genotype data)
#### 1. Subset PLINK Binaries + Convert to VCFs
#### 2. Output Frequency with PLINK
#### 3. Logistic Regression with PLINK
#### 4. Annotate VCFs with Annovar 
#### 5. RVTESTS Burden tests
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

    # ==================== 2. Subset PLINK Binaries ====================
    ## Remove NeuroX + keep males + genotype quality of 15% or less missingness

    plink --bfile $FILE_DIR/HARDCALLS_PD_september_2018_no_cousins \
    --remove-fam $FILE_DIR/NeuroX.fID.txt \
    --chr 14 --geno 0.15 --from-bp 56025790 --to-bp 56168244 \
    --make-bed --out hardcallsNoNeuroX/bin/KTN1.GWAS
        # Done 04.09.2020
        # PLINK output:
            # Total genotyping rate in remaining samples is 0.470074.
            # 929 variants removed due to missing genotype data (--geno).
            # 452 variants and 32338 people pass filters and QC.
            # Among remaining phenotypes, 14671 are cases and 17667 are controls.


    plink --bfile hardcallsNoNeuroX/bin/KTN1.GWAS --recode vcf --out hardcallsNoNeuroX/vcf/KTN1.GWAS
        # Done 04.09.2020
        # PLINK output:
            # Total genotyping rate is 0.988557.
            # 452 variants and 32338 people pass filters and QC.
            # Among remaining phenotypes, 14671 are cases and 17667 are controls.

    cd hardcallsNoNeuroX/vcf
    bgzip KTN1.GWAS.vcf
    tabix -f -p vcf KTN1.GWAS.vcf.gz

    # ==================== 2. Output Frequency ==================== 
    cd $KTN1_DIR/
    # Overall
    plink --bfile $FILE_DIR/HARDCALLS_PD_september_2018_no_cousins \
    --chr 14 --from-bp 56025790 --to-bp 56168244 \
    --remove-fam $FILE_DIR/NeuroX.fID.txt \
    --geno 0.15 --freq \
    --out hardcallsNoNeuroX/freq/KTN1
        # Done 04.09.2020


    # Case-control
    plink --bfile $FILE_DIR/HARDCALLS_PD_september_2018_no_cousins \
    --chr 14 --from-bp 56025790 --to-bp 56168244 \
    --remove-fam $FILE_DIR/NeuroX.fID.txt \
    --geno 0.15 --freq case-control  \
    --out hardcallsNoNeuroX/freq/KTN1
        # Done 04.09.2020

    # ==================== 3. Output Logistic Regression ====================
    plink --bfile $FILE_DIR/HARDCALLS_PD_september_2018_no_cousins \
    --remove-fam $FILE_DIR/NeuroX.fID.txt \
    --chr 14 --from-bp 56025790 --to-bp 56168244 \
    --geno 0.15 \
    --covar $FILE_DIR/IPDGC_all_samples_covariates.tab \
    --covar-name sex,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE \
    --assoc \
    --out hardcallsNoNeuroX/logistic/KTN1
        # Done 04.09.2020

    # ==================== 4. Annotate VCFs with Annovar ====================

    # Downloading databases 
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
    # annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/ 
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 humandb/ 
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a humandb/
    # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_genome humandb/

    table_annovar.pl $KTN1_DIR/hardcallsNoNeuroX/vcf/KTN1.GWAS.vcf.gz $FILE_DIR/annovar/humandb/ -buildver hg19 \
    --thread 16 \
    -out $KTN1_DIR/hardcallsNoNeuroX/annotation/KTN1.annovar \
    -remove -protocol avsnp147,refGene,ensGene,gnomad211_genome \
    -operation f,g,g,f \
    -nastring . \
    -vcfinput

    cd $KTN1_DIR/hardcallsNoNeuroX/annotation
    head -1 KTN1.annovar.hg19_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct KTN1.annovar.hg19_multianno.txt > KTN1.trimmed.annotation.txt
        # Done 04.09.2020

    generating list of coding variants
    awk '$6=="exonic" {print}' KTN1.trimmed.annotation.txt > KTN1.trimmed.annotation.coding.variants.txt
    #awk '{ print $1":"$2 }' KTN1.trimmed.annotation.coding.variants.WGS.txt > KTN1.trimmed.annotation.coding.variants.SNPs.txt
        # Done 04.09.2020
        # Found no coding variants (almost all intronic)

    ##generating list of nonsynonymous variants
    # awk '$9=="nonsynonymous" {print}' KTN1.trimmed.annotation.coding.variants.txt > KTN1.trimmed.annotation.coding.variants.nonsynonymous.txt
    # awk '{ print $1":"$2 }' KTN1.trimmed.annotation.coding.variants.nonsynonymous.txt > KTN1.trimmed.annotation.coding.variants.nonsynonymous.SNPs.txt

    ###check how many coding variants
    ## run burden test only on coding variants if theres enough 


    # ==================== 5. Burden analysis ====================
    # Note that here we do not adjust by AGE since the model does not fit in rvtests

    cd $KTN1_DIR

    rvtest --noweb --inVcf $KTN1_DIR/hardcallsNoNeuroX/vcf/KTN1.GWAS.vcf.gz \
    --pheno $FILE_DIR/IPDGC_all_samples_covariates.vcf.tab \
    --covar $FILE_DIR/IPDGC_all_samples_covariates.vcf.tab \
    --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,DUTCH,FINLAND,GERMANY,MCGILL,MF,NIA,OSLO,PROBAND,PROPARK,SHULMAN,SPAIN3,SPAIN4,TUBI,UK_GWAS,VANCE \
    --burden cmc,zeggini,mb,fp,cmcWald \
    --kernel skat,skato \
    --geneFile $FILE_DIR/refFlat_hg19.txt \
    --freqUpper 0.03 \
    --out hardcallsNoNeuroX/burden/BURDEN.KTN1.freq03


## 2. Whole Genome Analysis (AMP-PD)
#### 1. Extract KTN1 Region in terra
#### 2. Fisher test in PLINK
#### 3. Annotation
#### 4. Burden analysis
#### 5. Single Variant Wald and Score Tests via RVTESTS 
 
 
    #!/bin/bash
    
    ##set up environment
    module load plink samtools rvtests annovar

    ## setting up directories
    cd $KTN1_DIR/

    mkdir AMPPD

    # Initializing some variables 
    COV_NAMES="SEX,AGE,FAMILY_HISTORY,EDUCATION,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,STUDY_BioFIND,STUDY_PDBP,STUDY_PPMI"

    ### Prep file 
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


    # ==================== 1. Extract KTN1 Region  ====================

    ##extract KTN1 region
    # Extract the KTN1 region
        # These are hg38 co-ordinates 
    # no filter
    plink --bfile $KTN1_DIR/amppd/plink/WGS_AMP_PD_pheno_sex \
    --chr 14 --from-bp 55559072 --to-bp 55701526 \
    --make-bed --out $KTN1_DIR1/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex
        # Done 07.10.2020
        # PLINK output:
            # Total genotyping rate is 0.997671.
            # 1513 variants and 3074 people pass filters and QC.
            # Among remaining phenotypes, 1735 are cases and 1100 are controls.

    # --mac 1
    plink --bfile $KTN1_DIR/amppd/plink/WGS_AMP_PD_pheno_sex \
    --chr 14 --from-bp 55559072 --to-bp 55701526 \
    --mac 1 \
    --make-bed --out $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_mac1
        # Done 07.10.2020
        # PLINK output:
            # 1509 variants and 3074 people pass filters and QC.


    # --mac 3 
    plink --bfile $KTN1_DIR/amppd/plink/WGS_AMP_PD_pheno_sex \
    --chr 14 --from-bp 55559072 --to-bp 55701526 \
    --mac 3 \
    --make-bed --out $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_mac3 
            # Done 07.10.20
            # PLINK output:
                # 1089 variants and 3074 people pass filters and QC.

    # wc -l KTN1_WGS_AMP_PD_pheno_sex_mac3.bim
        # 1089 KTN1_WGS_AMP_PD_pheno_sex_mac3.bim


    # ==================== 2. Fisher test ====================
    plink --bfile $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex --fisher --pheno $FILE_DIR/AMP_PD/pheno_Mike.txt --covar $FILE_DIR/AMP_PD/covs_Mike.txt --covar-name SEX, AGE, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 --out KTN1_WGS_AMP_PD --ci 0.95


    # ==================== 3. Annotation ====================

    ## no filters
    plink --bfile KTN1_WGS_AMP_PD_pheno_sex --recode 'vcf-fid' --out KTN1_WGS_AMP_PD_pheno_sex
    bgzip KTN1_WGS_AMP_PD_pheno_sex.vcf
    tabix -f -p vcf KTN1_WGS_AMP_PD_pheno_sex.vcf.gz

    #Now using ANNOVAR

    table_annovar.pl $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex.vcf.gz $KTN1_DIR/amppd/humandb -buildver hg38 \
    --thread 16 \
    -out $KTN1_DIR/amppd/annotation/KTN1_WGS.annovar \
    -remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput

    head -1 $KTN1_DIR/amppd/annotation/KTN1_WGS.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct KTN1_WGS.annovar.hg38_multianno.t1xt > $KTN1_DIR/amppd/annotation/KTN1.trimmed.annotation.txt

    ##annotate mac3 vcf

    plink --bfile KTN1_WGS_AMP_PD_pheno_sex_mac3 --recode 'vcf-fid' --out KTN1_WGS_AMP_PD_pheno_sex_mac3
    bgzip KTN1_WGS_AMP_PD_pheno_sex_mac3.vcf
    tabix -f -p vcf KTN1_WGS_AMP_PD_pheno_sex_mac3.vcf.gz

    table_annovar.pl $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_mac3.vcf.gz $KTN1_DIR/amppd/humandb -buildver hg38 \
    --thread 16 \
    -out $KTN1_DIR/amppd/annotation/KTN1_WGS_mac3.annovar \
    -remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput

    head -1 $KTN1_DIR/amppd/annotation/KTN1_WGS_mac3.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct KTN1_WGS_mac3.annovar.hg38_multianno.txt > $KTN1_DIR/amppd/annotation/KTN1.mac3.trimmed.annotation.txt

    ##annotate mac1 vcf
    plink --bfile KTN1_WGS_AMP_PD_pheno_sex_mac1 --recode 'vcf-fid' --out KTN1_WGS_AMP_PD_pheno_sex_mac1
    bgzip KTN1_WGS_AMP_PD_pheno_sex_mac1.vcf
    tabix -f -p vcf KTN1_WGS_AMP_PD_pheno_sex_mac1.vcf.gz

    table_annovar.pl $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_mac1.vcf.gz $KTN1_DIR/amppd/humandb -buildver hg38 \
    --thread 16 \
    -out $KTN1_DIR/amppd/annotation/KTN1_WGS_mac1.annovar \
    -remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput

    head -1 $KTN1_DIR/amppd/annotation/KTN1_WGS_mac1.annovar.hg38_multianno.txt > header.txt
    colct="$(wc -w header.txt| cut -f1 -d' ')"
    cut -f1-$colct KTN1_WGS_mac1.annovar.hg38_multianno.txt > $KTN1_DIR/amppd/annotation/KTN1.mac1.trimmed.annotation.txt

    # ==================== 4. Burden anaysis ====================
    ##generating list of coding variants
    awk '$6=="exonic" {print}' KTN1.trimmed.annotation.txt > KTN1.trimmed.annotation.coding.variants.WGS.txt
    awk '{print $1" "$2" "$2" "$7}' KTN1.trimmed.annotation.coding.variants.WGS.txt > KTN1.trimmed.annotation.coding.variants.WGS.SNPs.txt
    awk '$9=="nonsynonymous" {print}' KTN1.trimmed.annotation.coding.variants.WGS.txt > KTN1.trimmed.annotation.coding.variants.WGS.nonsynonymous.txt

    plink --bfile $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex --extract range /data/LNG/anni/KTN1/amppd/annotation/KTN1.trimmed.annotation.coding.variants.WGS.SNPs.txt --recode 'vcf-fid' --out /data/LNG/anni/KTN1/amppd/plink/KTN1_CODING_AMP

    bgzip KTN1_CODING_AMP.vcf
    tabix -f -p vcf KTN1_CODING_AMP.vcf.gz

    #MAF < 0.03
    ##all variants
    rvtest --noweb --inVcf $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex.vcf.gz \
    --pheno $FILE_DIR/AMP_PD/covs_Mike.txt \
    --covar $FILE_DIR/AMP_PD/covs_Mike.txt \
    --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald \
    --geneFile $FILE_DIR/refFlat_hg38.txt \
    --freqUpper 0.03 \
    --out AMP_PD_BURDEN.KTN1.maf003

    ##coding variants
    rvtest --noweb --inVcf $KTN1_DIR/amppd/plink/KTN1_CODING_AMP.vcf.gz --pheno $FILE_DIR/AMP_PD/covs_Mike.txt --covar $FILE_DIR/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile $FILE_DIR/refFlat_hg38.txt --freqUpper 0.03 --out AMP_PD_BURDEN.KTN1.maf003_CODING


    ## with mac 1 annotated file
    #generating list of coding variants
    awk '$6=="exonic" {print}' KTN1.mac1.trimmed.annotation.txt > KTN1.mac1.trimmed.annotation.coding.variants.WGS.txt
    awk '{print $1" "$2" "$2" "$7}' KTN1.mac1.trimmed.annotation.coding.variants.WGS.txt > KTN1.mac1.trimmed.annotation.coding.variants.WGS.SNPs.txt
    awk '$9=="nonsynonymous" {print}' KTN1.mac1.trimmed.annotation.coding.variants.WGS.txt > KTN1.mac1.trimmed.annotation.coding.variants.WGS.nonsynonymous.txt

    plink --bfile $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_mac1 --extract range $KTN1_DIR/amppd/annotation/KTN1.mac1.trimmed.annotation.coding.variants.WGS.SNPs.txt --recode 'vcf-fid' --out $KTN1_DIR/amppd/plink/KTN1_mac1_CODING_AMP

    bgzip KTN1_mac1_CODING_AMP.vcf
    tabix -f -p vcf KTN1_mac1_CODING_AMP.vcf.gz

    ## run burden with no freq, 1%, 3%

    ##all variants
    ##freq none
    rvtest --noweb --inVcf $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_mac1.vcf.gz \
    --pheno $FILE_DIR/AMP_PD/covs_Mike.txt \
    --covar $FILE_DIR/AMP_PD/covs_Mike.txt \
    --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald \
    --geneFile $FILE_DIR/refFlat_hg38.txt \
    --out AMP_PD_BURDEN.KTN1.mac1


    ##all variants
    ##freq = 1%
    rvtest --noweb --inVcf $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_mac1.vcf.gz \
    --pheno $FILE_DIR/AMP_PD/covs_Mike.txt \
    --covar $FILE_DIR/AMP_PD/covs_Mike.txt \
    --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald \
    --geneFile $FILE_DIR/refFlat_hg38.txt \
    --freqUpper 0.01 \
    --out AMP_PD_BURDEN.KTN1.mac1.freq001


    ##all variants
    ##freq = 3%
    rvtest --noweb --inVcf $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_mac1.vcf.gz \
    --pheno $FILE_DIR/AMP_PD/covs_Mike.txt \
    --covar $FILE_DIR/AMP_PD/covs_Mike.txt \
    --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald \
    --geneFile $FILE_DIR/refFlat_hg38.txt \
    --freqUpper 0.03 \
    --out AMP_PD_BURDEN.KTN1.mac1.freq003

    ##coding variants
    rvtest --noweb --inVcf $KTN1_DIR/amppd/plink/KTN1_mac1_CODING_AMP.vcf.gz --pheno $FILE_DIR/AMP_PD/covs_Mike.txt --covar $FILE_DIR/AMP_PD/covs_Mike.txt --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --kernel skat,skato --burden cmc,zeggini,mb,fp,cmcWald --geneFile $FILE_DIR/refFlat_hg38.txt --freqUpper 0.01 --out AMP_PD_BURDEN.KTN1.mac1.maf001_CODING

    # ==================== 5. Single Variant Wald and Score Tests via RVTESTS ====================

    ## Generate vcf file for mac 3
    plink --bfile $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_mac3 \
    --mac 3 \
    --recode 'vcf-fid' --out $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_mac3
        # Done 07.10.2020

    bgzip $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_mac3.vcf
    tabix -f -p vcf $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_mac3.vcf.gz
        # Done 07.10.2020


    # For statistics, need --mac 3: Reason being if there is one allele, then not enough information to do stats, 
    # it is assigned a 0 and results in an NA


    # Tests done on single variants done on ALL variants 
            # Wald = same as logistic in PLINK, better for common variants 
            # score = better for rare variants
    rvtest --noweb --hide-covar --inVcf $KTN1_DIR/amppd/plink/KTN1_WGS_AMP_PD_pheno_sex_mac3.vcf.gz \
    --pheno $FILE_DIR/AMP_PD/covs_Mike.txt \
    --pheno-name PHENO \
    --covar $FILE_DIR/AMP_PD/covs_Mike.txt \
    --covar-name SEX,AAO,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --single wald,score --geneFile $FILE_DIR/refFlat_hg38.txt \
    --out AMP_PD_mac3.KTN1.LOGISTIC.allPCs
        # Files generated at AMP_PD_mac3.KTN1.LOGISTIC.allPCs.*.assoc
            # --hide-covar flag to remove each and every individual covariate test from the final output
            # The variant line reported includes all PCs 
        # Done 21.04.2020

# Summary Statistics 
### 1. GWAS Summary Statistics (Nalls 2019)
### 2. Age of Onset GWAS Summary Statistics

# KTN1 expression analysis
## 1. eQTL (GTEx v8)
### Putamen
#### Subset to KTN1 genes
    awk '$2 ~ /KTN1/ {print}' Brain_Putamen_basal_ganglia.v8.egenes.txt > Brain_Putamen_basal_ganglia.v8.egenes.KTN1.txt
    wc -l Brain_Putamen_basal_ganglia.v8.egenes.KTN1.txt 
        ## Output:
        #  2 Brain_Putamen_basal_ganglia.v8.egenes.KTN1.txt

    awk '$2 == "ENSG00000126777.17" {print}' Brain_Putamen_basal_ganglia.v8.signif_variant_gene_pairs.txt > Brain_Putamen_basal_ganglia.v8.signif_variant_gene_pairs_KTN1.txt
    wc -l Brain_Putamen_basal_ganglia.v8.signif_variant_gene_pairs_KTN1.txt
        ## Output:
        #  0 Brain_Putamen_basal_ganglia.v8.signif_variant_gene_pairs_KTN1.txt

### Substantia nigra
#### Subset to KTN1 genes
    awk '$2 ~ /KTN1/ {print}' Brain_Substantia_nigra.v8.egenes.txt > Brain_Substantia_nigra.v8.egenes.KTN1.txt
    wc -l Brain_Substantia_nigra.v8.egenes.KTN1.txt
        ## Output:
        #  2 Brain_Substantia_nigra.v8.egenes.KTN1.txt

    awk '$2 == "ENSG00000126777.17" {print}' Brain_Substantia_nigra.v8.signif_variant_gene_pairs.txt > Brain_Substantia_nigra.v8.signif_variant_gene_pairs_KTN1.txt
    wc -l Brain_Substantia_nigra.v8.signif_variant_gene_pairs_KTN1.txt
        ## Output:
        #  0 Brain_Substantia_nigra.v8.signif_variant_gene_pairs_KTN1.txt

## 2. Expression covariate correction + ANOVA (AMP-PD)
    all_ktn1_expression_notebook_github.ipynb
