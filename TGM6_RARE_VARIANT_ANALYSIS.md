# TGM6 Rare Variant Analysis
## Information

#### DESCRIPTION

 - **Author(s):** Ashley Hall, Kimberley Billingsley 
 -  **Project:**  TGM6 Letter (IPDGC)
 - **Date Last Updated:** 07.07.20
 -  This pipeline takes direction from [here](https://github.com/ipdgc/IPDGC-Trainees/blob/master/LIN28A.md)
#### BACKGROUD


-   **Aim:**   To identify if variants in TGM6 are associated with risk of Parkinson's disease. 
-   **Original Publication:**  https://jin.imrpress.com/article/2020/1757-448X/1757-448X-19-1-51.shtml 

#### WORKFLOW:

 **1.  Extract the TGM6 region from the AMP-PD WGS**
 ****2. Run --model and --assoc in PLINK**
3. Generate VCFs
 4. Annotate variants with ANNOVAR
  5. Burden Analysis with RVtests
   6. Single variant Wald and score test with RVtests**

## 1. EXTRACT TGM6 VARIANTS

    # Create a working directory for the TMG6 analyses 
    # Change into directory
    cd $WORKING_DIR
    
    module load PLINK
    
    # Extract the TGM6 variants
    	# These are hg38 co-ordinates 
    	
    # --mac 1
    plink --bfile $AMP_PLINK \
    --chr 20 --from-bp 2380901 --to-bp 2432753 \
    --mac 1 \
    --make-bed --out $WORKING_DIR/TGM6_WGS_AMP_PD_mac1
   
    # wc -l TGM6_WGS_AMP_PD_mac1.bim
	# 478 TGM6_WGS_AMP_PD_mac1.bim
	
    # --mac 3
    plink --bfile $AMP_PLINK \
    --chr 20 --from-bp 2380901 --to-bp 2432753 \
    --mac 3 \
    --make-bed --out $WORKING_DIR/TGM6_WGS_AMP_PD_mac3
    
	 # wc -l TGM6_WGS_AMP_PD_mac3.bim
	 # 362 TGM6_WGS_AMP_PD_mac3.bim

   ## 2. RUN --MODEL AND --ASSOC IN PLINK
  Run  --model and --assoc to obtain SNPs frequencies and the distribution of alleles of the dataset. 
#### Run with  `--model`

    plink --bfile $WORKING_DIR/TGM6_WGS_AMP_PD_mac3 --model \
    --out TGM6_WGS_AMP_PD_mac3_noCov
  #### Run `--assoc` with no covariates

    plink --bfile $WORKING_DIR/TGM6_WGS_AMP_PD_mac3 --assoc \
    --pheno $PHENO \
    --ci 0.95 \
    --out TGM6_WGS_AMP_PD_mac3_noCov
   ## 3. GENERATE VCFS
   #### Generate a VCF file for  `--mac 1`

    plink --bfile $WORKING_DIR/TGM6_WGS_AMP_PD_mac1 \
    --mac 1 \
    --recode 'vcf-fid' --out $WORKING_DIR/TGM6_WGS_AMP_PD_mac1
#### Generate a VCF file for  `--mac 3`

    plink --bfile $WORKING_DIR/TGM6_WGS_AMP_PD_mac3 \
    --mac 3 \
    --recode 'vcf-fid' --out $WORKING_DIR/TGM6_WGS_AMP_PD_mac3
#### Zip and tabix the VCF for  `--mac 1`

    bgzip $WORKING_DIR/TGM6_WGS_AMP_PD_mac1.vcf
    tabix -f -p vcf $WORKING_DIR/TGM6_WGS_AMP_PD_mac1.vcf.gz
#### Zip and tabix the VCF for  `--mac 3`
    bgzip $WORKING_DIR/TGM6_WGS_AMP_PD_mac3.vcf
    tabix -f -p vcf $WORKING_DIR/TGM6_WGS_AMP_PD_mac3.vcf.gz
## 4. ANNOTATE VARIANTS WITH ANNOVAR


  
    module load ANNOVAR
    table_annovar.pl $WORKING_DIR/TGM6_WGS_AMP_PD_mac1.vcf.gz \
    $ANNOVAR_DATA/hg38 --thread 20 -buildver hg38 \
    -out $WORKING_DIR/TGM6_WGS_AMP_PD_mac1.annovar \
    -arg '-splicing 15',,, \
    -remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
    -operation g,f,f,f -nastring . -vcfinput -polish

## 5. BURDEN ANALYSIS WITH RVTESTS

*For burdens, only need  `--mac 1`

**Extract  TGM6 coding variants** 

    grep "exonic"  $WORKING_DIR/TGM6_WGS_AMP_PD_mac1.annovar.hg38_multianno.txt > $WORKING_DIR/CODING_regions.txt
             cut -f 1,2,3,7 $WORKING_DIR/CODING_regions.txt > CODING.txt
#### Prep the files and generate PLINK binary with only coding variants

    plink --bfile $WORKING_DIR/TGM6_WGS_AMP_PD_mac1 --extract range CODING.txt --make-bed --out $WORKING_DIR/TGM6_exonic_AMP_PD
#### Create the VCF with only coding variants
    plink --bfile $WORKING_DIR/TGM6_exonic_AMP_PD --recode 'vcf-fid' --out $WORKING_DIR/TGM6_exonic_AMP_PD
    bgzip $WORKING_DIR/TGM6_exonic_AMP_PD.vcf
    tabix -f -p vcf $WORKING_DIR/TGM6_exonic_AMP_PD.vcf.gz
### Run RVTests
NO MAF THRESHOLD on ALL and CODING ONLY variants

        ## ALL VARIANTS ##
        module load rvtests
        
        COV_NAMES="SEX,AGE,FAMILY_HISTORY,EDUCATION,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,STUDY_BioFIND,STUDY_PDBP,STUDY_PPMI"
        
        rvtest --noweb --hide-covar --inVcf $WORKING_DIR/TGM6_WGS_AMP_PD_mac1.vcf.gz \
        --pheno covs_Mike.txt \
        --pheno-name PHENO \
        --covar covs_Mike.txt \
        --covar-name $COV_NAMES \
        --kernel skat,skato --geneFile refFlat_hg38.txt \
        --out AMP_PD_BURDEN.TGM6.NOFREQTHRESHOLD
        
	    ## CODING VARIANTS ##
	    rvtest --noweb --hide-covar --inVcf $WORKING_DIR/TGM6_exonic_AMP_PD.vcf.gz \
	    --pheno covs_Mike.txt \
	    --covar covs_Mike.txt \
	    --pheno-name PHENO \
	    --covar-name $COV_NAMES \
	    --kernel skat,skato --geneFile refFlat_hg38.txt \
	    --out AMP_PD_BURDEN.TGM6.NOFREQTHRESHOLD_CODING
#### MAF < 0.01 on ALL and CODING ONLY variants

    rvtest --noweb --hide-covar --inVcf $WORKING_DIR/TGM6_WGS_AMP_PD_mac1.vcf.gz \
    --pheno covs_Mike.txt \
    --pheno-name PHENO \
    --covar covs_Mike.txt \
    --covar-name $COV_NAMES \
    --kernel skat,skato --geneFile refFlat_hg38.txt \
    --freqUpper 0.01 \
    --out AMP_PD_BURDEN.TGM6.maf001
    
    ## CODING VARIANTS ##
    rvtest --noweb --hide-covar --inVcf $WORKING_DIR/TGM6_exonic_AMP_PD.vcf.gz \
    --pheno covs_Mike.txt \
    --covar covs_Mike.txt \
    --pheno-name PHENO \
    --covar-name $COV_NAMES \
    --kernel skat,skato --geneFile refFlat_hg38.txt \
    --freqUpper 0.01 \
    --out AMP_PD_BURDEN.TGM6_CODING.maf001
#### MAF < 0.03 on ALL and CODING ONLY variants
    rvtest --noweb --hide-covar --inVcf $WORKING_DIR/TGM6_WGS_AMP_PD_mac1.vcf.gz \
    --pheno covs_Mike.txt \
    --pheno-name PHENO \
    --covar covs_Mike.txt \
    --covar-name $COV_NAMES \
    --kernel skat,skato --geneFile refFlat_hg38.txt \
    --freqUpper 0.03 \
    --out AMP_PD_BURDEN.TGM6.maf003
    
    ## CODING VARIANTS ##
    rvtest --noweb --hide-covar --inVcf $WORKING_DIR/TGM6_exonic_AMP_PD.vcf.gz \
    --pheno covs_Mike.txt \
    --covar covs_Mike.txt \
    --pheno-name PHENO \
    --covar-name $COV_NAMES \
    --kernel skat,skato --geneFile refFlat_hg38.txt \
    --freqUpper 0.03 \
    --out AMP_PD_BURDEN.TGM6_CODING.maf003

## 6. SINGLE VARIANT WALD AND SCORE TESTS WITH RVTESTS
**For statistics, need --mac 3**: 

#### Single Variant Wald + Score Tests RVTests

rvtest --noweb --hide-covar --inVcf $WORKING_DIR/TGM6_WGS_AMP_PD_mac3.vcf.gz \
--pheno $COV_FILE \
--pheno-name PHENO \
--covar $COV_FILE \
--covar-name $COV_NAMES \
--single wald,score --geneFile refFlat_hg38.txt \
--out AMP_PD_mac3.TGM6.LOGISTIC.allPCs
