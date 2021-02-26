# Candidate PD-risk gene SLC25A46 Evaluation using AMP-PD data

by Jonggeol Jeff Kim and Hui Liu

----

Part of the effort by Hui Liu and Jonggeol Jeffrey Kim to confirm candidate PD-risk gene SLC25A46 based on this paper: https://www.ncbi.nlm.nih.gov/pubmed/32259769

----
Software used:

1. PLINK v2.00a2LM/1.9
2. BCFTools
3. Annovar
4. RVTESTS

----
Steps:
1. [Subset and convert data](#1)
2. [Annotation](#2)
3. [Allele frequency + count](#3)
4. [Burden analysis](#4)
5. [Score test](#5)
6. [Multiple test correction](#6)
7. [Compound Heterozygous Fisher Test](#7)

<a id="1"></a>

# 1. Subset and convert data

The data has already been cleaned (e.g. genotype quality, ancestry etc.) and subsetted into a PLINK binary format.

Subsetting log can be read here:


```bash
# kernel: bash
cat SLC25A46_WGS_AMP_PD.log
```

    PLINK v1.90b6.16 64-bit (19 Feb 2020)
    Options in effect:
      --bfile sample_filtered_geno05_euro
      --chr 5
      --from-bp 110738145
      --make-bed
      --out SLC25A46_WGS_AMP_PD
      --to-bp 110765157
    ...
    243 out of 28195229 variants loaded from .bim file.
    2927 people (0 males, 0 females, 2927 ambiguous) loaded from .fam.
    ...
    Total genotyping rate is 0.999947.
    243 variants and 2927 people pass filters and QC.
    Note: No phenotypes present.
    ...


## Filter data

Filter out participants with sex information missing, as well as variants 0 allele count.

### Pull out columns 1, 2, and 5 (FID, IID, SEX) without header.


```bash
module load plink/2.3-alpha
plink2 --bfile SLC25A46_WGS_AMP_PD \
       --update-sex data/AMP_PD_sex.txt \
       --remove-nosex \
       --mac 1 \
       --make-bed \
       --out data/SLC25A46.filtered
```
    [+] Loading plink  2.3-alpha 

    PLINK v2.00a2.3LM 64-bit Intel (24 Jan 2020)   www.cog-genomics.org/plink/2.0/
    (C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to data/SLC25A46.filtered.log.
    Options in effect:
      --bfile /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/SLC25A46_WGS_AMP_PD
      --mac 1
      --make-bed
      --out data/SLC25A46.filtered
      --remove-nosex
      --update-sex data/AMP_PD_sex.txt
    
    ...
    2927 samples (0 females, 0 males, 2927 ambiguous; 2927 founders) loaded from
    /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/SLC25A46_WGS_AMP_PD.fam.
    243 variants loaded from
    /data/LNG/saraB/AMP_PD/genesfortrainees_2ndpart/SLC25A46_WGS_AMP_PD.bim.
    Note: No phenotype data present.
    --update-sex: 2697 samples updated.
    230 samples removed due to sex filter(s).
    2697 samples (1188 females, 1509 males; 2697 founders) remaining after main
    filters.
    Calculating allele frequencies... done.
    59 variants removed due to allele frequency threshold(s)
    (--maf/--max-maf/--mac/--max-mac).
    184 variants remaining after main filters.
    ...


## Create VCFs


```bash
module load plink
# use plink 1.9 for this
plink --bfile data/SLC25A46.filtered \
      --recode vcf-fid bgz \
      --out data/vcf/SLC25A46.filtered.fid
```

    [-] Unloading plink  2.3-alpha 
    [+] Loading plink  1.9.0-beta4.4  on cn0927 
    
    The following have been reloaded with a version change:
      1) plink/2.3-alpha => plink/1.9.0-beta4.4
    
    PLINK v1.90b4.4 64-bit (21 May 2017)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to data/vcf/SLC25A46.filtered.fid.log.
    Options in effect:
      --bfile data/SLC25A46.filtered
      --out data/vcf/SLC25A46.filtered.fid
      --recode vcf-fid bgz
    
    ...
    184 variants loaded from .bim file.
    2697 people (1509 males, 1188 females) loaded from .fam.
    ...
    Total genotyping rate is 0.99994.
    184 variants and 2697 people pass filters and QC.
    Note: No phenotypes present.
    ...



```bash
cd data/vcf

module load samtools
# tabulate
tabix -f -p vcf SLC25A46.filtered.fid.vcf.gz

cd ../..
```

    [+] Loading samtools 1.10  ... 


<a id="2"></a>
# 2. Determine Allele Frequency + Count


```bash
# bring back plink2
module load plink/2.3-alpha

plink2 --bfile data/SLC25A46.filtered \
       --freq \
       --out result/freq/SLC25A46.filtered.af
plink2 --bfile data/SLC25A46.filtered \
       --geno-counts \
       --out result/freq/SLC25A46.filtered.ac
```

    [+] Loading plink  2.3-alpha 
    PLINK v2.00a2.3LM 64-bit Intel (24 Jan 2020)   www.cog-genomics.org/plink/2.0/
    (C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to result/freq/SLC25A46.filtered.af.log.
    Options in effect:
      --bfile data/SLC25A46.filtered
      --freq
      --out result/freq/SLC25A46.filtered.af
    
    ...
    2697 samples (1188 females, 1509 males; 2697 founders) loaded from
    data/SLC25A46.filtered.fam.
    184 variants loaded from data/SLC25A46.filtered.bim.
    Note: No phenotype data present.
    Calculating allele frequencies... done.
    --freq: Allele frequencies (founders only) written to222232324252526262727282829293030313232333334343535363637383839394040414142424344444545464647474848495050515152525353545455555657575858595960606161626363646465656666676768696970707171727273737475757676777778787979808081828283838484858586868788888989909091919292939494959596969797989899%
    result/freq/SLC25A46.filtered.af.afreq .
    ...
    PLINK v2.00a2.3LM 64-bit Intel (24 Jan 2020)   www.cog-genomics.org/plink/2.0/
    (C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to result/freq/SLC25A46.filtered.ac.log.
    Options in effect:
      --bfile data/SLC25A46.filtered
      --geno-counts
      --out result/freq/SLC25A46.filtered.ac
    
    ...
    2697 samples (1188 females, 1509 males; 2697 founders) loaded from
    data/SLC25A46.filtered.fam.
    184 variants loaded from data/SLC25A46.filtered.bim.
    Note: No phenotype data present.
    Calculating allele frequencies... done.
    --geno-counts: Genotype counts written to1717181919202021212222232324252526262727282829293030313232333334343535363637383839394040414142424344444545464647474848495050515152525353545455555657575858595960606161626363646465656666676768696970707171727273737475757676777778787979808081828283838484858586868788888989909091919292939494959596969797989899%
    result/freq/SLC25A46.filtered.ac.gcount .
    ...


<a id="3"></a>
# 3. Annotation


```bash
table_annovar.pl data/vcf/SLC25A46.filtered.fid.vcf.gz \
                 $ANNOVAR_DATA/hg38 --thread 20 -buildver hg38 \
                 -out result/anno/SLC25A46.filtered.fid.annovar \
                 -arg '-splicing 15',,, \
                 -remove -protocol refGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
                 -operation g,f,f,f -nastring . -vcfinput -polish
```


## Trim results

```bash
# remove unnecessary VCF file
rm result/anno/SLC25A46.filtered.fid.annovar.hg38_multianno.vcf
```


```bash
# truncate unnecessay genotype information
cd result/anno

head -1 SLC25A46.filtered.fid.annovar.hg38_multianno.txt > header.txt
colct="$(wc -w header.txt| cut -f1 -d' ')"
cut -f1-$colct SLC25A46.filtered.fid.annovar.hg38_multianno.txt > SLC25A46.filtered.fid.annovar.hg38_multianno.trimmed.txt
rm header.txt

cd ../..
```

<a id="4"></a>
# 4. Gene-Level Burden Analyses


```python
# makes rvtests more manageable
grep SLC25A46 /data/LNG/KIM/burden_resources/refFlat_hg38.txt > data/SLC25A46.hg38.genefile
```


```python
module load rvtests

rvtest --inVcf data/vcf/SLC25A46.filtered.fid.vcf.gz \
       --pheno data/covariateInfo_AMPPD.txt \
       --covar data/covariateInfo_AMPPD.txt \
       --covar-name SEX,AGE,EDUCATION,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,STUDY_BioFIND,STUDY_PDBP,STUDY_PPMI \
       --kernel skat,skato \
       --burden cmc,zeggini,mb,fp,cmcWald \
       --geneFile data/SLC25A46.hg38.genefile \
       --freqUpper 0.03 \
       --out result/burden/SLC25A46.filtered.burden.maf03


rvtest --inVcf data/vcf/SLC25A46.filtered.fid.vcf.gz \
       --pheno data/covariateInfo_AMPPD.txt \
       --covar data/covariateInfo_AMPPD.txt \
       --covar-name SEX,AGE,EDUCATION,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,STUDY_BioFIND,STUDY_PDBP,STUDY_PPMI \
       --kernel skat,skato \
       --burden cmc,zeggini,mb,fp,cmcWald \
       --geneFile data/SLC25A46.hg38.genefile \
       --freqUpper 0.01 \
       --out result/burden/SLC25A46.filtered.burden.maf01
```

    [+] Loading rvtests  2.1.0  on cn0929 
    Thank you for using rvtests (version: 20190205, git: c86e589efef15382603300dc7f4c3394c82d69b8)
      For documentations, refer to http://zhanxw.github.io/rvtests/
      For questions and comments, plase send to Xiaowei Zhan <zhanxw@umich.edu>
      For bugs and feature requests, please submit at: https://github.com/zhanxw/rvtests/issues
    ...
    [INFO]	Program version: 20190205
    [INFO]	Analysis started at: Thu Sep  3 18:55:49 2020
    [INFO]	Loaded [ 2697 ] samples from genotype files
    [INFO]	Loaded [ 2697 ] sample phenotypes
    [INFO]	Begin to read covariate file
    [INFO]	Loaded 1509 male, 1188 female and 0 sex-unknown samples from data/covariateInfo_AMPPD.txt
    [INFO]	Loaded 1647 cases, 1050 controls, and 0 missing phenotypes
    [WARN]	-- Enabling binary phenotype mode -- 
    [INFO]	Analysis begins with [ 2697 ] samples...
    [INFO]	MadsonBrowning test significance will be evaluated using 10000 permutations at alpha = 0.05
    [INFO]	SKAT test significance will be evaluated using 10000 permutations at alpha = 0.05 weight = Beta[beta1 = 1.00, beta2 = 25.00]
    [INFO]	SKAT-O test significance will be evaluated using weight = Beta[beta1 = 1.00, beta2 = 25.00]
    [INFO]	Loaded [ 1 ] genes.
    [INFO]	Impute missing genotype to mean (by default)
    [INFO]	Set upper minor allele frequency limit to 0.03
    [INFO]	Analysis started
    [INFO]	Analyzed [ 132 ] variants from [ 1 ] genes/regions
    [INFO]	Analysis ends at: Thu Sep  3 18:55:53 2020
    [INFO]	Analysis took 4 seconds
    RVTESTS finished successfully
    ...
    [INFO]	Program version: 20190205
    [INFO]	Analysis started at: Thu Sep  3 18:55:53 2020
    [INFO]	Loaded [ 2697 ] samples from genotype files
    [INFO]	Loaded [ 2697 ] sample phenotypes
    [INFO]	Begin to read covariate file
    [INFO]	Loaded 1509 male, 1188 female and 0 sex-unknown samples from data/covariateInfo_AMPPD.txt
    [INFO]	Loaded 1647 cases, 1050 controls, and 0 missing phenotypes
    [WARN]	-- Enabling binary phenotype mode -- 
    [INFO]	Analysis begins with [ 2697 ] samples...
    [INFO]	MadsonBrowning test significance will be evaluated using 10000 permutations at alpha = 0.05
    [INFO]	SKAT test significance will be evaluated using 10000 permutations at alpha = 0.05 weight = Beta[beta1 = 1.00, beta2 = 25.00]
    [INFO]	SKAT-O test significance will be evaluated using weight = Beta[beta1 = 1.00, beta2 = 25.00]
    [INFO]	Loaded [ 1 ] genes.
    [INFO]	Impute missing genotype to mean (by default)
    [INFO]	Set upper minor allele frequency limit to 0.01
    [INFO]	Analysis started
    [INFO]	Analyzed [ 115 ] variants from [ 1 ] genes/regions
    [INFO]	Analysis ends at: Thu Sep  3 18:55:56 2020
    [INFO]	Analysis took 3 seconds
    RVTESTS finished successfully


<a id="5"></a>
# 5. Score Test

```python
rvtest --inVcf data/vcf/SLC25A46.filtered.fid.vcf.gz \
       --pheno data/covariateInfo_AMPPD.txt \
       --covar data/covariateInfo_AMPPD.txt \
       --covar-name SEX,AGE,EDUCATION,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,STUDY_BioFIND,STUDY_PDBP,STUDY_PPMI \
       --single score \
       --geneFile data/SLC25A46.hg38.genefile \
       --out result/score/SLC25A46.filtered.score
```

    Thank you for using rvtests (version: 20190205, git: c86e589efef15382603300dc7f4c3394c82d69b8)
      For documentations, refer to http://zhanxw.github.io/rvtests/
      For questions and comments, plase send to Xiaowei Zhan <zhanxw@umich.edu>
      For bugs and feature requests, please submit at: https://github.com/zhanxw/rvtests/issues
    ...
    [INFO]	Program version: 20190205
    [INFO]	Analysis started at: Fri Sep  4 17:27:10 2020
    [INFO]	Loaded [ 2697 ] samples from genotype files
    [INFO]	Loaded [ 2697 ] sample phenotypes
    [INFO]	Begin to read covariate file
    [INFO]	Loaded 1509 male, 1188 female and 0 sex-unknown samples from data/covariateInfo_AMPPD.txt
    [INFO]	Loaded 1647 cases, 1050 controls, and 0 missing phenotypes
    [WARN]	-- Enabling binary phenotype mode -- 
    [INFO]	Analysis begins with [ 2697 ] samples...
    [INFO]	Loaded [ 1 ] genes.
    [INFO]	Impute missing genotype to mean (by default)
    [INFO]	Analysis started
    [INFO]	Analyzed [ 718 ] variants from [ 1 ] genes/regions
    [INFO]	Analysis ends at: Fri Sep  4 17:27:14 2020
    [INFO]	Analysis took 4 seconds
    RVTESTS finished successfully


<a id="6"></a>
# 6. Add multiple-test corrections to score results

We will use False Discovery Rate (FDR)


```python
# kernel: python 3.7
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.stats as sm
from patsy import dmatrices

os.chdir("/data/LNG/KIM/burdenAnalysis/SLC25A46-replication/")
```


```python
score = pd.read_csv("result/score/SLC25A46.filtered.score.SingleScore.assoc", sep = "\t")
```


```python
score
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Gene</th>
      <th>CHROM</th>
      <th>POS</th>
      <th>REF</th>
      <th>ALT</th>
      <th>N_INFORMATIVE</th>
      <th>AF</th>
      <th>U</th>
      <th>V</th>
      <th>STAT</th>
      <th>DIRECTION</th>
      <th>EFFECT</th>
      <th>SE</th>
      <th>PVALUE</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110738249</td>
      <td>T</td>
      <td>C</td>
      <td>2697</td>
      <td>0.000927</td>
      <td>1.069550</td>
      <td>0.908378</td>
      <td>1.259320</td>
      <td>+</td>
      <td>1.177430</td>
      <td>1.049220</td>
      <td>0.261780</td>
    </tr>
    <tr>
      <th>1</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110738407</td>
      <td>C</td>
      <td>T</td>
      <td>2697</td>
      <td>0.000371</td>
      <td>-0.045776</td>
      <td>0.388462</td>
      <td>0.005394</td>
      <td>-</td>
      <td>-0.117838</td>
      <td>1.604450</td>
      <td>0.941452</td>
    </tr>
    <tr>
      <th>2</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110738538</td>
      <td>A</td>
      <td>C</td>
      <td>2697</td>
      <td>0.002039</td>
      <td>0.600238</td>
      <td>2.072550</td>
      <td>0.173837</td>
      <td>+</td>
      <td>0.289614</td>
      <td>0.694621</td>
      <td>0.676724</td>
    </tr>
    <tr>
      <th>3</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110738643</td>
      <td>G</td>
      <td>A</td>
      <td>2697</td>
      <td>0.000185</td>
      <td>-0.109522</td>
      <td>0.097021</td>
      <td>0.123633</td>
      <td>-</td>
      <td>-1.128850</td>
      <td>3.210460</td>
      <td>0.725127</td>
    </tr>
    <tr>
      <th>4</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110738696</td>
      <td>T</td>
      <td>C</td>
      <td>2697</td>
      <td>0.021505</td>
      <td>-0.053915</td>
      <td>20.393400</td>
      <td>0.000143</td>
      <td>-</td>
      <td>-0.002644</td>
      <td>0.221440</td>
      <td>0.990474</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>713</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110764313</td>
      <td>T</td>
      <td>C</td>
      <td>2697</td>
      <td>0.001112</td>
      <td>-1.321120</td>
      <td>1.396020</td>
      <td>1.250240</td>
      <td>-</td>
      <td>-0.946349</td>
      <td>0.846359</td>
      <td>0.263506</td>
    </tr>
    <tr>
      <th>714</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110764363</td>
      <td>C</td>
      <td>T</td>
      <td>2697</td>
      <td>0.117723</td>
      <td>-3.186630</td>
      <td>119.149000</td>
      <td>0.085226</td>
      <td>-</td>
      <td>-0.026745</td>
      <td>0.091612</td>
      <td>0.770336</td>
    </tr>
    <tr>
      <th>715</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110764618</td>
      <td>G</td>
      <td>A</td>
      <td>2697</td>
      <td>0.099926</td>
      <td>-12.843500</td>
      <td>81.840200</td>
      <td>2.015590</td>
      <td>-</td>
      <td>-0.156934</td>
      <td>0.110539</td>
      <td>0.155691</td>
    </tr>
    <tr>
      <th>716</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110764674</td>
      <td>G</td>
      <td>GA</td>
      <td>2697</td>
      <td>0.000927</td>
      <td>0.547030</td>
      <td>0.464448</td>
      <td>0.644296</td>
      <td>+</td>
      <td>1.177810</td>
      <td>1.467340</td>
      <td>0.422159</td>
    </tr>
    <tr>
      <th>717</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110764969</td>
      <td>G</td>
      <td>T</td>
      <td>2697</td>
      <td>0.243789</td>
      <td>-14.176400</td>
      <td>139.849000</td>
      <td>1.437040</td>
      <td>-</td>
      <td>-0.101369</td>
      <td>0.084561</td>
      <td>0.230618</td>
    </tr>
  </tbody>
</table>
<p>718 rows × 14 columns</p>
</div>




```python
score = score.drop_duplicates()
```


```python
score
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Gene</th>
      <th>CHROM</th>
      <th>POS</th>
      <th>REF</th>
      <th>ALT</th>
      <th>N_INFORMATIVE</th>
      <th>AF</th>
      <th>U</th>
      <th>V</th>
      <th>STAT</th>
      <th>DIRECTION</th>
      <th>EFFECT</th>
      <th>SE</th>
      <th>PVALUE</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110738249</td>
      <td>T</td>
      <td>C</td>
      <td>2697</td>
      <td>0.000927</td>
      <td>1.069550</td>
      <td>0.908378</td>
      <td>1.259320</td>
      <td>+</td>
      <td>1.177430</td>
      <td>1.049220</td>
      <td>0.261780</td>
    </tr>
    <tr>
      <th>1</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110738407</td>
      <td>C</td>
      <td>T</td>
      <td>2697</td>
      <td>0.000371</td>
      <td>-0.045776</td>
      <td>0.388462</td>
      <td>0.005394</td>
      <td>-</td>
      <td>-0.117838</td>
      <td>1.604450</td>
      <td>0.941452</td>
    </tr>
    <tr>
      <th>2</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110738538</td>
      <td>A</td>
      <td>C</td>
      <td>2697</td>
      <td>0.002039</td>
      <td>0.600238</td>
      <td>2.072550</td>
      <td>0.173837</td>
      <td>+</td>
      <td>0.289614</td>
      <td>0.694621</td>
      <td>0.676724</td>
    </tr>
    <tr>
      <th>3</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110738643</td>
      <td>G</td>
      <td>A</td>
      <td>2697</td>
      <td>0.000185</td>
      <td>-0.109522</td>
      <td>0.097021</td>
      <td>0.123633</td>
      <td>-</td>
      <td>-1.128850</td>
      <td>3.210460</td>
      <td>0.725127</td>
    </tr>
    <tr>
      <th>4</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110738696</td>
      <td>T</td>
      <td>C</td>
      <td>2697</td>
      <td>0.021505</td>
      <td>-0.053915</td>
      <td>20.393400</td>
      <td>0.000143</td>
      <td>-</td>
      <td>-0.002644</td>
      <td>0.221440</td>
      <td>0.990474</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>179</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110764313</td>
      <td>T</td>
      <td>C</td>
      <td>2697</td>
      <td>0.001112</td>
      <td>-1.321120</td>
      <td>1.396020</td>
      <td>1.250240</td>
      <td>-</td>
      <td>-0.946349</td>
      <td>0.846359</td>
      <td>0.263506</td>
    </tr>
    <tr>
      <th>180</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110764363</td>
      <td>C</td>
      <td>T</td>
      <td>2697</td>
      <td>0.117723</td>
      <td>-3.186630</td>
      <td>119.149000</td>
      <td>0.085226</td>
      <td>-</td>
      <td>-0.026745</td>
      <td>0.091612</td>
      <td>0.770336</td>
    </tr>
    <tr>
      <th>181</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110764618</td>
      <td>G</td>
      <td>A</td>
      <td>2697</td>
      <td>0.099926</td>
      <td>-12.843500</td>
      <td>81.840200</td>
      <td>2.015590</td>
      <td>-</td>
      <td>-0.156934</td>
      <td>0.110539</td>
      <td>0.155691</td>
    </tr>
    <tr>
      <th>182</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110764674</td>
      <td>G</td>
      <td>GA</td>
      <td>2697</td>
      <td>0.000927</td>
      <td>0.547030</td>
      <td>0.464448</td>
      <td>0.644296</td>
      <td>+</td>
      <td>1.177810</td>
      <td>1.467340</td>
      <td>0.422159</td>
    </tr>
    <tr>
      <th>183</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110764969</td>
      <td>G</td>
      <td>T</td>
      <td>2697</td>
      <td>0.243789</td>
      <td>-14.176400</td>
      <td>139.849000</td>
      <td>1.437040</td>
      <td>-</td>
      <td>-0.101369</td>
      <td>0.084561</td>
      <td>0.230618</td>
    </tr>
  </tbody>
</table>
<p>184 rows × 14 columns</p>
</div>




```python
score['FDR-Corrected'] = sm.multitest.multipletests(score['PVALUE'])[1]
# alpha = 0.05
```

    /usr/local/Anaconda/envs/py3.7/lib/python3.7/site-packages/ipykernel_launcher.py:1: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      """Entry point for launching an IPython kernel.



```python
score['FDR-Corrected']
```




    0      1.0
    1      1.0
    2      1.0
    3      1.0
    4      1.0
          ... 
    179    1.0
    180    1.0
    181    1.0
    182    1.0
    183    1.0
    Name: FDR-Corrected, Length: 184, dtype: float64




```python
score.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Gene</th>
      <th>CHROM</th>
      <th>POS</th>
      <th>REF</th>
      <th>ALT</th>
      <th>N_INFORMATIVE</th>
      <th>AF</th>
      <th>U</th>
      <th>V</th>
      <th>STAT</th>
      <th>DIRECTION</th>
      <th>EFFECT</th>
      <th>SE</th>
      <th>PVALUE</th>
      <th>FDR-Corrected</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110738249</td>
      <td>T</td>
      <td>C</td>
      <td>2697</td>
      <td>0.000927</td>
      <td>1.069550</td>
      <td>0.908378</td>
      <td>1.259320</td>
      <td>+</td>
      <td>1.177430</td>
      <td>1.049220</td>
      <td>0.261780</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110738407</td>
      <td>C</td>
      <td>T</td>
      <td>2697</td>
      <td>0.000371</td>
      <td>-0.045776</td>
      <td>0.388462</td>
      <td>0.005394</td>
      <td>-</td>
      <td>-0.117838</td>
      <td>1.604450</td>
      <td>0.941452</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110738538</td>
      <td>A</td>
      <td>C</td>
      <td>2697</td>
      <td>0.002039</td>
      <td>0.600238</td>
      <td>2.072550</td>
      <td>0.173837</td>
      <td>+</td>
      <td>0.289614</td>
      <td>0.694621</td>
      <td>0.676724</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110738643</td>
      <td>G</td>
      <td>A</td>
      <td>2697</td>
      <td>0.000185</td>
      <td>-0.109522</td>
      <td>0.097021</td>
      <td>0.123633</td>
      <td>-</td>
      <td>-1.128850</td>
      <td>3.210460</td>
      <td>0.725127</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110738696</td>
      <td>T</td>
      <td>C</td>
      <td>2697</td>
      <td>0.021505</td>
      <td>-0.053915</td>
      <td>20.393400</td>
      <td>0.000143</td>
      <td>-</td>
      <td>-0.002644</td>
      <td>0.221440</td>
      <td>0.990474</td>
      <td>1.0</td>
    </tr>
  </tbody>
</table>
</div>




```python
score.to_csv("result/score/SLC25A46.filtered.score.SingleScore.FDR.dup_cleaned.assoc", sep="\t", index=False)
```

Let's also filter out the non-synonymous variants of interest. They are:

* 5:110739267
* 5:110739354
* 5:110746300
* 5:110756712
* 5:110761217
* 5:110761271
* 5:110761292
* 5:110761301
* 5:110761662


```python
anno = pd.read_csv("result/anno/SLC25A46.filtered.fid.annovar.hg38_multianno.trimmed.txt", sep = "\t")
anno.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Chr</th>
      <th>Start</th>
      <th>End</th>
      <th>Ref</th>
      <th>Alt</th>
      <th>Func.refGene</th>
      <th>Gene.refGene</th>
      <th>GeneDetail.refGene</th>
      <th>ExonicFunc.refGene</th>
      <th>AAChange.refGene</th>
      <th>...</th>
      <th>non_topmed_AF_popmax</th>
      <th>non_neuro_AF_popmax</th>
      <th>non_cancer_AF_popmax</th>
      <th>controls_AF_popmax</th>
      <th>CLNALLELEID</th>
      <th>CLNDN</th>
      <th>CLNDISDB</th>
      <th>CLNREVSTAT</th>
      <th>CLNSIG</th>
      <th>Otherinfo</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>5</td>
      <td>110738249</td>
      <td>110738249</td>
      <td>T</td>
      <td>C</td>
      <td>splicing</td>
      <td>SLC25A46</td>
      <td>NM_001303250:exon1:c.10+2T&gt;C</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>0.0006</td>
      <td>0.0005</td>
      <td>.</td>
      <td>0.0008</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>0.000927</td>
    </tr>
    <tr>
      <th>1</th>
      <td>5</td>
      <td>110738407</td>
      <td>110738407</td>
      <td>C</td>
      <td>T</td>
      <td>intronic</td>
      <td>SLC25A46</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>0.0267</td>
      <td>0.0269</td>
      <td>.</td>
      <td>0.0292</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>0.000371</td>
    </tr>
    <tr>
      <th>2</th>
      <td>5</td>
      <td>110738538</td>
      <td>110738538</td>
      <td>A</td>
      <td>C</td>
      <td>intronic</td>
      <td>SLC25A46</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>0.0037</td>
      <td>0.0036</td>
      <td>.</td>
      <td>0.0081</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>0.002039</td>
    </tr>
    <tr>
      <th>3</th>
      <td>5</td>
      <td>110738643</td>
      <td>110738643</td>
      <td>G</td>
      <td>A</td>
      <td>intronic</td>
      <td>SLC25A46</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>0.000185</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>110738696</td>
      <td>110738696</td>
      <td>T</td>
      <td>C</td>
      <td>intronic</td>
      <td>SLC25A46</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>0.0302</td>
      <td>0.0264</td>
      <td>.</td>
      <td>0.0285</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>0.021510</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 58 columns</p>
</div>




```python
is_nonsyn = anno['ExonicFunc.refGene']=="nonsynonymous SNV"
anno_nonsyn = anno[is_nonsyn]
anno_nonsyn
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Chr</th>
      <th>Start</th>
      <th>End</th>
      <th>Ref</th>
      <th>Alt</th>
      <th>Func.refGene</th>
      <th>Gene.refGene</th>
      <th>GeneDetail.refGene</th>
      <th>ExonicFunc.refGene</th>
      <th>AAChange.refGene</th>
      <th>...</th>
      <th>non_topmed_AF_popmax</th>
      <th>non_neuro_AF_popmax</th>
      <th>non_cancer_AF_popmax</th>
      <th>controls_AF_popmax</th>
      <th>CLNALLELEID</th>
      <th>CLNDN</th>
      <th>CLNDISDB</th>
      <th>CLNREVSTAT</th>
      <th>CLNSIG</th>
      <th>Otherinfo</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>10</th>
      <td>5</td>
      <td>110739267</td>
      <td>110739267</td>
      <td>C</td>
      <td>A</td>
      <td>exonic</td>
      <td>SLC25A46</td>
      <td>.</td>
      <td>nonsynonymous SNV</td>
      <td>SLC25A46:NM_001303249:exon1:c.C148A:p.P50T,SLC...</td>
      <td>...</td>
      <td>0.0012</td>
      <td>0.0018</td>
      <td>.</td>
      <td>0.0004</td>
      <td>520610</td>
      <td>NEUROPATHY,_HEREDITARY_MOTOR_AND_SENSORY,_TYPE...</td>
      <td>MedGen:C4225302,OMIM:616505</td>
      <td>criteria_provided,_single_submitter</td>
      <td>Uncertain_significance</td>
      <td>0.000742</td>
    </tr>
    <tr>
      <th>11</th>
      <td>5</td>
      <td>110739354</td>
      <td>110739354</td>
      <td>G</td>
      <td>A</td>
      <td>exonic</td>
      <td>SLC25A46</td>
      <td>.</td>
      <td>nonsynonymous SNV</td>
      <td>SLC25A46:NM_001303249:exon1:c.G235A:p.E79K,SLC...</td>
      <td>...</td>
      <td>0.0014</td>
      <td>0.0015</td>
      <td>.</td>
      <td>0.0011</td>
      <td>520315</td>
      <td>NEUROPATHY,_HEREDITARY_MOTOR_AND_SENSORY,_TYPE...</td>
      <td>MedGen:C4225302,OMIM:616505</td>
      <td>criteria_provided,_single_submitter</td>
      <td>Likely_benign</td>
      <td>0.002225</td>
    </tr>
    <tr>
      <th>61</th>
      <td>5</td>
      <td>110746300</td>
      <td>110746300</td>
      <td>C</td>
      <td>A</td>
      <td>exonic</td>
      <td>SLC25A46</td>
      <td>.</td>
      <td>nonsynonymous SNV</td>
      <td>SLC25A46:NM_001303249:exon4:c.C416A:p.T139N,SL...</td>
      <td>...</td>
      <td>0.0011</td>
      <td>0.0010</td>
      <td>.</td>
      <td>0.0011</td>
      <td>453904</td>
      <td>NEUROPATHY,_HEREDITARY_MOTOR_AND_SENSORY,_TYPE...</td>
      <td>MedGen:C4225302,OMIM:616505</td>
      <td>criteria_provided,_single_submitter</td>
      <td>Likely_benign</td>
      <td>0.001112</td>
    </tr>
    <tr>
      <th>129</th>
      <td>5</td>
      <td>110756712</td>
      <td>110756712</td>
      <td>G</td>
      <td>A</td>
      <td>exonic</td>
      <td>SLC25A46</td>
      <td>.</td>
      <td>nonsynonymous SNV</td>
      <td>SLC25A46:NM_001303249:exon7:c.G631A:p.V211M,SL...</td>
      <td>...</td>
      <td>0.0035</td>
      <td>0.0030</td>
      <td>.</td>
      <td>0.0035</td>
      <td>453906</td>
      <td>NEUROPATHY,_HEREDITARY_MOTOR_AND_SENSORY,_TYPE...</td>
      <td>MedGen:C4225302,OMIM:616505</td>
      <td>criteria_provided,_single_submitter</td>
      <td>Likely_benign</td>
      <td>0.003337</td>
    </tr>
    <tr>
      <th>156</th>
      <td>5</td>
      <td>110761217</td>
      <td>110761217</td>
      <td>G</td>
      <td>A</td>
      <td>exonic</td>
      <td>SLC25A46</td>
      <td>.</td>
      <td>nonsynonymous SNV</td>
      <td>SLC25A46:NM_001303250:exon8:c.G419A:p.R140Q,SL...</td>
      <td>...</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>0.000185</td>
    </tr>
    <tr>
      <th>158</th>
      <td>5</td>
      <td>110761271</td>
      <td>110761271</td>
      <td>G</td>
      <td>A</td>
      <td>exonic</td>
      <td>SLC25A46</td>
      <td>.</td>
      <td>nonsynonymous SNV</td>
      <td>SLC25A46:NM_001303250:exon8:c.G473A:p.G158D,SL...</td>
      <td>...</td>
      <td>0.0004</td>
      <td>0.0002</td>
      <td>.</td>
      <td>0.0004</td>
      <td>359143</td>
      <td>NEUROPATHY,_HEREDITARY_MOTOR_AND_SENSORY,_TYPE...</td>
      <td>MedGen:C4225302,OMIM:616505</td>
      <td>no_assertion_criteria_provided</td>
      <td>Pathogenic</td>
      <td>0.000185</td>
    </tr>
    <tr>
      <th>159</th>
      <td>5</td>
      <td>110761292</td>
      <td>110761292</td>
      <td>A</td>
      <td>G</td>
      <td>exonic</td>
      <td>SLC25A46</td>
      <td>.</td>
      <td>nonsynonymous SNV</td>
      <td>SLC25A46:NM_001303250:exon8:c.A494G:p.K165R,SL...</td>
      <td>...</td>
      <td>0.0041</td>
      <td>0.0037</td>
      <td>.</td>
      <td>0.0053</td>
      <td>454455</td>
      <td>NEUROPATHY,_HEREDITARY_MOTOR_AND_SENSORY,_TYPE...</td>
      <td>MedGen:C4225302,OMIM:616505</td>
      <td>criteria_provided,_single_submitter</td>
      <td>Likely_benign</td>
      <td>0.002225</td>
    </tr>
    <tr>
      <th>160</th>
      <td>5</td>
      <td>110761301</td>
      <td>110761301</td>
      <td>T</td>
      <td>G</td>
      <td>exonic</td>
      <td>SLC25A46</td>
      <td>.</td>
      <td>nonsynonymous SNV</td>
      <td>SLC25A46:NM_001303250:exon8:c.T503G:p.L168R,SL...</td>
      <td>...</td>
      <td>.</td>
      <td>7.348e-05</td>
      <td>.</td>
      <td>.</td>
      <td>454461</td>
      <td>NEUROPATHY,_HEREDITARY_MOTOR_AND_SENSORY,_TYPE...</td>
      <td>MedGen:C4225302,OMIM:616505</td>
      <td>criteria_provided,_single_submitter</td>
      <td>Uncertain_significance</td>
      <td>0.000185</td>
    </tr>
    <tr>
      <th>161</th>
      <td>5</td>
      <td>110761662</td>
      <td>110761662</td>
      <td>G</td>
      <td>T</td>
      <td>exonic</td>
      <td>SLC25A46</td>
      <td>.</td>
      <td>nonsynonymous SNV</td>
      <td>SLC25A46:NM_001303249:exon8:c.G894T:p.E298D,SL...</td>
      <td>...</td>
      <td>0.1177</td>
      <td>0.1124</td>
      <td>.</td>
      <td>0.1157</td>
      <td>454465</td>
      <td>NEUROPATHY,_HEREDITARY_MOTOR_AND_SENSORY,_TYPE...</td>
      <td>MedGen:C4225302,OMIM:616505</td>
      <td>criteria_provided,_single_submitter</td>
      <td>Benign</td>
      <td>0.000371</td>
    </tr>
  </tbody>
</table>
<p>9 rows × 58 columns</p>
</div>




```python
nonsym_POS = anno_nonsyn.rename(columns={"Start": "POS"})[['POS']]
score_nonsyn = pd.merge(nonsym_POS, score, on="POS")
score_nonsyn
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>POS</th>
      <th>Gene</th>
      <th>CHROM</th>
      <th>REF</th>
      <th>ALT</th>
      <th>N_INFORMATIVE</th>
      <th>AF</th>
      <th>U</th>
      <th>V</th>
      <th>STAT</th>
      <th>DIRECTION</th>
      <th>EFFECT</th>
      <th>SE</th>
      <th>PVALUE</th>
      <th>FDR-Corrected</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>110739267</td>
      <td>SLC25A46</td>
      <td>5</td>
      <td>C</td>
      <td>A</td>
      <td>2697</td>
      <td>0.000742</td>
      <td>0.823836</td>
      <td>0.599113</td>
      <td>1.132850</td>
      <td>+</td>
      <td>1.375090</td>
      <td>1.291950</td>
      <td>0.287168</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>110739354</td>
      <td>SLC25A46</td>
      <td>5</td>
      <td>G</td>
      <td>A</td>
      <td>2697</td>
      <td>0.002225</td>
      <td>-0.172238</td>
      <td>1.856260</td>
      <td>0.015982</td>
      <td>-</td>
      <td>-0.092788</td>
      <td>0.733974</td>
      <td>0.899401</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>110746300</td>
      <td>SLC25A46</td>
      <td>5</td>
      <td>C</td>
      <td>A</td>
      <td>2697</td>
      <td>0.001112</td>
      <td>0.048595</td>
      <td>0.935050</td>
      <td>0.002526</td>
      <td>+</td>
      <td>0.051971</td>
      <td>1.034150</td>
      <td>0.959920</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>110756712</td>
      <td>SLC25A46</td>
      <td>5</td>
      <td>G</td>
      <td>A</td>
      <td>2697</td>
      <td>0.003337</td>
      <td>2.870400</td>
      <td>3.130230</td>
      <td>2.632130</td>
      <td>+</td>
      <td>0.916992</td>
      <td>0.565213</td>
      <td>0.104721</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>110761217</td>
      <td>SLC25A46</td>
      <td>5</td>
      <td>G</td>
      <td>A</td>
      <td>2697</td>
      <td>0.000185</td>
      <td>-0.048355</td>
      <td>0.045875</td>
      <td>0.050969</td>
      <td>-</td>
      <td>-1.054060</td>
      <td>4.668870</td>
      <td>0.821385</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>5</th>
      <td>110761271</td>
      <td>SLC25A46</td>
      <td>5</td>
      <td>G</td>
      <td>A</td>
      <td>2697</td>
      <td>0.000185</td>
      <td>0.259504</td>
      <td>0.191297</td>
      <td>0.352031</td>
      <td>+</td>
      <td>1.356550</td>
      <td>2.286360</td>
      <td>0.552966</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>6</th>
      <td>110761292</td>
      <td>SLC25A46</td>
      <td>5</td>
      <td>A</td>
      <td>G</td>
      <td>2697</td>
      <td>0.002225</td>
      <td>-2.064850</td>
      <td>2.313330</td>
      <td>1.843060</td>
      <td>-</td>
      <td>-0.892586</td>
      <td>0.657478</td>
      <td>0.174593</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>7</th>
      <td>110761301</td>
      <td>SLC25A46</td>
      <td>5</td>
      <td>T</td>
      <td>G</td>
      <td>2697</td>
      <td>0.000185</td>
      <td>-0.141477</td>
      <td>0.120936</td>
      <td>0.165506</td>
      <td>-</td>
      <td>-1.169850</td>
      <td>2.875560</td>
      <td>0.684137</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>8</th>
      <td>110761662</td>
      <td>SLC25A46</td>
      <td>5</td>
      <td>G</td>
      <td>T</td>
      <td>2697</td>
      <td>0.000371</td>
      <td>-0.008878</td>
      <td>0.130762</td>
      <td>0.000603</td>
      <td>-</td>
      <td>-0.067895</td>
      <td>2.765410</td>
      <td>0.980413</td>
      <td>1.0</td>
    </tr>
  </tbody>
</table>
</div>




```python
# reorder column
score_nonsyn_cols = list(score_nonsyn.columns)
score_nonsyn_new_cols = ["Gene", "CHROM", "POS"] + score_nonsyn_cols[3:]
score_nonsyn[score_nonsyn_new_cols]
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Gene</th>
      <th>CHROM</th>
      <th>POS</th>
      <th>REF</th>
      <th>ALT</th>
      <th>N_INFORMATIVE</th>
      <th>AF</th>
      <th>U</th>
      <th>V</th>
      <th>STAT</th>
      <th>DIRECTION</th>
      <th>EFFECT</th>
      <th>SE</th>
      <th>PVALUE</th>
      <th>FDR-Corrected</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110739267</td>
      <td>C</td>
      <td>A</td>
      <td>2697</td>
      <td>0.000742</td>
      <td>0.823836</td>
      <td>0.599113</td>
      <td>1.132850</td>
      <td>+</td>
      <td>1.375090</td>
      <td>1.291950</td>
      <td>0.287168</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110739354</td>
      <td>G</td>
      <td>A</td>
      <td>2697</td>
      <td>0.002225</td>
      <td>-0.172238</td>
      <td>1.856260</td>
      <td>0.015982</td>
      <td>-</td>
      <td>-0.092788</td>
      <td>0.733974</td>
      <td>0.899401</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110746300</td>
      <td>C</td>
      <td>A</td>
      <td>2697</td>
      <td>0.001112</td>
      <td>0.048595</td>
      <td>0.935050</td>
      <td>0.002526</td>
      <td>+</td>
      <td>0.051971</td>
      <td>1.034150</td>
      <td>0.959920</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110756712</td>
      <td>G</td>
      <td>A</td>
      <td>2697</td>
      <td>0.003337</td>
      <td>2.870400</td>
      <td>3.130230</td>
      <td>2.632130</td>
      <td>+</td>
      <td>0.916992</td>
      <td>0.565213</td>
      <td>0.104721</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110761217</td>
      <td>G</td>
      <td>A</td>
      <td>2697</td>
      <td>0.000185</td>
      <td>-0.048355</td>
      <td>0.045875</td>
      <td>0.050969</td>
      <td>-</td>
      <td>-1.054060</td>
      <td>4.668870</td>
      <td>0.821385</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>5</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110761271</td>
      <td>G</td>
      <td>A</td>
      <td>2697</td>
      <td>0.000185</td>
      <td>0.259504</td>
      <td>0.191297</td>
      <td>0.352031</td>
      <td>+</td>
      <td>1.356550</td>
      <td>2.286360</td>
      <td>0.552966</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>6</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110761292</td>
      <td>A</td>
      <td>G</td>
      <td>2697</td>
      <td>0.002225</td>
      <td>-2.064850</td>
      <td>2.313330</td>
      <td>1.843060</td>
      <td>-</td>
      <td>-0.892586</td>
      <td>0.657478</td>
      <td>0.174593</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>7</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110761301</td>
      <td>T</td>
      <td>G</td>
      <td>2697</td>
      <td>0.000185</td>
      <td>-0.141477</td>
      <td>0.120936</td>
      <td>0.165506</td>
      <td>-</td>
      <td>-1.169850</td>
      <td>2.875560</td>
      <td>0.684137</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>8</th>
      <td>SLC25A46</td>
      <td>5</td>
      <td>110761662</td>
      <td>G</td>
      <td>T</td>
      <td>2697</td>
      <td>0.000371</td>
      <td>-0.008878</td>
      <td>0.130762</td>
      <td>0.000603</td>
      <td>-</td>
      <td>-0.067895</td>
      <td>2.765410</td>
      <td>0.980413</td>
      <td>1.0</td>
    </tr>
  </tbody>
</table>
</div>


```python
score_nonsyn[score_nonsyn_new_cols].to_csv("result/score/SLC25A46.filtered.score.SingleScore.FDR.dup_cleaned.nonsyn.assoc", sep="\t", index=False)
```

<a id="7"></a>

# 7. Compound Heterozygous Fisher Test

1 control and 1 case was found with p.E79K/p.V211M compound heterozygosity (no homozygous found in either groups).

Generate 2x2 matrix

```R
# kernel: R 3.6
data <- matrix(c(1,1049,1,1046), ncol=2)
data
```


         [,1] [,2]
    [1,]    1    1
    [2,] 1049 1046

```R
# run fisher test
fisher.test(data)
```

    	Fisher's Exact Test for Count Data
    
    data:  data
    p-value = 1
    alternative hypothesis: true odds ratio is not equal to 1
    95 percent confidence interval:
      0.01269833 78.30123556
    sample estimates:
    odds ratio 
     0.9971414 
