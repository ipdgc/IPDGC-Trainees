# The genetic architecture of resilience in Parkinson's disease

Sara Bandres-Ciga, Hui Liu, Mohammad Dehestani

#### Load modules and set up working directory
```
module load plink/2.0_alpha_1_final
module load plink
module load R

cd /data/LNG/saraB/PRS_resilience/FINAL/
```

### 1. IPDGC data - 7204 are cases and 9412 are controls
#### Remove problematic regions from Chang 2017 summary stats + filter by MAF 0.01
```
MHC region (hg19): chr6:28477797-33448354
plink --bfile HARDCALLS_PD_september_2018_no_cousins --keep noChangIDs.txt --make-bed --out noChangIDs_HARDCALLS_PD_september_2018_no_cousins
plink --bfile noChangIDs_HARDCALLS_PD_september_2018_no_cousins --remove-fam /data/LNG/saraB/NeuroX.fID.txt --make-bed --out noChangIDs_noNeuroX_HARDCALLS_PD_september_2018_no_cousins
plink --bfile noChangIDs_noNeuroX_HARDCALLS_PD_september_2018_no_cousins --exclude range /data/LNG/saraB/PRS_resilience/MHC.txt --maf 0.01 --make-bed --out ALL_no_MHC_noNeurox
```

#### Format fam file to match cov file
```
R
library("data.table")
test <- fread("ALL_no_MHC_noNeurox.fam", header = F)
colnames(test) <- c("FID", "IID", "MAT", "PAT", "SEX", "PHENO")
test$ID <- paste(test$FID, test$IID, sep = "_")
test2 <- test[,c("ID","ID","MAT", "PAT", "SEX", "PHENO")]
write.table(test2, file = "ALL_no_MHC_noNeurox.fam", quote = F, row.names = F, col.names = F, sep = "\t")
```

### GenoML to set up best threshold to define risk 
#### a) Munge datasets.
```
genoml discrete supervised munge --p 0.01 --prefix resilience_p1E2 --geno maf05_no_MHC --pheno pheno.csv --gwas gwas.csv --impute mean --skip_prune no --r2_cutoff 0.1 --feature_selection 500 --target_features snps.txt --confounders confounders.csv --adjust_data yes --adjust_normalize yes --umap_reduce no
genoml discrete supervised munge --p 0.001 --prefix resilience_p1E3 --geno maf05_no_MHC --pheno pheno.csv --gwas gwas.csv --impute mean --skip_prune no --r2_cutoff 0.1 --feature_selection 500 --target_features snps.txt --confounders confounders.csv --adjust_data yes --adjust_normalize yes --umap_reduce no
genoml discrete supervised munge --p 0.0001 --prefix resilience_p1E4 --geno maf05_no_MHC --pheno pheno.csv --gwas gwas.csv --impute mean --skip_prune no --r2_cutoff 0.1 --feature_selection 500 --target_features snps.txt --confounders confounders.csv --adjust_data yes --adjust_normalize yes --umap_reduce no
genoml discrete supervised munge --p 0.00001 --prefix resilience_p1E5 --geno maf05_no_MHC --pheno pheno.csv --gwas gwas.csv --impute mean --skip_prune no --r2_cutoff 0.1 --feature_selection 500 --target_features snps.txt --confounders confounders.csv --adjust_data yes --adjust_normalize yes --umap_reduce no
genoml discrete supervised munge --p 0.000001 --prefix resilience_p1E6 --geno maf05_no_MHC --pheno pheno.csv --gwas gwas.csv --impute mean --skip_prune no --r2_cutoff 0.1 --feature_selection 500 --target_features snps.txt --confounders confounders.csv --adjust_data yes --adjust_normalize yes --umap_reduce no
genoml discrete supervised munge --p 0.0000001 --prefix resilience_p1E7 --geno maf05_no_MHC --pheno pheno.csv --gwas gwas.csv --impute mean --skip_prune no --r2_cutoff 0.1 --feature_selection 500 --target_features snps.txt --confounders confounders.csv --adjust_data yes --adjust_normalize yes --umap_reduce no
genoml discrete supervised munge --p 0.00000001 --prefix resilience_p1E8 --geno maf05_no_MHC --pheno pheno.csv --gwas gwas.csv --impute mean --skip_prune no --r2_cutoff 0.1 --feature_selection 500 --target_features snps.txt --confounders confounders.csv --adjust_data yes --adjust_normalize yes --umap_reduce no
```

#### b) Train datasets.
```
genoml discrete supervised train --prefix risk_p1E2
genoml discrete supervised train --prefix risk_p1E3
genoml discrete supervised train --prefix risk_p1E4
genoml discrete supervised train --prefix risk_p1E5
genoml discrete supervised train --prefix risk_p1E6
genoml discrete supervised train --prefix risk_p1E7
genoml discrete supervised train --prefix risk_p1E8
```

#### c) Tune and cross-validate models.
```
genoml discrete supervised tune --prefix risk_p1E2 --max_tune 25
genoml discrete supervised tune --prefix risk_p1E3 --max_tune 10 --n_cv 3
genoml discrete supervised tune --prefix risk_p1E4 --max_tune 25
genoml discrete supervised tune --prefix risk_p1E5 --max_tune 25
genoml discrete supervised tune --prefix risk_p1E6 --max_tune 25
genoml discrete supervised tune --prefix risk_p1E7 --max_tune 25
genoml discrete supervised tune --prefix risk_p1E8 --max_tune 25
```

### Polygenic risk score
#### a) Extract risk clumped SNPs considering the bin 1e-3 nominated by GenoML - (selecting the bin that accounts for the most variance and performs better)
```
sed 's/SNP/MarkerName/g' risk_p1E3.variants_and_alleles.tab > risk_p1E3.variants_and_alleles_rename.tab

R
library(data.table)
Chang <- fread("nochr_METAANALYSIS_PdGeneAnd23AndMe1_Chang2017.tbl", header =T)
vars_1e <- fread("risk_p1E3.variants_and_alleles_rename.tab", header =F)
colnames(vars_1e) <- c("MarkerName", "Allele_final1")
total <- merge(Chang, vars_1e, by="MarkerName")
head(total)
outPut <- total[,c("MarkerName","Allele_final1","Effect")]
write.table(outPut, file = "Chang_1e3.toscore.txt", quote = F, row.names = F, col.names = F, sep = "\t")
```

#### b) Calculate scores in plink
```
plink --bfile maf05_no_MHC --score Chang_1e3.toscore.txt --out Chang_1e3.Training
```

#### c) Convert plink scores into z-scores and divide our data in percentiles - The goal is to identify affected and unaffected individuals at the upper 25th percentile
```
R
library("tidyr")
library("data.table")
library("dplyr")
library("plyr")
library("ggplot2")
temp_data <- read.table("Chang_1e3.Training.profile", header = T) 
temp_covs <- read.table("/data/LNG/saraB/WGS/toPRSice_phenosAndCovs_renamed.tab", header = T)
data <- merge(temp_data, temp_covs, by = "FID")
data$CASE <- data$PHENO.x - 1
meanControls <- mean(data$SCORE[data$CASE == 0])
sdControls <- sd(data$SCORE[data$CASE == 0])
data$zSCORE <- (data$SCORE - meanControls)/sdControls
Model <- glm(CASE ~ zSCORE + AGE + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = data, family = 'binomial')
summary(Model)

# > summary(Model)
# Call:
# glm(formula = CASE ~ zSCORE + AGE + sex + PC1 + PC2 + PC3 + PC4 + 
#     PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial", 
#     data = data)
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -2.1209  -1.2272   0.7537   0.9678   1.6546  
# 
# Coefficients:
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  2.170938   0.118050  18.390  < 2e-16 ***
# zSCORE       0.353827   0.020366  17.373  < 2e-16 ***
# AGE         -0.011371   0.001535  -7.410 1.26e-13 ***
# sex         -0.688730   0.040754 -16.900  < 2e-16 ***
# PC1         11.911682   4.512153   2.640  0.00829 ** 
# PC2         -3.604525   3.983284  -0.905  0.36551    
# PC3          2.898928   4.129909   0.702  0.48272    
# PC4         -8.850907   4.179341  -2.118  0.03419 *  
# PC5         -3.096120   3.991377  -0.776  0.43792    
# PC6         -8.713972   4.034589  -2.160  0.03079 *  
# PC7          3.966482   4.056270   0.978  0.32814    
# PC8          3.211387   4.020397   0.799  0.42442    
# PC9         -1.822958   4.083083  -0.446  0.65526    
# PC10         3.386459   3.959063   0.855  0.39235    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 14616  on 11032  degrees of freedom
# Residual deviance: 13887  on 11019  degrees of freedom
#   (5583 observations deleted due to missingness)
# AIC: 13915
# 
# Number of Fisher Scoring iterations: 4
```

#### d) Make a density plot
```
data$probDisease <- predict(Model, data, type = c("response"))
data$predicted <- ifelse(data$probDisease > 0.5, "DISEASE", "CONTROL")
data$reported <- ifelse(data$CASE == 1, "DISEASE","CONTROL")
densPlot <- ggplot(data, aes(probDisease, fill = reported, color = reported)) + geom_density(alpha = 0.5) + theme_bw()
ggsave(plot = densPlot, filename = "plotDensity.png", width = 10, height = 5, units = "in", dpi = 300)
```

#### e) Detect controls in the upper quantile of risk (25%) and based on those intervals we select cases with the same risk range
```
R
library("data.table")
library("dplyr")
temp_data <- fread("Chang_1e3.Training.profile", header = T) 
temp_covs <- fread("/data/LNG/saraB/WGS/toPRSice_phenosAndCovs_renamed.tab", header = T)
data <- merge(temp_data, temp_covs, by = "FID")
data$CASE <- data$PHENO - 1
meanControls <- mean(data$SCORE[data$CASE == 0])
sdControls <- sd(data$SCORE[data$CASE == 0])
data$zSCORE <- (data$SCORE - meanControls)/sdControls
max_Z_controls = max(data$zSCORE[data$CASE == 0])
decile_Z_controls = quantile(data$zSCORE[data$CASE == 0], probs = c(0.75))
training_data = subset(data, zSCORE >= decile_Z_controls & zSCORE <= max_Z_controls)
write.table(training_data, file = "Training_individuals_highriskquantile.txt", quote = F, row.names = F, sep = "\t")
```

#### f) Extract IPDGC individuals (cases/controls) in the upper quantile (25%) from plink binaries -  3011 are cases and 2353 are controls.
```
plink --bfile /data/LNG/saraB/PRS_resilience/ALL_no_MHC_noNeurox --keep Training_individuals_highriskquantile.txt --make-bed --out Training_cases_controls_noNeuroX_highriskquantile
```

#### g) Exclude 90 PD hits + 1Mb upstream/downstream to avoid following analyses
```
plink --bfile Training_cases_controls_noNeuroX_highriskquantile --exclude range /data/LNG/saraB/PRS_resilience/listofregions.txt --make-bed --out Training_cases_controls_highriskquantile_no_regions
```

#### h) Run a GWAS in the top 25% of risk, this will give you betas that are specific to the Training set. 
```
module load plink/2.0_alpha_1_final
plink2 --bfile Training_cases_controls_highriskquantile_no_regions --maf 0.05 --hwe 0.00001 --covar /data/LNG/saraB/WGS/noage_toPRSice_phenosAndCovs_renamed.tab --glm hide-covar --out training_resilience_GWAS 
```

#### i) Any hits? What's the lambda? Make a MH plot of pvalues per genomic region on the -log10 scale.
```
conda activate basicML
python
import pandas as pd
results = pd.read_csv("training_resilience_GWAS.PHENO1.glm.logistic", sep = "\t", engine='c') ## 4909240 vars remaining
results.describe()
hits = results[results['P'] < 0.00001]
print(hits)
```
```
R
library(data.table)
data <- fread("training_resilience_GWAS.PHENO1.glm.logistic", header = T)
n <- length(data$"P")
p <- data$"P"
x2obs <- qchisq(as.numeric(as.character(p)), 1, lower.tail = FALSE)
x2exp <- qchisq((1:n)/n, 1, lower.tail = FALSE)
lambda <- median(x2obs)/median(x2exp)
lambda -> ###[1] 1.03252 
```
```
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
df = results[results['P'] < 0.05]
df.rename(columns={'#CHROM':'CHR', 'POS':'BP'}, inplace=True)
df['-log10p'] = -np.log10(df['P'])
df['-log10p'] = np.where((df['-log10p'] > 20), 20, df['-log10p']) # Anything with -log10P > 20 is reduced to 20.
df = df.sort_values(['CHR','BP'])
df.reset_index(inplace=True, drop=True); df['i'] = df.index
plot = sns.relplot(data=df, x='i', y='-log10p', aspect=3, palette='bright', hue='CHR', legend=None) 
plt.axhline(8, ls='--', linewidth=3, color='red', alpha=0.25)
chrom_df=df.groupby('CHR')['i'].median()
plot.ax.set_xlabel('CHR'); plot.ax.set_xticks(chrom_df)
plot.ax.set_xticklabels(chrom_df.index)
plot.fig.suptitle('Resilience Manhattan (plot)')
plot.savefig("resilience_training_data_manhattan.png")
```

#### j) Train this model on GenoML - Split the data in 70-30% -> 70% to generate weights and 30% to fit the model
```
library(data.table)
fam <- fread("Training_cases_controls_highriskquantile_no_regions.fam", header = F)
colnames(fam) <- c("FID", "IID", "MAT", "PAT", "SEX", "PHENO")
fam$random <- rep(1:10, length = length(fam$PHENO))
fam$random <- rep(1:10)
fam_subsetted1 <- subset(fam, random > 7)
write.table(fam_subsetted1, "subsetted_Training_cases_controls_highriskquantile_no_regions.txt", quote = F, sep = "\t", row.names = F, col.names = F)
```
### Polygenic resilience
#### a) Extract 70% of the data - 2124 cases and 1632 controls remaining after main filters = 3756 samples
```
plink --bfile Training_cases_controls_highriskquantile_no_regions --remove subsetted_Training_cases_controls_highriskquantile_no_regions.txt --make-bed --out SEVENTY_training_resilience_GWAS
plink2 --bfile Training_cases_controls_highriskquantile_no_regions --remove subsetted_Training_cases_controls_highriskquantile_no_regions.txt --maf 0.05 --hwe 0.00001 --covar /data/LNG/saraB/WGS/noage_toPRSice_phenosAndCovs_renamed.tab --glm hide-covar --out SEVENTY_training_resilience_GWAS 
```

#### b) GenoML - Munge datasets.
```
genoml discrete supervised munge --p 0.01 --prefix resilience_p1E2_seventy --geno maf05_no_MHC_seventy --pheno pheno_seventy.csv --gwas gwas_seventy.csv --impute mean --skip_prune no --r2_cutoff 0.1 --feature_selection 500 --target_features ../snps.txt --addit confounders_seventy.csv
genoml discrete supervised munge --p 0.001 --prefix resilience_p1E3_seventy --geno maf05_no_MHC_seventy --pheno pheno_seventy.csv --gwas gwas_seventy.csv --impute mean --skip_prune no --r2_cutoff 0.1 --feature_selection 500 --target_features ../snps.txt --addit confounders_seventy.csv
genoml discrete supervised munge --p 0.0001 --prefix resilience_p1E4_seventy --geno maf05_no_MHC_seventy --pheno pheno_seventy.csv --gwas gwas_seventy.csv --impute mean --skip_prune no --r2_cutoff 0.1 --feature_selection 500 --target_features ../snps.txt --addit confounders_seventy.csv
genoml discrete supervised munge --p 0.00001 --prefix resilience_p1E5_seventy --geno maf05_no_MHC_seventy --pheno pheno_seventy.csv --gwas gwas_seventy.csv --impute mean --skip_prune no --r2_cutoff 0.1 --feature_selection 500 --target_features ../snps.txt --addit confounders_seventy.csv
genoml discrete supervised munge --p 0.000001 --prefix resilience_p1E6_seventy --geno maf05_no_MHC_seventy --pheno pheno_seventy.csv --gwas gwas_seventy.csv --impute mean --skip_prune no --r2_cutoff 0.1 --feature_selection 500 --target_features ../snps.txt --addit confounders_seventy.csv
genoml discrete supervised munge --p 0.0000001 --prefix resilience_p1E7_seventy --geno maf05_no_MHC_seventy --pheno pheno_seventy.csv --gwas gwas_seventy.csv --impute mean --skip_prune no --r2_cutoff 0.1 --feature_selection 500 --target_features ../snps.txt --addit confounders_seventy.csv
genoml discrete supervised munge --p 0.00000001 --prefix resilience_p1E8_seventy --geno maf05_no_MHC_seventy --pheno pheno_seventy.csv --gwas gwas_seventy.csv --impute mean --skip_prune no --r2_cutoff 0.1 --feature_selection 500 --target_features ../snps.txt --addit confounders_seventy.csv
```

#### c) Train datasets.
```
genoml discrete supervised train --prefix resilience_p1E2_seventy
genoml discrete supervised train --prefix resilience_p1E3_seventy
genoml discrete supervised train --prefix resilience_p1E4_seventy
genoml discrete supervised train --prefix resilience_p1E5_seventy
genoml discrete supervised train --prefix resilience_p1E6_seventy
genoml discrete supervised train --prefix resilience_p1E7_seventy
genoml discrete supervised train --prefix resilience_p1E8_seventy
```

#### d) Tune and cross-validate models.
```
genoml discrete supervised tune --prefix resilience_p1E2_seventy --max_tune 25
genoml discrete supervised tune --prefix resilience_p1E3_seventy --max_tune 25
genoml discrete supervised tune --prefix resilience_p1E4_seventy --max_tune 25
genoml discrete supervised tune --prefix resilience_p1E5_seventy --max_tune 25
genoml discrete supervised tune --prefix resilience_p1E6_seventy --max_tune 25
genoml discrete supervised tune --prefix resilience_p1E7_seventy --max_tune 25
genoml discrete supervised tune --prefix resilience_p1E8_seventy --max_tune 25
```

#### e) Run Polygenic Resilience Score in 30% of data using p-value threshold 1e-3 obtained in 70 % of the data (best performance in last model)
```
awk '{print $3, $6, $9}' SEVENTY_training_resilience_GWAS.PHENO1.glm.logistic > temp
```

```
R
library(data.table)
temp <- fread("temp", header = T)
temp$BETA <- log(temp$OR)
head(temp)
outPut <- temp[,c("ID","A1","BETA")]
temp2 <- fread("/data/LNG/saraB/PRS_resilience/ML/seventy/resilience_p1E3_seventy.variants_and_alleles.tab", header = F)
colnames(temp2) <- c("ID", "ALLELE")
total <- merge(outPut,temp2, by="ID")
total2 <- total[,c("ID","A1","BETA")]
write.table(total2, "SEVENTY_IPDGC_toscore_1e3.txt", quote = F, sep = "\t", row.names = F, col.names = F)
```

```
## 887 are cases and 721 are controls.
plink --bfile Training_cases_controls_highriskquantile_no_regions --keep subsetted_Training_cases_controls_highriskquantile_no_regions.txt --score SEVENTY_IPDGC_toscore_1e3.txt --out IPDGC_resilience_THIRTY 
```

```
python
import pandas as pd
import statsmodels.api as sm
prs_df = pd.read_csv("IPDGC_resilience_THIRTY.profile", delim_whitespace=True)
prs_df.head()
prs_reduced_df = prs_df[['IID','SCORE','PHENO']]
prs_reduced_mean = prs_reduced_df['SCORE'].mean()
prs_reduced_std = prs_reduced_df['SCORE'].std()
prs_reduced_df.rename(columns={'SCORE':'PRS'}, inplace=True)
prs_reduced_df['PRS_Z'] = (prs_reduced_df['PRS'] - prs_reduced_mean)/prs_reduced_std
cov_df = pd.read_csv("/data/LNG/saraB/WGS/toPRSice_phenosAndCovs_renamed.tab", delim_whitespace=True)
analysis_df = prs_reduced_df.merge(cov_df, on='IID', how='inner')
analysis_df['TEST'] = analysis_df['PHENO'] - 1
this_formula = "TEST ~ PRS_Z + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + sex"
res = sm.formula.glm(formula=this_formula, family=sm.families.Binomial(), data=analysis_df).fit() 
res.summary()
# 
#                  Generalized Linear Model Regression Results                  
# ==============================================================================
# Dep. Variable:                   TEST   No. Observations:                 1608
# Model:                            GLM   Df Residuals:                     1595
# Model Family:                Binomial   Df Model:                           12
# Link Function:                  logit   Scale:                          1.0000
# Method:                          IRLS   Log-Likelihood:                -1066.6
# Date:                Fri, 25 Jun 2021   Deviance:                       2133.1
# Time:                        14:52:56   Pearson chi2:                 1.61e+03
# No. Iterations:                     4                                         
# Covariance Type:            nonrobust                                         
# ==============================================================================
#                  coef    std err          z      P>|z|      [0.025      0.975]
# ------------------------------------------------------------------------------
# Intercept      1.0973      0.161      6.816      0.000       0.782       1.413
# PRS_Z          0.2101      0.053      3.993      0.000       0.107       0.313
# PC1           -8.0702     10.341     -0.780      0.435     -28.339      12.198
# PC2           33.7868     10.897      3.101      0.002      12.429      55.145
# PC3            9.9624     10.580      0.942      0.346     -10.775      30.700
# PC4            8.7896     10.701      0.821      0.411     -12.183      29.763
# PC5            0.0060     10.614      0.001      1.000     -20.796      20.808
# PC6          -30.7987     10.416     -2.957      0.003     -51.214     -10.384
# PC7            9.3624     10.294      0.909      0.363     -10.814      29.539
# PC8            2.2187     11.329      0.196      0.845     -19.987      24.424
# PC9           18.8174     10.727      1.754      0.079      -2.207      39.842
# PC10          -5.3844     10.162     -0.530      0.596     -25.301      14.532
# sex           -0.6016      0.103     -5.815      0.000      -0.804      -0.399
# ==============================================================================
```

#### f) Violin plot
```
R
library(ggplot2)
data <- read.table("IPDGC_resilience_THIRTY.profile", header = T) 
data$CASE <- data$PHENO - 1
meanControls <- mean(data$SCORE[data$CASE == 0])
sdControls <- sd(data$SCORE[data$CASE == 0])
data$zSCORE <- (data$SCORE - meanControls)/sdControls
data$CASE[data$CASE ==0] <- "Controls"
data$CASE[data$CASE ==1] <- "PD"
p <- ggplot(data, aes(x= reorder(as.factor(CASE), zSCORE), y=zSCORE, fill=as.factor(CASE))) +
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("lightblue", "orange")) + theme_bw() + ylab("Resilience IPDGC data (Z-transformed)") +xlab("") + theme(legend.position = "none")
ggsave("resilience_IPDGC.jpeg", dpi = 600, units = "in", height = 6, width = 6)
```

#############################################################################################################################

### 2. UKBB data - 2639 are cases and 14301 are controls

#### Format cov file 
```
Cov file -> /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_GENO_PD_CASE_CONTROL_with_PC.txt
awk '{print $1, $2}' /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_GENO_PD_CASE_CONTROL_with_PC.txt > indvs_to_keep.txt
```

#### Extract clumped SNPs from Phase I (RISK)
```
sed 's/:/ /g' resilience_p1E3.variants_and_alleles_rename.tab  > SNPS_risk_p1E3.variants_and_alleles_rename.tab
awk '{print $1, $2, $2, $3}' SNPS_risk_p1E3.variants_and_alleles_rename.tab > SNPS_risk_p1E3.variants_and_alleles_rename_extracted.tab
```
```
cat /data/LNG/saraB/PRS_resilience/chrnum.txt| while read LINE 
do
echo "THIS IS LINE" $LINE
plink2 --pfile /data/CARD/UKBIOBANK/FILTER_IMPUTED_DATA/chr$LINE.UKBB.EU.filtered --extract range /data/LNG/saraB/PRS_resilience/ML/SNPS_risk_p1E3.variants_and_alleles_rename_extracted.tab --keep indvs_to_keep.txt --make-bed --out  /data/LNG/saraB/PRS_resilience/FINAL/Testing.UKBB.clumped.chr$LINE
done
```
#### Merge data
```
plink --bfile Testing.UKBB.clumped.chr1 --merge-list list.txt --make-bed --out Testing.UKBB.clumped_ALL
```

#### Update pheno 
```
## 2639 are cases and 14301 are controls - European ancestry + unrelated
plink --bfile Testing.UKBB.clumped_ALL --pheno pheno_UKBB_V2.txt --out Testing.UKBB.clumped_ALL_updated --make-bed
sed 's/-9/1/g' Testing.UKBB.clumped_ALL_updated.fam > temp
mv temp Testing.UKBB.clumped_ALL_updated.fam 
```

### Polygenic risk score
#### a) Format bim file and calculate risk score using "Chang_1e3.toscore.txt" estimates.
```  
library(data.table)
bim <- fread("Testing.UKBB.clumped_ALL_updated.bim", header =F)
colnames(bim) <- c("chr", "chr_pos", "third", "pos", "a1", "a2")
bim$chr_conca_pos <- paste(bim$chr, bim$pos, sep = ":")
head(bim)
bim_final <- bim[,c("chr","chr_conca_pos","third", "pos", "a1", "a2")]
write.table(bim_final, file = "FINAL_Testing.UKBB.clumped_ALL_updated.bim", quote = F, row.names = F, col.names = F, sep = "\t")

```

```
cp Testing.UKBB.clumped_ALL_updated.fam FINAL_Testing.UKBB.clumped_ALL_updated.fam
cp Testing.UKBB.clumped_ALL_updated.bed FINAL_Testing.UKBB.clumped_ALL_updated.bed
```

```
## 1030 predictors loaded
plink --bfile FINAL_Testing.UKBB.clumped_ALL_updated --score Chang_1e3.toscore.txt --out PRS_risk_threshold_3e1.snps.Testing.UKBB
```

#### b) Convert plink scores into z-scores and divide data in percentiles - The goal is to identify affected and unaffected individuals at the upper 25th percentile. Compare ZSCORE in cases and controls by running glm
```
R
library("tidyr")
library("data.table")
library("dplyr")
library("plyr")
library("ggplot2")
temp_data <- read.table("PRS_risk_threshold_3e1.snps.Testing.UKBB.profile", header = T) 
temp_covs <- read.table("/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_GENO_PD_CASE_CONTROL_with_PC.txt", header = T)
data <- merge(temp_data, temp_covs, by = "FID")
data$CASE <- data$PHENO.x - 1
meanControls <- mean(data$SCORE[data$CASE == 0])
sdControls <- sd(data$SCORE[data$CASE == 0])
data$zSCORE <- (data$SCORE - meanControls)/sdControls
Model <- glm(CASE ~ zSCORE + AGE_OF_RECRUIT + GENETIC_SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = data, family = 'binomial')
summary(Model)

# Coefficients:
#                 Estimate Std. Error z value Pr(>|z|)    
# (Intercept)     3.953843   0.384230  10.290   <2e-16 ***
# zSCORE          0.215365   0.021678   9.935   <2e-16 ***
# AGE_OF_RECRUIT -0.094690   0.006069 -15.603   <2e-16 ***
# GENETIC_SEX     0.619085   0.044266  13.985   <2e-16 ***
# PC1            -1.088319   1.233096  -0.883   0.3775    
# PC2             1.645119   1.422492   1.157   0.2475    
# PC3             2.770907   1.426747   1.942   0.0521 .  
# PC4            -0.525241   1.423121  -0.369   0.7121    
# PC5             0.343819   1.422733   0.242   0.8090    
# PC6             1.355441   1.429791   0.948   0.3431    
# PC7            -0.114277   1.428759  -0.080   0.9363    
# PC8            -1.712985   1.429464  -1.198   0.2308    
# PC9             1.071797   1.430367   0.749   0.4537    
# PC10            0.217127   1.430843   0.152   0.8794    
```

#### c) Make density plot
```
data$probDisease <- predict(Model, data, type = c("response"))
data$predicted <- ifelse(data$probDisease > 0.5, "DISEASE", "CONTROL")
data$reported <- ifelse(data$CASE == 1, "DISEASE","CONTROL")
densPlot <- ggplot(data, aes(probDisease, fill = reported, color = reported)) + geom_density(alpha = 0.5) + theme_bw()
ggsave(plot = densPlot, filename = "plotDensity_UKBB.png", width = 10, height = 5, units = "in", dpi = 300)
```

#### d) Detect controls in the upper quantile of risk (25%) and based on those intervals we select cases with the same risk range
```
R
library("data.table")
library("dplyr")
temp_data <- fread("PRS_risk_threshold_3e1.snps.Testing.UKBB.profile", header = T) 
temp_covs <- fread("/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_GENO_PD_CASE_CONTROL_with_PC.txt", header = T)
data <- merge(temp_data, temp_covs, by = "FID")
data$CASE <- data$PHENO.x - 1
meanControls <- mean(data$SCORE[data$CASE == 0])
sdControls <- sd(data$SCORE[data$CASE == 0])
data$zSCORE <- (data$SCORE - meanControls)/sdControls
max_Z_controls = max(data$zSCORE[data$CASE == 0])
decile_Z_controls = quantile(data$zSCORE[data$CASE == 0], probs = c(0.75))
training_data = subset(data, zSCORE >= decile_Z_controls & zSCORE <= max_Z_controls)
write.table(training_data, file = "Testing_UKBB_individuals_highriskquantile.txt", quote = F, row.names = F, sep = "\t")
```

#### e) Extract individuals (cases/controls) in the in the upper quantile (25%) from plink binaries
```
## 847 are cases and 3576 are controls.
plink --bfile Testing.UKBB.clumped_ALL_updated --keep Testing_UKBB_individuals_highriskquantile.txt --make-bed --out Testing_UKBB_cases_controls_highriskquantile
```

#### f) Run a GWAS in this dataset of top 25% of risk, which will give us betas that are specific to the UKBB set 
```
awk '{print $1, $2, $2, $5}' training_resilience_GWAS.PHENO1.glm.logistic > SNPS_training_resilience_GWAS.PHENO1.glm.logistic
awk '{print $1, $2}' Testing_UKBB_individuals_highriskquantile.txt > extracted_subsetted_controls_and_cases_Testing_UKBB_cases_controls_highriskquantile.txt
```
```
cat /data/LNG/saraB/PRS_resilience/chrnum.txt| while read LINE 
do
echo "THIS IS LINE" $LINE
plink2 --pfile /data/CARD/UKBIOBANK/FILTER_IMPUTED_DATA/chr$LINE.UKBB.EU.filtered --keep /data/LNG/saraB/PRS_resilience/FINAL/extracted_subsetted_controls_and_cases_Testing_UKBB_cases_controls_highriskquantile.txt --extract range /data/LNG/saraB/PRS_resilience/ML/SNPS_training_resilience_GWAS.PHENO1.glm.logistic --exclude /data/LNG/saraB/PRS_resilience/FINAL/Testing.UKBB.resilience_GWAS_subset-merge.missnp --make-bed --out subsetted_Testing.UKBB.resilience_GWAS.chr$LINE
done
```

#### g) Merge data + update pheno
```
plink --bfile subsetted_Testing.UKBB.resilience_GWAS.chr1 --merge-list list2.txt --pheno pheno_UKBB_V2.txt --snps-only --make-bed --out Testing.UKBB.resilience_GWAS_subset
```
#### h) Run GWAS removing list of regions
```
plink2 --bfile Testing.UKBB.resilience_GWAS_subset --maf 0.05 --hwe 0.00001 --exclude range /data/LNG/saraB/PRS_resilience/listofregions.txt --covar /data/CARD/UKBIOBANK/PROJECTS/drug_mine_2020/covariates_phenome_final.txt --covar-name AGE_OF_RECRUIT,BATCH,GENETIC_SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 --glm hide-covar --out Testing_UKBB_resilience_GWAS_subset
```

#### i) Any hits on UKB dataset? What's the lambda? Make a MH plot of pvalues per genomic region on the -log10 scale.
```
conda activate basicML
python
import pandas as pd
results = pd.read_csv("Testing_UKBB_resilience_GWAS_subset.PHENO1.glm.logistic", sep = "\t", engine='c') ## 4754912 vars remaining
results.describe()
hits = results[results['P'] < 0.00001]
print(hits)
```
```
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
df = results[results['P'] < 0.05]
df.rename(columns={'#CHROM':'CHR', 'POS':'BP'}, inplace=True)
df['-log10p'] = -np.log10(df['P'])
df['-log10p'] = np.where((df['-log10p'] > 20), 20, df['-log10p']) # Anything with -log10P > 20 is reduced to 20.
df = df.sort_values(['CHR','BP'])
df.reset_index(inplace=True, drop=True); df['i'] = df.index
plot = sns.relplot(data=df, x='i', y='-log10p', aspect=3, palette='bright', hue='CHR', legend=None) 
plt.axhline(8, ls='--', linewidth=3, color='red', alpha=0.25)
chrom_df=df.groupby('CHR')['i'].median()
plot.ax.set_xlabel('CHR'); plot.ax.set_xticks(chrom_df)
plot.ax.set_xticklabels(chrom_df.index)
plot.fig.suptitle('Resilience Manhattan (plot) - UKB')
plot.savefig("resilience_UKBB_data_manhattan.png")
```
```
R
library(data.table)
data <- fread("Testing_UKBB_resilience_GWAS_subset.PHENO1.glm.logistic", header = T)
n <- length(data$"P")
p <- data$"P"
x2obs <- qchisq(as.numeric(as.character(p)), 1, lower.tail = FALSE)
x2exp <- qchisq((1:n)/n, 1, lower.tail = FALSE)
lambda <- median(x2obs)/median(x2exp)
### lambda -> [1] 1.022911
```

### Polygenic resilience
#### a) Run Polygenic Resilience Score in 30% of data using p-value threshold 1e-3 obtained in 70 % of the data (best performance in last model)
```
library(data.table)
bim <- fread("Testing.UKBB.resilience_GWAS_subset.bim", header = F)
colnames(bim) <- c("col1", "col2", "col3", "col4", "col5", "col6")
bim$col7 <- paste(bim$col1, bim$col4, sep = ":")
outPut <- bim[,c("col1","col7","col3", "col4", "col5", "col6")]
write.table(outPut, "Testing.UKBB.resilience_GWAS_subset_REFORMATTED.bim", quote = F, sep = "\t", row.names = F, col.names = F)
```

```
cp Testing.UKBB.resilience_GWAS_subset.fam Testing.UKBB.resilience_GWAS_subset_REFORMATTED.fam
cp Testing.UKBB.resilience_GWAS_subset.bed Testing.UKBB.resilience_GWAS_subset_REFORMATTED.bed
```

```
## 847 are cases and 3576 are controls. 802 valid predictors loaded.
plink --bfile Testing.UKBB.resilience_GWAS_subset_REFORMATTED --score /data/LNG/saraB/PRS_resilience/ML/SEVENTY_IPDGC_toscore_1e3.txt --out UKBB_resilience 
```

```
python
import pandas as pd
import statsmodels.api as sm
prs_df = pd.read_csv("UKBB_resilience.profile", delim_whitespace=True)
prs_df.head()
prs_reduced_df = prs_df[['IID','SCORE']]
prs_reduced_mean = prs_reduced_df['SCORE'].mean()
prs_reduced_std = prs_reduced_df['SCORE'].std()
prs_reduced_df.rename(columns={'SCORE':'PRS'}, inplace=True)
prs_reduced_df['PRS_Z'] = (prs_reduced_df['PRS'] - prs_reduced_mean)/prs_reduced_std
pheno_df = pd.read_csv("/data/LNG/saraB/PRS_resilience/pheno.ukbb.txt", delim_whitespace=True, names=['FID','IID','DUECES'])
pheno_reduced_df = pheno_df[['IID']]
pheno_reduced_df['PHENO'] = 1
cov_df = pd.read_csv("/data/CARD/UKBIOBANK/PROJECTS/drug_mine_2020/covariates_phenome_final.txt", delim_whitespace=True)
temp_df = prs_reduced_df.merge(pheno_reduced_df, on='IID', how='left')
temp_df['PHENO'].fillna(0, inplace=True)
analysis_df = temp_df.merge(cov_df, on='IID', how='inner')
this_formula = "PHENO ~ PRS_Z + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + AGE_OF_RECRUIT + TOWNSEND + GENETIC_SEX"
res = sm.formula.glm(formula=this_formula, family=sm.families.Binomial(), data=analysis_df).fit() 
res.summary()

# """
#                  Generalized Linear Model Regression Results                  
# ==============================================================================
# Dep. Variable:                  PHENO   No. Observations:                 4423
# Model:                            GLM   Df Residuals:                     4408
# Model Family:                Binomial   Df Model:                           14
# Link Function:                  logit   Scale:                          1.0000
# Method:                          IRLS   Log-Likelihood:                -1545.2
# Date:                Sun, 26 Sep 2021   Deviance:                       3090.4
# Time:                        16:40:33   Pearson chi2:                 4.51e+03
# No. Iterations:                     5                                         
# Covariance Type:            nonrobust                                         
# ==================================================================================
#                      coef    std err          z      P>|z|      [0.025      0.975]
# ----------------------------------------------------------------------------------
# Intercept          2.2312      0.871      2.562      0.010       0.525       3.938
# PRS_Z              0.0675      0.048      1.416      0.157      -0.026       0.161
# PC1               -0.0270      0.032     -0.851      0.395      -0.089       0.035
# PC2               -0.0568      0.032     -1.759      0.079      -0.120       0.006
# PC3                0.0100      0.031      0.319      0.750      -0.051       0.071
# PC4               -0.0154      0.023     -0.667      0.505      -0.061       0.030
# PC5                0.0142      0.010      1.403      0.161      -0.006       0.034
# PC6                0.0030      0.030      0.098      0.922      -0.057       0.063
# PC7                0.0134      0.027      0.499      0.618      -0.039       0.066
# PC8               -0.0253      0.027     -0.926      0.354      -0.079       0.028
# PC9                0.0074      0.011      0.668      0.504      -0.014       0.029
# PC10               0.0228      0.023      0.982      0.326      -0.023       0.068
# AGE_OF_RECRUIT    -0.0731      0.012     -5.852      0.000      -0.098      -0.049
# TOWNSEND           0.0052      0.017      0.305      0.761      -0.028       0.038
# GENETIC_SEX        0.5214      0.097      5.360      0.000       0.331       0.712
# ==================================================================================
```

#### b) Violin plot
```
R
library(ggplot2)
data <- read.table("UKBB_resilience.profile", header = T) 
data$CASE <- data$PHENO - 1
meanControls <- mean(data$SCORE[data$CASE == 0])
sdControls <- sd(data$SCORE[data$CASE == 0])
data$zSCORE <- (data$SCORE - meanControls)/sdControls
data$CASE[data$CASE ==0] <- "Controls"
data$CASE[data$CASE ==1] <- "PD"
p <- ggplot(data, aes(x= reorder(as.factor(CASE), zSCORE), y=zSCORE, fill=as.factor(CASE))) +
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("lightgreen", "red")) + theme_bw() + ylab("Resilience UKB data (Z-transformed)") +xlab("") + theme(legend.position = "none")
ggsave("resilience_UKB.jpeg", dpi = 600, units = "in", height = 6, width = 6)
```
################################################################################################################

### 3. AMP_PD data - version 2.5 - 2248 are cases and 2817 are controls

#### Format cov file
```
binary files-> /data/CARD/PD/AMP-PD/Plink/2021_v2_5release/rsID_pfiles_fromAMP/all_chrs_merged_maf_1e-3_rsIDs_callrate_sex_ancestry_EUR_related_het
cov file -> /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt
awk '{print $1, $2}'/data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt > tokeep.txt
awk '{print $1, $2, $10}'/data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt > pheno_ampv2.5.txt
```

#### Keep the right individuals and update pheno
```
plink --bfile /data/CARD/PD/AMP-PD/Plink/2021_v2_5release/rsID_pfiles_fromAMP/all_chrs_merged_maf_1e-3_rsIDs_callrate_sex_ancestry_EUR_related_het --keep tokeep.txt --make-bed --out AMP_PD_v2_clean_unrelated_eur
plink --bfile AMP_PD_v2_clean_unrelated_eur --pheno pheno_ampv2.5.txt --make-bed --out updated_AMP_PD_v2_clean_unrelated_eur
```

#### Reformat score file to run PRS (AMP-PD is on hg38)
```
sed 's/:/_/g' Chang_1e3.toscore.txt > temp_Chang_1e3.toscore.txt
```

```
R
library("data.table")
temp_data <- fread("/Users/bandrescigas/Desktop/PROJECTS/Alastair2019/mendelRandoAuto/SummaryStats/PD/chrPosRs.tab", header = F) 
colnames(temp_data) <- c("POS", "SNP")
temp_PRS <- fread("temp_Chang_1e3.toscore.txt", header = T)
colnames(temp_PRS) <- c("POS", "A1", "BETA")
data <- merge(temp_data, temp_PRS, by = "POS")
data2 <- subset(data, SNP != ".")
outPut <- data2[,c("SNP","A1","BETA")]
write.table(outPut, file = "AMP_formatted_PRS_risk_threshold_1e3.snps.toscore.txt", quote = F, row.names = F, col.names = F, sep = "\t")
```

```
awk '{print $1}' AMP_formatted_PRS_risk_threshold_1e3.snps.toscore.txt > SNPs_Chang2017_1e3.txt
```
#### Extract SNPs of interest to run PRS at 1E-3
```
988 variants remaining
plink --bfile updated_AMP_PD_v2_clean_unrelated_eur --extract SNPs_Chang2017_1e3.txt  --make-bed --out Testing.AMP_PD.clumped.ALL
```

### Polygenic risk score
#### a) Calculate score using estimates from Chang2017
```
988 valid predictors loaded.
plink --bfile Testing.AMP_PD.clumped.ALL --score AMP_formatted_PRS_risk_threshold_1e3.snps.toscore.txt --out PRS_risk_threshold_3e1.snps.Testing.AMP_PD
```

#### b) Convert plink scores into z-scores and divide our data in percentiles - The goal is to identify affected and unaffected individuals at the upper 25th percentile. Here we compare ZSCORE in cases and controls by running glm
```
R
library("tidyr")
library("data.table")
library("dplyr")
library("plyr")
library("ggplot2")
temp_data <- read.table("PRS_risk_threshold_3e1.snps.Testing.AMP_PD.profile", header = T) 
temp_covs <- read.table("/data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt", header = T)
data1 <- merge(temp_data, temp_covs, by = "FID")
temp_data2 <- read.table("Testing.AMP_PD.clumped.ALL.fam", header = F) 
colnames(temp_data2) <- c("FID", "IID", "MAT", "PAT", "SEX", "STATUS")
data2 <- merge(data1, temp_data2, by = "FID")
data <- subset(data2, STATUS != "-9")
data$CASE <- data$PHENO - 1
meanControls <- mean(data$SCORE[data$CASE == 0])
sdControls <- sd(data$SCORE[data$CASE == 0])
data$zSCORE <- (data$SCORE - meanControls)/sdControls
Model <- glm(CASE ~ zSCORE + SEX.x + PC1 + PC2 + PC3 + PC4 + PC5, data = data, family = 'binomial')
summary(Model)

# Call:
# glm(formula = CASE ~ zSCORE + SEX.x + PC1 + PC2 + PC3 + PC4 + 
#     PC5, family = "binomial", data = data)
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -1.7279  -1.0651  -0.6375   1.1285   2.1950  
# 
# Coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.68506    0.10312   6.643 3.07e-11 ***
# zSCORE       0.28342    0.03007   9.424  < 2e-16 ***
# SEX.x       -0.61364    0.05997 -10.232  < 2e-16 ***
# PC1         -1.91955    2.14263  -0.896    0.370    
# PC2          0.62414    0.97015   0.643    0.520    
# PC3         16.45497    2.94593   5.586 2.33e-08 ***
# PC4          0.45109    1.12928   0.399    0.690    
# PC5          0.45375    1.27301   0.356    0.722    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 6957.5  on 5064  degrees of freedom
# Residual deviance: 6473.8  on 5057  degrees of freedom
# AIC: 6489.8
# 
# Number of Fisher Scoring iterations: 4
```

#### c) Make a density plot
```
data$probDisease <- predict(Model, data, type = c("response"))
data$predicted <- ifelse(data$probDisease > 0.5, "DISEASE", "CONTROL")
data$reported <- ifelse(data$CASE == 1, "DISEASE","CONTROL")
densPlot <- ggplot(data, aes(probDisease, fill = reported, color = reported)) + geom_density(alpha = 0.5) + theme_bw()
ggsave(plot = densPlot, filename = "plotDensity_AMP_PD.png", width = 10, height = 5, units = "in", dpi = 300)
```

#### d) Detect controls in the upper quantile of risk (25%) and based on those intervals we select cases with the same risk range
```
R
library("data.table")
library("dplyr")
temp_data <- fread("PRS_risk_threshold_3e1.snps.Testing.AMP_PD.profile", header = T) 
data$CASE <- data$PHENO - 1
meanControls <- mean(data$SCORE[data$CASE == 0])
sdControls <- sd(data$SCORE[data$CASE == 0])
data$zSCORE <- (data$SCORE - meanControls)/sdControls
max_Z_controls = max(data$zSCORE[data$CASE == 0])
decile_Z_controls = quantile(data$zSCORE[data$CASE == 0], probs = c(0.75))
training_data = subset(data, zSCORE >= decile_Z_controls & zSCORE <= max_Z_controls)
write.table(training_data, file = "Testing_AMP_PD_individuals_highriskquantile.txt", quote = F, row.names = F, sep = "\t")
```

#### e) Prep to extract GWAS vars from training to speed up analysis
```
awk '{print $3}' training_resilience_GWAS.PHENO1.glm.logistic > temp
sed 's/:/_/g' temp > vars_GWAS_training.txt
```

```
R
library("data.table")
temp_data <- fread("/Users/bandrescigas/Desktop/PROJECTS/Alastair2019/mendelRandoAuto/SummaryStats/PD/chrPosRs.tab", header = F) 
colnames(temp_data) <- c("POS", "SNP")
temp_PRS <- fread("vars_GWAS_training.txt", header = T)
colnames(temp_PRS) <- c("POS")
data <- merge(temp_data, temp_PRS, by = "POS")
data2 <- subset(data, SNP != ".")
outPut <- data2[,c("SNP","POS")]
write.table(outPut, file = "GWAS_vars_toextract.txt", quote = F, row.names = F, col.names = F, sep = "\t")
```

```
awk '{print $1}' GWAS_vars_toextract.txt > SNPS_GWAS_vars_toextract.txt
```

#### f) Extract individuals (cases/controls) in the in the upper quantile (25%) from plink binaries 
```
## 798 are cases and 705 are controls
plink --bfile updated_AMP_PD_v2_clean_unrelated_eur --keep Testing_AMP_PD_individuals_highriskquantile.txt --extract SNPS_GWAS_vars_toextract.txt --make-bed --out Testing.AMP_PD.toGWAS
```

#### g) Run a GWAS
```
plink2 --bfile Testing.AMP_PD.toGWAS --maf 0.05 --hwe 0.00001 --covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt --covar-name SEX,PC1,PC2,PC3,PC4,PC5 --glm hide-covar --out Testing_AMP_PD_resilience_GWAS_subset
```

#### h) Any hits on AMP-PD dataset? What's the lambda? Make a MH plot of pvalues per genomic region on the -log10 scale.
```
conda activate basicML
python
import pandas as pd
results = pd.read_csv("Testing_AMP_PD_resilience_GWAS_subset.PHENO1.glm.logistic", sep = "\t", engine='c') ## 4462801 vars remaining
results.describe()
hits = results[results['P'] < 0.00001]
print(hits)
```

```
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
df = results[results['P'] < 0.05]
df.rename(columns={'#CHROM':'CHR', 'POS':'BP'}, inplace=True)
df['-log10p'] = -np.log10(df['P'])
df['-log10p'] = np.where((df['-log10p'] > 20), 20, df['-log10p']) # Anything with -log10P > 20 is reduced to 20.
df = df.sort_values(['CHR','BP'])
df.reset_index(inplace=True, drop=True); df['i'] = df.index
plot = sns.relplot(data=df, x='i', y='-log10p', aspect=3, palette='bright', hue='CHR', legend=None) 
plt.axhline(8, ls='--', linewidth=3, color='red', alpha=0.25)
chrom_df=df.groupby('CHR')['i'].median()
plot.ax.set_xlabel('CHR'); plot.ax.set_xticks(chrom_df)
plot.ax.set_xticklabels(chrom_df.index)
plot.fig.suptitle('Resilience Manhattan (plot) - AMP-PD')
plot.savefig("resilience_AMP_PD_data_manhattan.png")
```

```
R
library(data.table)
data <- fread("Testing_AMP_PD_resilience_GWAS_subset.PHENO1.glm.logistic", header = T)
data2 <- subset(data, P != "NA")
n <- length(data2$"P")
p <- data2$"P"
x2obs <- qchisq(as.numeric(as.character(p)), 1, lower.tail = FALSE)
x2exp <- qchisq((1:n)/n, 1, lower.tail = FALSE)
lambda <- median(x2obs)/median(x2exp)
lambda -> ###[[1] 1.011361
```

### Polygenic resilience score
#### a) Reformat resilience weights to rsIDs (AMP-PD is in hg38)
```
sed 's/:/_/g' SEVENTY_IPDGC_toscore_1e3.txt > temp_SEVENTY_IPDGC_toscore_1e3.txt
```
```
R
library("data.table")
temp_data <- fread("/Users/bandrescigas/Desktop/PROJECTS/Alastair2019/mendelRandoAuto/SummaryStats/PD/chrPosRs.tab", header = F) 
colnames(temp_data) <- c("POS", "SNP")
temp_PRS <- fread("temp_SEVENTY_IPDGC_toscore_1e3.txt", header = F)
colnames(temp_PRS) <- c("POS", "A1", "BETA")
data <- merge(temp_data, temp_PRS, by = "POS")
data2 <- subset(data, SNP != ".")
outPut <- data2[,c("SNP","A1","BETA")]
write.table(outPut, file = "AMP_formatted_PRS_resilience.snps.toscore.txt", quote = F, row.names = F, col.names = F, sep = "\t")
```

#### b) Run resilience profile 
```
## 798 are cases and 705 are controls, 794 valid predictors loaded.
plink --bfile Testing.AMP_PD.toGWAS --score AMP_formatted_PRS_resilience.snps.toscore.txt --out AMP_PD_resilience
```

```
python
import pandas as pd
import statsmodels.api as sm
prs_df = pd.read_csv("AMP_PD_resilience.profile", delim_whitespace=True)
prs_df.head()
prs_reduced_df = prs_df[['IID','SCORE']]
prs_reduced_mean = prs_reduced_df['SCORE'].mean()
prs_reduced_std = prs_reduced_df['SCORE'].std()
prs_reduced_df.rename(columns={'SCORE':'PRS'}, inplace=True)
prs_reduced_df['PRS_Z'] = (prs_reduced_df['PRS'] - prs_reduced_mean)/prs_reduced_std
cov_df = pd.read_csv("COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt", delim_whitespace=True)
pheno_df = pd.read_csv("pheno_ampv2.5.txt", delim_whitespace=True, header=None, names=['FID','IID','PHENO'])
temp_df = cov_df.merge(pheno_df, on='IID', how='inner')
analysis_df = prs_reduced_df.merge(temp_df, on='IID', how='left')
analysis_df['TEST'] = analysis_df['PHENO'] - 1
this_formula = "TEST ~ PRS_Z + PC1 + PC2 + PC3 + PC4 + PC5 + SEX"
res = sm.formula.glm(formula=this_formula, family=sm.families.Binomial(), data=analysis_df).fit() 
res.summary()

#                  Generalized Linear Model Regression Results                  
# ==============================================================================
# Dep. Variable:                   TEST   No. Observations:                 1503
# Model:                            GLM   Df Residuals:                     1495
# Model Family:                Binomial   Df Model:                            7
# Link Function:                  logit   Scale:                          1.0000
# Method:                          IRLS   Log-Likelihood:                -973.90
# Date:                Sun, 26 Sep 2021   Deviance:                       1947.8
# Time:                        18:37:45   Pearson chi2:                 1.50e+03
# No. Iterations:                     4                                         
# Covariance Type:            nonrobust                                         
# ==============================================================================
#                  coef    std err          z      P>|z|      [0.025      0.975]
# ------------------------------------------------------------------------------
# Intercept      1.3782      0.193      7.155      0.000       1.001       1.756
# PRS_Z          0.0232      0.054      0.427      0.669      -0.083       0.130
# PC1            4.5168      4.092      1.104      0.270      -3.504      12.537
# PC2            0.8375      1.777      0.471      0.637      -2.646       4.321
# PC3           25.1833      5.692      4.424      0.000      14.027      36.340
# PC4            1.6135      2.151      0.750      0.453      -2.602       5.829
# PC5           -0.7417      2.284     -0.325      0.745      -5.218       3.735
# SEX           -0.7258      0.109     -6.681      0.000      -0.939      -0.513
# ==============================================================================

## META-ANALYSIS (RESILIENCE SCORE)
# Make a table with estimates:
# COHORT	BETA	SE	P
# AMP_PD  	value	value	value
# UKB	value	0.052	value
# IPDGC	value	value	value
```

#### c) Violin plot
```
R
library(ggplot2)
data <- read.table("AMP_PD_resilience.profile", header = T) 
data$CASE <- data$PHENO - 1
meanControls <- mean(data$SCORE[data$CASE == 0])
sdControls <- sd(data$SCORE[data$CASE == 0])
data$zSCORE <- (data$SCORE - meanControls)/sdControls
data$CASE[data$CASE ==0] <- "Controls"
data$CASE[data$CASE ==1] <- "PD"
p <- ggplot(data, aes(x= reorder(as.factor(CASE), zSCORE), y=zSCORE, fill=as.factor(CASE))) +
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("yellow", "purple")) + theme_bw() + ylab("Resilience AMP-PD data (Z-transformed)") +xlab("") + theme(legend.position = "none")
ggsave("resilience_AMP_PD.jpeg", dpi = 600, units = "in", height = 6, width = 6)

```
### PEARSON correlations
#### IPDGC 
```
## 7204 are cases and 9412 are controls. 1060 valid predictors loaded.
plink --bfile /data/LNG/saraB/PRS_resilience/ALL_no_MHC_noNeurox --maf 0.05 --geno 0.01 --score /data/LNG/saraB/PRS_resilience/ML/Chang_1e3.toscore.txt --out ALL_IPDGC_RISK

## 7204 are cases and 9412 are controls. 863 valid predictors loaded.
plink --bfile /data/LNG/saraB/PRS_resilience/ALL_no_MHC_noNeurox --maf 0.05 --geno 0.01 --score /data/LNG/saraB/PRS_resilience/ML/SEVENTY_IPDGC_toscore_1e3.txt --out ALL_IPDGC_RESILIENCE
```

#### UKBB
```
## 2639 are cases and 14301 are controls.1030 valid predictors loaded.
plink --bfile FINAL_Testing.UKBB.clumped_ALL_updated --score Chang_1e3.toscore.txt --out ALL_UKBB_RISK

## 2639 are cases and 14301 are controls. 802 valid predictors loaded.
module load plink/2.0_alpha_1_final
cat /data/LNG/saraB/PRS_resilience/chrnum.txt| while read LINE 
do
echo "THIS IS LINE" $LINE
plink2 --pfile /data/CARD/UKBIOBANK/FILTER_IMPUTED_DATA/chr$LINE.UKBB.EU.filtered --keep /data/LNG/saraB/PRS_resilience/indvs_to_keep.txt --extract range /data/LNG/saraB/PRS_resilience/ML/SNPS_training_resilience_GWAS.PHENO1.glm.logistic --exclude /data/LNG/saraB/PRS_resilience/FINAL/Testing.UKBB.resilience_GWAS_subset-merge.missnp --make-bed --out Testing.UKBB.resilience_GWAS.$LINE.ALLindvs
done

plink --bfile Testing.UKBB.resilience_GWAS.1.ALLindvs --merge-list list_final.txt --pheno pheno_UKBB_V2.txt --snps-only --make-bed --out Testing.UKBB.resilience_GWAS_ALL
cp Testing.UKBB.resilience_GWAS_subset_REFORMATTED.bim Testing.UKBB.resilience_GWAS_ALL.bim

plink --bfile Testing.UKBB.resilience_GWAS_ALL --score /data/LNG/saraB/PRS_resilience/ML/SEVENTY_IPDGC_toscore_1e3.txt --out ALL_UKBB_RESILIENCE
```

#### AMP-PD
```
## 2248 are cases and 2817 are controls. 977 valid predictors loaded
plink --bfile updated_AMP_PD_v2_clean_unrelated_eur --maf 0.05 --geno 0.01 --score AMP_formatted_PRS_risk_threshold_1e3.snps.toscore.txt --out ALL_AMP_PD_RISK

## 2248 are cases and 2817 are controls. 793 valid predictors loaded.
plink --bfile updated_AMP_PD_v2_clean_unrelated_eur --maf 0.05 --geno 0.01 --score AMP_formatted_PRS_resilience.snps.toscore.txt --out ALL_AMP_PD_RESILIENCE
```

##############################################################################################

#### META for all cohorts
```
R
install.packages("rmeta")
library(rmeta)
library(data.table)
sumstats <- fread("summary.txt", header = T)
met <- meta.summaries(d = sumstats$BETA, se = sumstats$SE, method = c("fixed"), logscale = F, names = sumstats$COHORT) # logscale option set for OR.
met$test # gives you the Z and P for the meta-analysis.
```

###############################################################################################

#### Polygenic resilience score accounting for LRRK2 G2019S and GBA N370S status

```
## Identify LRRK2 G2019S and GBA N370S carriers
library(data.table)
fam <- fread("SOFTCALLS_PD_unfiltered_no_duplicates.fam", header =F)
colnames(fam) <- c("FID", "IID", "a", "b", "c", "d")
fam$conca <- paste(fam$FID, fam$IID, sep = "_")
head(fam)
fam_final <- fam[,c("conca","conca","a", "b", "c", "d")]
write.table(fam_final, file = "concat_SOFTCALLS_PD_unfiltered_no_duplicates.fam", quote = F, row.names = F, col.names = F, sep = "\t")

## IPDGC LRRK2 ## hg19
plink --bfile /data/LNG/saraB/important_files/concat_SOFTCALLS_PD_unfiltered_no_duplicates --chr 12 --from-bp 40734202 --to-bp 40734202 --recode A --out LRRK2_IPDGC
## IPDGC GBA - N370S carriers ## hg19
plink --bfile /data/LNG/saraB/important_files/concat_SOFTCALLS_PD_unfiltered_no_duplicates --chr 1 --from-bp 155205634 --to-bp 155205634 --recode A --out GBA_IPDGC 

## Generate cov file
library(data.table)
A <- fread("/data/LNG/saraB/important_files/phenosAndCovs_renamed.tab", header =T)
B <- fread("/data/LNG/saraB/PRS_resilience/LRRK2_IPDGC.raw", header =T)
total <- merge(A, B, by="FID")
C <- fread("/data/LNG/saraB/PRS_resilience/GBA_IPDGC.raw", header =T)
total_final <- merge(total, C, by="FID")
head(total_final)
write.table(total_final, file = "GBA_LRRK2_phenosAndCovs_renamed.tab", quote = F, row.names = F, sep = "\t")

## UKBB LRRK2 ## hg19
plink --bfile /data/CARD/UKBIOBANK/IMPUTED_DATA/LRRK2_area_v2 --chr 12 --from-bp 40734202 --to-bp 40734202 --recode A --out LRRK2_UKBB_indvs 
## UKBB GBA - N370S carriers ## hg19 - NO VARS FOUND!
plink --bfile /data/CARD/UKBIOBANK/IMPUTED_DATA/GBAv2 --chr 1 --from-bp 155205634 --to-bp 155205634 --recode A --out GBA_UKBB_indvs 

## Generate cov file
library(data.table)
A <- fread("/data/CARD/UKBIOBANK/PROJECTS/drug_mine_2020/covariates_phenome_final.txt", header =T)
B <- fread("/data/LNG/saraB/PRS_resilience/LRRK2_UKBB_indvs.raw", header =T)
total <- merge(A, B, by="FID")
head(total)
write.table(total, file = "LRRK2_covariates_phenome_final.txt", quote = F, row.names = F, sep = "\t")

## AMP-PD ## hg38
plink --bfile updated_AMP_PD_v2_clean_unrelated_eur --chr 12 --from-bp 40340400 --to-bp 40340400 --recode A --out LRRK2_AMP_highriskquantile 
## AMP-PD GBA - N370S carriers ## hg38
plink --bfile updated_AMP_PD_v2_clean_unrelated_eur --chr 1 --from-bp 155235843  --to-bp 155235843  --recode A --out GBA_AMP_highriskquantile 

library(data.table)
A <- fread("/data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt", header =T)
B <- fread("/data/LNG/saraB/PRS_resilience/FINAL/LRRK2_AMP_highriskquantile.raw", header =T)
total <- merge(A, B, by="FID")
C <- fread("/data/LNG/saraB/PRS_resilience/FINAL/GBA_AMP_highriskquantile.raw", header =T)
total_final <- merge(total, C, by="FID")
head(total_final)
write.table(total_final, file = "GBA_LRRK2_AMPv2.5.tab", quote = F, row.names = F, sep = "\t")

#### COURAGE-PD ## hg19
## COURAGE - PD LRRK2 - G2019S carriers ## hg19
plink --bfile COURAGE_EU --chr 12 --from-bp 40734202 --to-bp 40734202 --recode A --out LRRK2_GBA/LRRK2_COURAGE
## COURAGE - PD - N370S carriers ## hg19
plink --bfile COURAGE_EU --chr 1 --from-bp 155205634 --to-bp 155205634 --recode A --out LRRK2_GBA/GBA_COURAGE

## Generate cov file
A <- fread("euCOURAGE_COV", header =T)
B <- fread("LRRK2_COURAGE.raw", header =T)
total <- merge(A, B, by="FID")
C <- fread("GBA_COURAGE.raw", header =T)
total_final <- merge(total, C, by="FID")
colnames(total_final)[21] <- "rs34637584_A"
colnames(total_final)[27] <- "rs76763715_C"
colnames(total_final)[4] <- "IID.X1"
colnames(total_final)[25] <- "SEX.x"
write.table(total_final, file = "GBA_LRRK2_COURAGE_Covs.tab", quote = F, row.names = F, sep = "\t")

## Incorporate LRRK2 G2019S status in the polygenic resilience model as a covariate 

## IPDGC
python
import pandas as pd
import statsmodels.api as sm
prs_df = pd.read_csv("IPDGC_resilience_THIRTY.profile", delim_whitespace=True)
prs_df.head()
prs_reduced_df = prs_df[['IID','SCORE','PHENO']]
prs_reduced_mean = prs_reduced_df['SCORE'].mean()
prs_reduced_std = prs_reduced_df['SCORE'].std()
prs_reduced_df.rename(columns={'SCORE':'PRS'}, inplace=True)
prs_reduced_df['PRS_Z'] = (prs_reduced_df['PRS'] - prs_reduced_mean)/prs_reduced_std
cov_df = pd.read_csv("GBA_LRRK2_phenosAndCovs_renamed.tab", delim_whitespace=True)
analysis_df = prs_reduced_df.merge(cov_df, on='IID', how='inner')
analysis_df['TEST'] = analysis_df['PHENO'] - 1
this_formula = "TEST ~ PRS_Z + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + sex + rs34637584_A"
res = sm.formula.glm(formula=this_formula, family=sm.families.Binomial(), data=analysis_df).fit() 
res.summary()

## UKBB
python
import pandas as pd
import statsmodels.api as sm
prs_df = pd.read_csv("UKBB_resilience.profile", delim_whitespace=True)
prs_df.head()
prs_reduced_df = prs_df[['IID','SCORE']]
prs_reduced_mean = prs_reduced_df['SCORE'].mean()
prs_reduced_std = prs_reduced_df['SCORE'].std()
prs_reduced_df.rename(columns={'SCORE':'PRS'}, inplace=True)
prs_reduced_df['PRS_Z'] = (prs_reduced_df['PRS'] - prs_reduced_mean)/prs_reduced_std
pheno_df = pd.read_csv("/data/LNG/saraB/PRS_resilience/pheno.ukbb.txt", delim_whitespace=True, names=['FID','IID','DUECES'])
pheno_reduced_df = pheno_df[['IID']]
pheno_reduced_df['PHENO'] = 1
cov_df = pd.read_csv("LRRK2_covariates_phenome_final.txt", delim_whitespace=True)
temp_df = prs_reduced_df.merge(pheno_reduced_df, on='IID', how='left')
temp_df['PHENO'].fillna(0, inplace=True)
analysis_df = temp_df.merge(cov_df, on='IID', how='inner')
this_formula = "PHENO ~ PRS_Z + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + AGE_OF_RECRUIT + TOWNSEND + GENETIC_SEX + rs34637584_A"
res = sm.formula.glm(formula=this_formula, family=sm.families.Binomial(), data=analysis_df).fit() 
res.summary()

## AMP-PD
python
import pandas as pd
import statsmodels.api as sm
prs_df = pd.read_csv("AMP_PD_resilience.profile", delim_whitespace=True)
prs_df.head()
prs_reduced_df = prs_df[['IID','SCORE']]
prs_reduced_mean = prs_reduced_df['SCORE'].mean()
prs_reduced_std = prs_reduced_df['SCORE'].std()
prs_reduced_df.rename(columns={'SCORE':'PRS'}, inplace=True)
prs_reduced_df['PRS_Z'] = (prs_reduced_df['PRS'] - prs_reduced_mean)/prs_reduced_std
cov_df = pd.read_csv("GBA_LRRK2_AMPv2.5.tab", delim_whitespace=True)
pheno_df = pd.read_csv("pheno_ampv2.5.txt", delim_whitespace=True, header=None, names=['FID','IID','PHENO'])
temp_df = cov_df.merge(pheno_df, on='IID', how='inner')
analysis_df = prs_reduced_df.merge(temp_df, on='IID', how='left')
analysis_df['TEST'] = analysis_df['PHENO'] - 1
this_formula = "TEST ~ PRS_Z + PC1 + PC2 + PC3 + PC4 + PC5 + SEX + rs34637584_A"
res = sm.formula.glm(formula=this_formula, family=sm.families.Binomial(), data=analysis_df).fit() 
res.summary()

## COURAGE - PD
temp_data <- fread("COURAGE_resilience_flipped.profile")
temp_covs <- fread("GBA_LRRK2_COURAGE_Covs.tab")
data <- merge(temp_data, temp_covs, by = c("FID"="FID"))
data$CASE <- data$PHENO.y - 1
meanSCORE <- mean(data$SCORE)
sdSCORE <- sd(data$SCORE)
data$zSCORE <- (data$SCORE - meanSCORE)/sdSCORE
Model_LRRK2 <- glm(CASE ~ zSCORE + AGE+ SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + rs34637584_A, data = data, family = 'binomial')
summary(Model_LRRK2)

## Incorporate GBA - N370S status in the polygenic resilience model as a covariate 

## IPDGC
python
import pandas as pd
import statsmodels.api as sm
prs_df = pd.read_csv("IPDGC_resilience_THIRTY.profile", delim_whitespace=True)
prs_df.head()
prs_reduced_df = prs_df[['IID','SCORE','PHENO']]
prs_reduced_mean = prs_reduced_df['SCORE'].mean()
prs_reduced_std = prs_reduced_df['SCORE'].std()
prs_reduced_df.rename(columns={'SCORE':'PRS'}, inplace=True)
prs_reduced_df['PRS_Z'] = (prs_reduced_df['PRS'] - prs_reduced_mean)/prs_reduced_std
cov_df = pd.read_csv("/data/LNG/saraB/WGS/toPRSice_phenosAndCovs_renamed.tab", delim_whitespace=True)
analysis_df = prs_reduced_df.merge(cov_df, on='IID', how='inner')
analysis_df['TEST'] = analysis_df['PHENO'] - 1
this_formula = "TEST ~ PRS_Z + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + sex + rs76763715_C"
res = sm.formula.glm(formula=this_formula, family=sm.families.Binomial(), data=analysis_df).fit() 
res.summary()

## AMP-PD
python
import pandas as pd
import statsmodels.api as sm
prs_df = pd.read_csv("AMP_PD_resilience.profile", delim_whitespace=True)
prs_df.head()
prs_reduced_df = prs_df[['IID','SCORE']]
prs_reduced_mean = prs_reduced_df['SCORE'].mean()
prs_reduced_std = prs_reduced_df['SCORE'].std()
prs_reduced_df.rename(columns={'SCORE':'PRS'}, inplace=True)
prs_reduced_df['PRS_Z'] = (prs_reduced_df['PRS'] - prs_reduced_mean)/prs_reduced_std
cov_df = pd.read_csv("COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt", delim_whitespace=True)
pheno_df = pd.read_csv("pheno_ampv2.5.txt", delim_whitespace=True, header=None, names=['FID','IID','PHENO'])
temp_df = cov_df.merge(pheno_df, on='IID', how='inner')
analysis_df = prs_reduced_df.merge(temp_df, on='IID', how='left')
analysis_df['TEST'] = analysis_df['PHENO'] - 1
this_formula = "TEST ~ PRS_Z + PC1 + PC2 + PC3 + PC4 + PC5 + SEX + rs76763715_C"
res = sm.formula.glm(formula=this_formula, family=sm.families.Binomial(), data=analysis_df).fit() 
res.summary()

## COURAGE - PD
temp_data <- fread("COURAGE_resilience_flipped.profile")
temp_covs <- fread("GBA_LRRK2_COURAGE_Covs.tab")
data <- merge(temp_data, temp_covs, by = c("FID"="FID"))
data$CASE <- data$PHENO.y - 1
meanSCORE <- mean(data$SCORE)
sdSCORE <- sd(data$SCORE)
data$zSCORE <- (data$SCORE - meanSCORE)/sdSCORE
Model_GBA <- glm(CASE ~ zSCORE + AGE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + rs76763715_C, data = data, family = 'binomial')
summary(Model_GBA)

```
#############################################################################################################################

#### Heritability analyses through LDSC

```
## Reformat summary statistics for LDSC
library("data.table")
library("dplyr")
library("tidyr")
library("plyr")
PD <- fread("Meta_GWAS_4cohorts_LD", header = T)
TEMP <- fread("test_chrPosRs.tab", header = F) 
colnames(TEMP) <- c("SNP", "snpid")
merged <- merge(TEMP, PD, by="SNP")
merged$OR <- as.numeric(merged$OR)
merged$BETA <- log(merged$OR)
merged$SE_fixed_qnorm<-abs(log(merged$OR)/qnorm(merged$P/2))
merged$Zscore <- merged$BETA/merged$SE_fixed_qnorm
merged$P.value <- merged$P
merged$N <- 16292
outPut <- merged[,c("snpid","A1","A2","Zscore", "N", "P.value")]
outPut2 <- subset(outPut, snpid != "<NA>")
write.table(outPut2, file = "PD_Resilience_toLDhub.txt", quote = F, row.names = F, sep = "\t")

## Copy to biowulf
scp /Users/bandrescigas/Desktop/PD_Resilience_toLDhub.txt bandrescigas@biowulf2.nih.gov:/data/LNG/saraB/

## Run LDSC package to calculate heritability
munge_sumstats.py --out Resilience_all_LDSC --sumstats /data/LNG/saraB/PD_Resilience_toLDhub.txt.gz \
--merge-alleles /data/LNG/saraB/1kg_eur/eur_w_ld_chr/w_hm3.snplist \
--snp snpid --a1 A1 --a2 A2 --p P.value 

ldsc.py \
--h2 /data/LNG/saraB/PD_Resilience_toLDhub.txt.gz \
--ref-ld-chr /data/LNG/saraB/1kg_eur/eur_w_ld_chr/ \
--w-ld-chr /data/LNG/saraB/1kg_eur/eur_w_ld_chr/ \
--out Resilience_all_LDSC
less Resilience_all_LDSC.log 

```

#############################################################################################################################

#### Data visualization - Functional Enrichment

```
module load R/4.0
R --vanilla --no-save
library(data.table)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)

## enrichment is calculated by this formula -> (intersection_size/term_size)*100
data <- fread("gProfiler_hsapiens_intersections.csv",header=T)
data$newOrder <- seq(1:30)
data$source <- factor(data$source, levels=c("GO:BP","GO:CC","GO:MF"))

bs=seq(1,4,0.5)

p1 = ggplot(data[data$source == "GO:BP"],) +
    geom_point(aes(x = enrichment,y = reorder(term_name, newOrder), size = `negative_log10_of_adjusted_p_value`), color= "#d55e00",alpha=0.8) +
    facet_grid(source ~ .) +
    scale_color_brewer(palette="Set1") +
scale_size(range = c(0.3, 3), breaks=as.vector(bs), labels=c("1.0","1.5","2.0","2.5","3.0","3.5","4.0")) +
    scale_x_continuous(limits = c(0,16), breaks = seq(0, 15, by = 5)) +
    xlab("Enrichment (%)") +
    ylab("Pathways") +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black", size = 6, face = "plain"),
        axis.text.y = element_text(color = "black", size = 6, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 8, face = "plain"),
        axis.title.y = element_text(color = "black", size = 8, face = "plain")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(strip.background =element_rect(fill="#d55e00")) +
    theme(strip.text = element_text(colour = 'white', size=8,face="bold")) +
    theme(legend.title = element_text(color = "black", size = 6,face="bold"),
          legend.text = element_text(color = "black", size = 6)) +
    theme(plot.margin = unit(c(1, 3, 1, 1), "cm"))
ggsave(paste("PathwayPlot1",".jpeg", sep = ""), plot = p1, device = "jpeg", scale = 1, width = 7, height = 3.3, units = "in", dpi = 300, limitsize = TRUE)
ggsave(paste("PathwayPlot1",".pdf", sep = ""), plot = p1, device = "pdf", scale = 1, width = 7, height = 3.3, units = "in", dpi = 300, limitsize = TRUE)

p2 = ggplot(data[data$source == "GO:CC"],) +
  geom_point(aes(x = enrichment,y = reorder(term_name, newOrder), size = `negative_log10_of_adjusted_p_value`), color= "#0072b2",alpha=0.8) +
  facet_grid(source ~ .) +
  scale_color_brewer(palette="Set1") +
  scale_size(range = c(0.3, 3), breaks=as.vector(bs), labels=c("1.0","1.5","2.0","2.5","3.0","3.5","4.0")) +
  scale_x_continuous(limits = c(0,16), breaks = seq(0, 15, by = 5)) +
  xlab("Enrichment (%)") +
  ylab("Pathways") +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 6, face = "plain"),
        axis.text.y = element_text(color = "black", size = 6, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 8, face = "plain"),
        axis.title.y = element_text(color = "black", size = 8, face = "plain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.background =element_rect(fill="#0072b2")) +
  theme(strip.text = element_text(colour = 'white', size=8,face="bold")) +
  theme(legend.title = element_text(color = "black", size = 6,face="bold"),
        legend.text = element_text(color = "black", size = 6)) +
  theme(plot.margin = unit(c(1, 3, 1, 1), "cm"))
ggsave(paste("PathwayPlot2",".jpeg", sep = ""), plot = p2, device = "jpeg", scale = 1, width = 7, height = 3.3, units = "in", dpi = 300, limitsize = TRUE)
ggsave(paste("PathwayPlot2",".pdf", sep = ""), plot = p2, device = "pdf", scale = 1, width = 7, height = 3.3, units = "in", dpi = 300, limitsize = TRUE)

p3 = ggplot(data[data$source == "GO:MF"],) +
  geom_point(aes(x = enrichment,y = reorder(term_name, newOrder), size = `negative_log10_of_adjusted_p_value`), color= "#009e73",alpha=0.8) +
  facet_grid(source ~ .) +
  scale_color_brewer(palette="Set1") +
  scale_size(range = c(0.3, 3), breaks=as.vector(bs), labels=c("1.0","1.5","2.0","2.5","3.0","3.5","4.0")) +
  scale_x_continuous(limits = c(0,16), breaks = seq(0, 15, by = 5)) +
  xlab("Enrichment (%)") +
  ylab("Pathways") +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 6, face = "plain"),
        axis.text.y = element_text(color = "black", size = 6, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 8, face = "plain"),
        axis.title.y = element_text(color = "black", size = 8, face = "plain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.background =element_rect(fill="#009e73")) +
  theme(strip.text = element_text(colour = 'white', size=8,face="bold")) +
  theme(legend.title = element_text(color = "black", size = 6,face="bold"),
        legend.text = element_text(color = "black", size = 6)) +
  theme(plot.margin = unit(c(1, 3, 1, 1), "cm"))
ggsave(paste("PathwayPlot3",".jpeg", sep = ""), plot = p3, device = "jpeg", scale = 1, width = 6.5, height = 3.3, units = "in", dpi = 300, limitsize = TRUE)
ggsave(paste("PathwayPlot3",".pdf", sep = ""), plot = p3, device = "pdf", scale = 1, width = 6.5, height = 3.3, units = "in", dpi = 300, limitsize = TRUE)

```

#### Data visualization - Forest plot

```
dat <- fread("Meta_forest.txt")
level_order <- c('META','COURAGE-PD','UKBB',  'AMP-PD','IPDGC')
resilience_plot <- ggplot(data=dat, aes(x = factor(COHORT,levels = level_order),
                         y = BETA,
                         ymin = L95,
                         ymax = U95)
) +
  geom_pointrange(
    aes(ymin = L95,
        ymax = U95),
    cex = 0.7
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 2) +
  theme(plot.title = element_text(size = 20,
                                  face = "bold"),
        axis.text.y = element_text(size = 8,
                                   face = 'bold'),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(face="bold"),
        axis.title = element_text(size = 16,
                                  face="bold"),
        legend.position = "none"
  ) +
  xlab(' ') +
  ylab("Beta coefficient") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Polygenic resilience score") +
  theme(plot.title = element_text(hjust=0.5))


ggsave("Forestplot-resilience.jpg", plot = resilience_plot,
       width = 8, height = 5, 
       units = "in", # other options c("in", "cm", "mm"), 
       dpi = 500)

```

