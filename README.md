# high-dimensional Iterative Causal Forest (hdiCF)
----------------------------------------------------------------------------------
The hdiCF algorithm identifies important subgroups with heterogeneous treatment effects without prior knowledge of treatment-covariate interactions and predefined covariates

<img src = images/FIG1_19Mar2023.jpg width=1000>

**Citation**

Wang T, Pate V, Wyss R, Buse JB, Kosorok MR, Stürmer T. High-dimensional Iterative Causal Forest (hdiCF) for Subgroup Identification Using Health Care Claims Data. Am J Epidemiol. 2023 (In Preparation).

**1. R packages recommended**
```{r packages, include=FALSE}
library(MASS)
library(grf)
library(tidyverse)
library(rlang)
library(rlist)
library(plyr)
library(caret)
library(caTools)
library(listdtr)
library(randomForest)
library(ggplot2)
library(ggridges)
library(data.table)
library(grid)
library(broom)
library(rstatix)
library(DMwR)
library(knitr)
library(Rfast)
library(spaMM)
```
**2. Installation**

Please following the information for [iCF instaliation](https://github.com/tianshengwang/iCF).

**3. Run hdiCF on real-world claims data**

We compared the two-year risk difference (RD) of hospitalized heart failure (HHF) of initiating any sodium-glucose cotransporter-2 inhibitors (SGLT2i) versus glucagon-like peptide-1 receptor agonists (GLP1RA) using a 20% random sample of all fee-for-service U.S. Medicare beneficiaries who had parts A (inpatient), B (outpatient physician services), and D (dispensed prescription drugs) coverage for at least one month from January 2012 to December 2017. The details of the cohort are available in the mehtod paper (Wang et al.) 

***Step 1. Run raw causal forest to predict outcome (Y.hat), propensity score (W.hat), and select variables***
```{}
 vars_forest = colnames( dat %>% dplyr::select(-c("Y", "W" ))  )
 X <- dat[,vars_forest]
 Y <- as.vector( as.numeric( dat[,"Y"] ) )
 W <- as.vector( as.numeric( dat[,"W"] ) )
 cf_raw_key.tr <- CF_RAW_key(dat, 1, "non-hd", hdpct=0.90) 
 Y.hat  <<- cf_raw_key.tr$Y.hat                 
 W.hat  <<- cf_raw_key.tr$W.hat                 
 HTE_P_cf.raw <<- cf_raw_key.tr$HTE_P_cf.raw    
 varimp_cf  <- cf_raw_key.tr$varimp_cf          
 selected_cf.idx <<- cf_raw_key.tr$selected_cf.idx 
 GG_VI(varimp_cf, "Variable importance" )
 ```
 <img src = images/VI_SGLTvGLP_HHF2y.jpeg width=800>
 
 ***Step 2. To tune the leaf size, use different values for the minimum leaf size (MLS) to grow forests at various depths (D).***
 ```{}
#Specify the decimal position for continuous variables in the subgroup definition.
split_val_round_posi=0
#For real-world projects (not simulations where we know the truth), the truth is set as "Unknown".
truth.list <<- TRUTH("Unknown")
#Define categorical variables with more than two levels.
vars_catover2 <<- NA  
```
```{}
D2_MLS=MinLeafSizeTune(dat=dat, denominator=25, treeNo = 1000, iterationNo=100, split_val_round_posi=0, "D2")
D2_MLS$depth_mean
D2_MLS$depth_gg
```
<img src = images/D2_MLS_tune.png width=350>

```{}
D3_MLS=MinLeafSizeTune(dat=dat, denominator=45, treeNo = 1000, iterationNo=100, split_val_round_posi=0, "D3")
D3_MLS$depth_mean
D3_MLS$depth_gg
```
<img src = images/D3_MLS_tune.png width=350>

```{}
D4_MLS=MinLeafSizeTune(dat=dat, denominator=65, treeNo = 1000, iterationNo=100, split_val_round_posi=0, "D4")
D4_MLS$depth_mean
D4_MLS$depth_gg
```
<img src = images/D4_MLS_tune.png width=350>

```{}
D5_MLS=MinLeafSizeTune(dat=dat, denominator=85, treeNo = 1000, iterationNo=100, split_val_round_posi=0, "D5")
D5_MLS$depth_mean
D5_MLS$depth_gg
```
<img src = images/D5_MLS_tune.png width=350>

*Note that despite using a smaller minimum leaf size (MLS), the causal forests do not grow deeper due to the presence of a strong 3-way interaction (W:X1:X3) in the simulated data set.* 

***Step 3. Implement iCF on Medicare SGLT2i vs GLP1RA new user cohort***
```{}
leafsize <<- list(D5=85, D4=65, D3=45, D2=25)

iCFCV_B1000_i200_rwd <- iCFCV(dat=dat,K=5, treeNo=1000, iterationNo=100, min.split.var=4,
                              split_val_round_posi=0, P_threshold=0.1, variable_type = "non-HD", 
                              hdpct= 0.95, HTE_P_cf.raw = HTE_P_cf.raw) 

iCFCV_B1000_i200_rwd
```

If you have further questions or comments, please contact Dr. Tiansheng Wang: tianwang@unc.edu
