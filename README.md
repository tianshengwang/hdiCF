# high-dimensional Iterative Causal Forest (hdiCF)

## hdiCF is free for nonprofit use

The hdiCF algorithm identifies important subgroups with heterogeneous treatment effects without prior knowledge of treatment-covariate interactions and predefined covariates

<img src = images/FIGURE1_hdiCF_full.jpg width=1000>

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

Please following the information for [iCF installation](https://github.com/tianshengwang/iCF).

**3. Run hdiCF on real-world claims data**

For simplicity, we focused on the ICD-10 era, included patients who initiated SGLT2i or GLP1RA treatment from October 2016 and followed them until December 2019. We compared the two-year risk difference (RD) of hospitalized heart failure (HHF) of initiating any sodium-glucose cotransporter-2 inhibitors (SGLT2i) versus glucagon-like peptide-1 receptor agonists (GLP1RA) using a 20% random sample of all fee-for-service U.S. Medicare beneficiaries who had parts A (inpatient), B (outpatient physician services), and D (dispensed prescription drugs) coverage for at least one month from October 2015 to December 2019. The details of the cohort are available in the mehtod paper (Wang et al.) 

***Step 1. High-dimensional feature identification***

***Step 2. Variable selection***
```{}
 PREPARE_HD <-function(train, dxgroup, atcgroup){
 # if (outcome=="adrd"){
    train0 <- train %>%  dplyr::mutate(Y = HHF_2yr_2yr,
                                       W = ifelse(SGLT==1,1,0),     
                                       sex=as.numeric(sex),
                                       race=as.numeric(race))
    train00 <- train0 %>% 
             dplyr::select(BENE_ID, Y, W, age, sex , race,
                               starts_with(c(paste0("dx",dxgroup) , 
                                             "cpt5", 
                                             paste0("atc", atcgroup)
                                             ),                  
                                             )
               ) %>% as.data.frame.matrix() 
  #remove columns with only one level
  train00 <- train00[, sapply(train00, function(col) length(unique(col))) > 1]   
  return(train00)
}

Train_BENEID_all <<- PREPARE_HD(Train2, 3, 4)
dat <- Train_BENEID_all %>% select(-c("BENE_ID"))
ID <-1:nrow(Train)
Train_ID <- cbind(Train, as.vector(ID)) %>% dplyr::rename (ID=`as.vector(ID)`)

vars_forest = colnames( dat %>% dplyr::select(-c("Y", "W" ))  )
 X <- dat[,vars_forest]
 Y <- as.vector( as.numeric( dat[,"Y"] ) )
 W <- as.vector( as.numeric( dat[,"W"] ) )

 cf_raw_key.tr <- CF_RAW_key(dat, 1, "hd", hdpct=0.90) 
 Y.hat  <<- cf_raw_key.tr$Y.hat                 
 W.hat  <<- cf_raw_key.tr$W.hat                 
 HTE_P_cf.raw <<- cf_raw_key.tr$HTE_P_cf.raw    
 varimp_cf  <- cf_raw_key.tr$varimp_cf          
 selected_cf.idx <<- cf_raw_key.tr$selected_cf.idx 
 GG_VI(varimp_cf, "Variable importance" )
 ```
 <img src = images/VI_allHD.png width=800>


 ```{}
colnames(X[,c(selected_cf.idx)])
X <<- X[,c(selected_cf.idx)]
dat <<- dat[,c("Y", "W", colnames(X)) ]
vars_forest = colnames( dat %>% dplyr::select(-c("Y", "W" ))  )
vars_catover2 <<- colnames(X[, -which(names(X) %in% c("sex")) ])

cf_raw_key.tr <- CF_RAW_key(Train, 1, "hd", hdpct=0.90) #use all selected variable in the 1st step
Y.hat  <<- cf_raw_key.tr$Y.hat
W.hat  <<- cf_raw_key.tr$W.hat
HTE_P_cf.raw <<- cf_raw_key.tr$HTE_P_cf.raw
HTE_P_cf.raw
varimp_cf  <- cf_raw_key.tr$varimp_cf
selected_cf.idx <<- cf_raw_key.tr$selected_cf.idx #MUST reselect important covaraites to run CF!!!!
length(selected_cf.idx)
colnames(X[,c(selected_cf.idx)])
GG_VI(varimp_cf, 'Variable importance for high-dimensional feature identification \n in SGLT2i vs GLP1RA cohort for HHF', colnames(X) )

 ```
 <img src = images/VI90percentile.png width=800>
 
 ***Step 3. Iplementation of iCF.***
 
 For details of [iCF algorithm: https://github.com/tianshengwang/iCF](https://github.com/tianshengwang/iCF) 
 
 ***To tune the leaf size, use different values for the minimum leaf size (MLS) to grow forests at various depths (D).***
 ```{}
#Specify the decimal position for continuous variables in the subgroup definition.
split_val_round_posi=0
#For real-world projects (not simulations where we know the truth), the truth is set as "Unknown".
truth.list <<- TRUTH("Unknown")
#Define categorical variables with more than two levels.
vars_catover2 <<- NA  
```
```{}
D2_MLS=MinLeafSizeTune(dat=dat, denominator=25, treeNo = 1000, iterationNo=100, split_val_round_posi=0, "D2", "steelblue4")
D2_MLS$depth_mean
D2_MLS$depth_gg
```
<img src = images/D2_MLS_tune_rwdHD.png width=350>

Notably, if you got this message "_Error: Can't subset columns that don't exist. x Column `parent_sign` doesn't exist._", it suggests the denominator used for developing is too small, leading to a too large MLS for D2 forest so that the tree does not even split (the node does not have a parent node). In this scenario, increasing the denominator will solve the problem.

```{}
D3_MLS=MinLeafSizeTune(dat=dat, denominator=45, treeNo = 1000, iterationNo=100, split_val_round_posi=0, "D3", "steelblue4")
D3_MLS$depth_mean
D3_MLS$depth_gg
```
<img src = images/D3_MLS_tune_rwdHD.png width=350>

```{}
D4_MLS=MinLeafSizeTune(dat=dat, denominator=65, treeNo = 1000, iterationNo=100, split_val_round_posi=0, "D4", "steelblue4")
D4_MLS$depth_mean
D4_MLS$depth_gg
```
<img src = images/D4_MLS_tune_rwdHD.png width=350>

```{}
D5_MLS=MinLeafSizeTune(dat=dat, denominator=85, treeNo = 1000, iterationNo=100, split_val_round_posi=0, "D5", "steelblue4")
D5_MLS$depth_mean
D5_MLS$depth_gg
```
<img src = images/D5_MLS_tune_rwdHD.png width=350>


***Run iCF on Medicare SGLT2i vs GLP1RA new user cohort***
```{}
leafsize <<- list(D5=D5_MLS$denominator, D4=D4_MLS$denominator, D3=D3_MLS$denominator, D2=D2_MLS$denominator)

iCFCV_B1000_i200_rwdHD <- iCFCV(dat=dat,K=5, treeNo=1000, iterationNo=100, min.split.var=4,
                              split_val_round_posi=0, P_threshold=0.1, variable_type = "hd", 
                              hdpct= 0.95, HTE_P_cf.raw = HTE_P_cf.raw) 

iCFCV_B1000_i200_rwdHD
```

If you have further questions or comments, please contact Dr. Tiansheng Wang: tianwang@unc.edu
