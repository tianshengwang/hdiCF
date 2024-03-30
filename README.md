# high-dimensional Iterative Causal Forest (hdiCF)

© 2024 Tiansheng Wang. This work is openly licensed via [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/deed.en)


**Citation**

**Wang T, Pate V, Wyss R, Buse JB, Kosorok MR, Stürmer T. High-dimensional Iterative Causal Forest (hdiCF): a Novel Algorithm for Subgroup Identification in Claims Data. Am J Epidemiol (In Press) March 25, 2024.**

The hdiCF algorithm identifies important subgroups with heterogeneous treatment effects without prior knowledge of treatment-covariate interactions and predefined covariates. The Step 3 is Implementation of [iCF](https://github.com/tianshengwang/iCF).

<img src = images/FIGURE1_hdiCF_3KEYsteps.jpg width=1000>


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

Create analytic cohort and ordinal HD variables in SAS.

***Step 2. Propensity score trimming and HD features preparation***
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


Train <- Train_BENEID_all %>% select(-c("BENE_ID", "IndexDate"))

vars_forest = colnames( Train %>% dplyr::select(-c("Y", "W" ))  )

XYW <- function(TrainDat){
  X <- TrainDat[,vars_forest]
  Y <- as.vector( as.numeric( TrainDat[,"Y"] ) )
  W <- as.vector( as.numeric( TrainDat[,"W"] ) )
  return(list(x=X, y=Y, w=W))
}
#all patients
X <<- XYW(Train)$x
Y <<- XYW(Train)$y
W <<- XYW(Train)$w

ncol(X); nrow(X); length(Y); length(W)
#Z<-Train[,vars_IV]
cf_raw_key.tr <- CF_RAW_key(Train, 1, "hd", hdPctTop=pct_inter) 
#==============================================#==============================================
Y.hat  <<- cf_raw_key.tr$Y.hat                 #
W.hat  <<- cf_raw_key.tr$W.hat  
HTE_P_cf.raw <<- cf_raw_key.tr$HTE_P_cf.raw    # run for overall population or each subgroup
varimp_cf  <- cf_raw_key.tr$varimp_cf          #
#==============================================#==============================================

#W.hat    <- predict(grf::regression_forest(X, W))$predictions

PSplot_allV <- GG_PS(Train, W.hat, "Propensity Score", "PS_allV")

```

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
 <img src = images/VI_allHD.png width=500>


 ```{}
colnames(X[,c(selected_cf.idx)])
X <<- X[,c(selected_cf.idx)]
dat <<- dat[,c("Y", "W", colnames(X)) ]
vars_forest = colnames( dat %>% dplyr::select(-c("Y", "W" ))  )


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
 <img src = images/VI90percentile.png width=500>
 
 ***Step 3. Iplementation of iCF.***
 
 For details of [iCF algorithm: https://github.com/tianshengwang/iCF](https://github.com/tianshengwang/iCF) 
 
 ***To tune the leaf size, use different values for the minimum leaf size (MLS) to grow forests at various depths (D).***
 ```{}
#Specify the decimal position for continuous variables in the subgroup definition.
split_val_round_posi=0
#Define categorical variables with more than two levels:
vars_catover2 <<- find_level_over2(X)  
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

c1_n200_I3_A4_K5_B1000_i100_Tc_L20_V95 <- iCFCV(dat=Train,K=5, treeNo=1000, iterationNo=100, min.split.var=4,
                                                split_val_round_posi=0, P_threshold=0.1, variable_type = "hd", 
                                                hdpct= 0.95, HTE_P_cf.raw = HTE_P_cf.raw) 

#Subgroup decision:
c1_n200_I3_A4_K5_B1000_i100_Tc_L20_V95$selectedSG_ori
```

Propensity score distribution across training and testing sets
```{}
GG_CV_Dx_PS(c1_n200_I3_A4_K5_B1000_i100_Tc_L20_V95, 5)
```
<img src = images/PS.png width=1350>

Inverse probability distribution across training and testing sets
```{}
GG_CV_Dx_iptw(c1_n200_I3_A4_K5_B1000_i100_Tc_L20_V95, 5)
```
<img src = images/IPTW.png width=1350>

transformed outcome distribution across training and testing sets
```{}
GG_CV_Dx_Ystar(c1_n200_I3_A4_K5_B1000_i100_Tc_L20_V95, 5)
```
<img src = images/Ystar.png width=1350>

non-Zero transformed outcome distribution across training and testing sets
```{}
GG_CV_Dx_YstarNo0(c1_n200_I3_A4_K5_B1000_i100_Tc_L20_V95, 5)
```
<img src = images/YstarNo0.png width=1350>

If you have further questions or comments, please contact Dr. Tiansheng Wang: tianwang@unc.edu
