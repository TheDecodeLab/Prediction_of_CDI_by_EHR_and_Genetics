#---- This R code was written by Jiang Li (lijiang1999@gmail.com) for the following manuscript entitled 
#"Prediction of Clostridioides difficile infection using a genetic variant from IL8 and clinical risk factors from EHR"

library(dplyr)
library(caret) # Max Kuhn
library(lattice)
library(ggplot2)
library(RColorBrewer)
library(rattle)
library(rpart)
library(randomcoloR)
library(ROCR)
library(splitstackshape)
library(e1071)
library(cobalt)
library(MatchIt)
library(DMwR)
library(pryr)
library(broom)
library(SDMTools)
library(ROCR)
library(pROC)
library(ResourceSelection)
library(bdpv)

###########################   Section Title #########################

# 1. example of SQL syntax for EHR inquiry
# 2. extraction of genetic dosage file for rs2227306
# 3. propensity score matching for sex and age
# 3.1 get the summary statistics after matchit
# 4. preprocessing the input file for ML
# 5. randomly split cohort into training and testing dataset
# 5.1 calculate MAF for rs2227306 in each subgroups
# 6. simulation of genotypes based on the MAF estimated from MyCODE sample
# 7. chi-square tests
# 8. Oversampling strategies (ROSE vs. SMOTE)
# 9. logistic regression analysis controlled for confounding factors
# 10. density plot after propensity score matching
# 11. create chi-square matrix for input features and create heatmaps
# 12. ML process with customized hyperparameter tunning
# 13. ML process with automatic hyperparameter tunning
# 14. get the summary statistics after matchit
# 15. Delong's test for AUC
# 15. comparing model with or without genetic feature
# 16. create feature importance bar and raddar plots
# 17. AUC/accuracy forest plot
# 18. make figure for f1 score comparison across different oversampling

############################### end ################################


# example of sql syntax for EHR inquiry ---------------------------------------------------


# Convert ICD-9 (787.91) to ICD-10
# The following crosswalk between ICD-9 to ICD-10 is based based on the General Equivalence Mappings (GEMS) information:
#   
# K52.2 - Allergic and dietetic gastroenteritis and colitis 
# K52.89 - Other specified noninfective gastroenteritis and colitis 
# R19.7 - Diarrhea, unspecified 

# SELECT *
#   FROM [PIDB].[dbo].[ENCOUNTERS] 
# INNER JOIN [PIDB].[dbo].[ENCOUNTERS_DX] ON ENCOUNTERS.ENC_ID = ENCOUNTERS_DX.ENC_ID
# WHERE ENC_ADM_DTTM < '08/30/2016 00:00:00.000'
# AND DX_CD in ('787.91', 'R19.7');


# SELECT *
#   FROM [PIDB].[dbo].[LAB_RESULTS]
# WHERE LAB_TKN_DTTM >= '01/01/2000 00:00:00.000' AND LAB_TKN_DTTM < '02/01/2019 00:00:00.000'
# AND LAB_LOINC_CD IN ('13957-6',
#                      '20762-1',
#                      '34713-8',
#                      '54067-4',
#                      '6362-8',
#                      '6362-6',
#                      '83087-7',
#                      '20761-3',
#                      '562-9',
#                      '82197-5',
#                      '61367-9',
#                      '34712-0',
#                      '6364-4',
#                      '6365-1',
#                      '87755-5',
#                      '74822-8',
#                      '46131-9',
#                      '57901-1');



# extraction of genetic dosage file for rs2227306 -------------------------

#input genetic dosage file created by PLINK
genetic <- read.table("GHSf90k_clumped_SNPLIST_recodeA.raw", header = T, stringsAsFactors = F)
head(genetic)
library(splitstackshape)
genetic_split <- cSplit(genetic, "IID", "_")  #
head(genetic_split)  #select IID_2
colnames(genetic_split)[colnames(genetic_split) == 'IID_2'] <- 'PT_ID' 
str(genetic_split)
genetic_split$PT_ID <- as.character(genetic_split$PT_ID)
length(unique(genetic_split$PT_ID)) #85580

phenotype <- read.table("allcov_with_INDEXage_binarycov.txt", header = T, fill = NA, stringsAsFactors = F)
head(phenotype)
phenotype_split <- cSplit(phenotype, "IID", "_")
colnames(phenotype_split)[colnames(phenotype_split) == 'IID_2'] <- 'PT_ID'
str(phenotype_split)
phenotype_split$PT_ID <- as.character(phenotype_split$PT_ID)
length(unique(phenotype_split$PT_ID)) #15304
genetic_phenotype <- merge(genetic_split[, c(6:17, 19)], phenotype_split[, c(2:8, 10:13, 15)], by = "PT_ID")
str(genetic_phenotype)
colnames(genetic_phenotype)[which(names(genetic_phenotype) == "PHENO")] <- "CDIFF"
genetic_phenotype$CDIFF[genetic_phenotype$CDIFF == 2] <- 0 #change noCDIFF into 0
genetic_phenotype_IBD <- merge(genetic_phenotype, ICD9_IBD[, 2:3], by = "PT_ID")
genetic_phenotype_notIBD <- genetic_phenotype[!(genetic_phenotype$PT_ID %in% genetic_phenotype_IBD$PT_ID), ]
genetic_phenotype_notIBD$IBD_risk <- 0
genetic_phenotype <- rbind(genetic_phenotype_IBD, genetic_phenotype_notIBD)
genetic_phenotype %>% group_by(CDIFF, AB_risk, IBD_risk, PPI_risk) %>% summarise(n = n()) %>% ungroup
genetic_phenotype$IBD_risk <- as.integer(genetic_phenotype$IBD_risk)
genetic_phenotype[is.na(genetic_phenotype)] <- 0

save(genetic_phenotype, file = "genetic_phenotype.RData", version = 2)


# propensity score matching for sex and age -----------------------------------------------

load("genetic_phenotype.RData")
head(genetic_phenotype)
merge <- genetic_phenotype
str(merge)
merge$CDIFF <- as.factor(merge$CDIFF)
merge$PT_SEX <- as.factor(merge$PT_SEX)

set.seed(1234)
m.out <- matchit(CDIFF ~ INDEX_AGE + PT_SEX, data = merge, method="nearest", ratio=10) #ratio change from 1 to 5 or 1 to 10
m.out
a <- summary(m.out)
a

bal.tab(m.out)


love.plot(m.out, binary = "std", thresholds = c(m = .1))

summary (m.out,standardize = T)

dim(m.out$match.matrix) #4755*100  #5911*100  #1156*20
match.matrix<-m.out[["match.matrix"]]
index <- unique(as.integer(as.vector(match.matrix)))  #17400 ()  #2726
index
merge.control <- merge[rownames(merge) %in% index, ]
merge.case <- merge[merge$CDIFF == 1, ]
merge.matched <- rbind(merge.case, merge.control)
merge.matched %>%
  dplyr::group_by(CDIFF) %>%
  dplyr::summarise(mean = mean(INDEX_AGE), SD = sd(INDEX_AGE)) %>% 
  ungroup()

est = merge.matched %>% do(tidy(t.test(INDEX_AGE ~ CDIFF, data = .)))
est
est_SEX = merge.matched %>% group_by(PT_SEX) %>% do(tidy(t.test(INDEX_AGE ~ CDIFF, data = .)))
est_SEX

merge %>%
  dplyr::group_by(CDIFF) %>%
  dplyr::summarise(mean = mean(INDEX_AGE), SD = sd(INDEX_AGE)) %>% 
  ungroup()

est = merge %>% do(tidy(t.test(INDEX_AGE ~ CDIFF, data = .)))
est
est_SEX = merge %>% group_by(PT_SEX) %>% do(tidy(t.test(INDEX_AGE ~ CDIFF, data = .)))
est_SEX

merge.matched %>% group_by(CDIFF, PT_SEX) %>% summarise(mean = mean(INDEX_AGE, na.rm = TRUE)) %>% ungroup 
merge.matched %>% group_by(CDIFF, PT_SEX) %>% summarise(n = n()) %>% ungroup 

#Getting strata from nearest neighbor matching
m.data$subclass <- vapply(rownames(m.data), function(x) {
  out <- which(rownames(m1$match.matrix) == x | apply(m1$match.matrix, 1, function(y) x %in% y))
  if (length(out) == 0) out <- NA_integer_
  out
}, integer(1L))



# get the summary statistics after matchit --------------------------------
  
  
load("merge_training.RData")  #n = 41786
load("merge_testing.RData")   #n = 17907
load("merge_smotetrain.RData")  #n = 13252


temp <- rbind(training, testing)
colnames(temp)
# colnames(temp)
# [1] "CDIFF"           "X4.74607055_T"   "X6.32177263_G"   "AB_risk"         "chemo_risk"      "PPI_risk"        "steroid_risk"   
# [8] "TNF_risk"        "transplant_risk" "HIV"             "T2DM"            "PT_SEX.Female"   "INDEX_AGE"       "IBD_risk" 
temp %>%
  #distinct(PT_ID, .keep_all = T) %>%
  dplyr::group_by(CDIFF) %>%
  dplyr::summarise(n = n()) %>% 
  ungroup()

temp %>%
  #distinct(PT_ID, .keep_all = T) %>%
  dplyr::count(CDIFF, X4.74607055_T) %>%
  dplyr::group_by(CDIFF) %>%
  dplyr::mutate(prop = prop.table(n)) %>%
  ungroup()
summary(table(temp$CDIFF, temp$X4.74607055_T))


temp %>%
  dplyr::count(CDIFF, AB_risk) %>%
  dplyr::group_by(CDIFF) %>%          # now required with changes to dplyr::count()
  dplyr::mutate(prop = prop.table(n)) %>%
  ungroup()
summary(table(temp$CDIFF, temp$AB_risk))
temp %>%
  dplyr::count(CDIFF, AB) %>%
  dplyr::group_by(CDIFF) %>%          # now required with changes to dplyr::count()
  dplyr::mutate(prop = prop.table(n)) %>%
  ungroup()
summary(table(temp$CDIFF, temp$AB))

temp %>%
  dplyr::count(CDIFF, chemo_risk) %>%
  dplyr::group_by(CDIFF) %>%          # now required with changes to dplyr::count()
  dplyr::mutate(prop = prop.table(n)) %>%
  ungroup()
summary(table(temp$CDIFF, temp$chemo_risk))

temp %>%
  dplyr::count(CDIFF, TNF_risk) %>%
  dplyr::group_by(CDIFF) %>%          # now required with changes to dplyr::count()
  dplyr::mutate(prop = prop.table(n)) %>%
  ungroup()
summary(table(temp$CDIFF, temp$TNF_risk))
fisher.test(temp$CDIFF,temp$TNF_risk)

temp %>%
  dplyr::count(CDIFF, transplant_risk) %>%
  dplyr::group_by(CDIFF) %>%          # now required with changes to dplyr::count()
  dplyr::mutate(prop = prop.table(n)) %>%
  ungroup()
summary(table(temp$CDIFF, temp$transplant_risk))
fisher.test(temp$CDIFF,temp$transplant_risk)

# preprocessing the input file for ML -------------------------------------

head(merge1)
merge1 %>% group_by(CDIFF) %>% tally()
merge1 %>% group_by(CDIFF, AB, IBD, PPI) %>% summarise(n = n()) %>% ungroup
# CDIFF    AB   IBD   PPI     n
# <int> <int> <int> <int> <int>
#   1     0     0     0     0 45063
# 2     0     0     0     1  3991
# 3     0     0     1     0   433
# 4     0     0     1     1    81
# 5     0     1     0     0  6788
# 6     0     1     0     1  1642
# 7     0     1     1     0    70
# 8     0     1     1     1    30
# 9     1     0     0     0  2054
# 10     1     0     0     1   537
# 11     1     0     1     0   176
# 12     1     0     1     1    51
# 13     1     1     0     0  1180
# 14     1     1     0     1   760
# 15     1     1     1     0    67
# 16     1     1     1     1    48

pp <- c('center', 'scale')
str(merge1)
preProcValues <- preProcess(merge1[, c(1, 11)], method = pp) #create preProcessing object for AGE (No 11 column)
merge1_impute <- predict(preProcValues, merge1) #apply preProcessing to dataset
head(merge1_impute)

# randomly split cohort into training and testing datasets ----------------



set.seed(42)
train <- createDataPartition(merge1_impute$CDIFF, #target variable vector
                             p = .70,  # % of data for training
                             list = FALSE, #should result be a list (T/F)
                             times = 1)  #Num of partitions to create
train
training <- merge1_impute[train,]
training
testing <- merge1_impute[-train,]
testing

dummies <- dummyVars(CDIFF ~ ., data = training[, 2:ncol(training)])
training_new <- data.frame(predict(dummies, newdata = training[, 2:ncol(training)]))
names(training_new)
str(training_new)
training <- cbind(training$CDIFF, training_new)
names(training)[1] <- "CDIFF"
head(training)
training$CDIFF <- as.factor(training$CDIFF) #target variable has to be as.factor but not INT

levels(training$CDIFF) <- c("ACDIFF", "CDIFF")

training$transplant <- as.factor(training$transplant)
training$AB <- as.factor(training$AB)
training$CHEMO <- as.factor(training$CHEMO)
training$PPI <- as.factor(training$PPI)
training$steroid <- as.factor(training$steroid)
training$TNF <- as.factor(training$TNF)
training$HIV <- as.factor(training$HIV)
training$T2D <- as.factor(training$T2D)
training$IBD <- as.factor(training$IBD)
training$PT_SEXFemale <- as.factor(training$PT_SEXFemale)

summary(training)

dummies <- dummyVars(CDIFF ~ ., data = testing[, 2:ncol(testing)])
testing_new <- data.frame(predict(dummies, newdata = testing[, 2:ncol(testing)]))
names(testing_new)
str(testing_new)
testing <- cbind(testing$CDIFF, testing_new)
names(testing)[1] <- "CDIFF"
head(testing)
summary(testing)
testing$CDIFF <- as.factor(testing$CDIFF) #target variable has to be as.factor but not INT

levels(testing$CDIFF) <- c("ACDIFF", "CDIFF")

testing$transplant <- as.factor(testing$transplant)
testing$AB <- as.factor(testing$AB)
testing$CHEMO <- as.factor(testing$CHEMO)
testing$PPI <- as.factor(testing$PPI)
testing$steroid <- as.factor(testing$steroid)
testing$TNF <- as.factor(testing$TNF)
testing$HIV <- as.factor(testing$HIV)
testing$T2D <- as.factor(testing$T2D)
testing$IBD <- as.factor(testing$IBD)
testing$PT_SEXFemale <- as.factor(testing$PT_SEXFemale)

summary(testing)

# calculate MAF for rs2227306 in each subgroups ---------------------------



training_testing_young_female <- training_testing_young[training_testing_young$PT_SEXFemale == 1, ] #5214
training_testing_young_male <- training_testing_young[training_testing_young$PT_SEXFemale == 0, ] #1803

training_testing_old_female <- training_testing_old[training_testing_old$PT_SEXFemale == 1, ] #4412
training_testing_old_male <- training_testing_old[training_testing_old$PT_SEXFemale == 0, ] #3875

a <- table(training_testing_young_female$CDIFF, training_testing_young_female$X4.74607055_T)
# 0    1    2
# ACDIFF 1769 2303  856  
# CDIFF    86  151   49

a[1,3]/(a[1,1] + a[1,2] + a[1,3]) #0.1737013
a[1,2]/(a[1,1] + a[1,2] + a[1,3]) #0.4673295
a[1,1]/(a[1,1] + a[1,2] + a[1,3]) #0.3589692

#MAF 
sqrt(a[1,3]/(a[1,1] + a[1,2] + a[1,3])) #0.4167749 theoretical

(a[1,2] + 2*a[1,3])/(2*(a[1,1] + a[1,2] + a[1,3])) #0.4073661 real

a[2,3]/(a[2,1] + a[2,2] + a[2,3]) #0.1713287
a[2,2]/(a[2,1] + a[2,2] + a[2,3]) #0.527972
a[2,1]/(a[2,1] + a[2,2] + a[2,3]) #0.3006993

sqrt(a[2,3]/(a[2,1] + a[2,2] + a[2,3])) #0.4139187 theoretical

(a[2,2] + 2*a[2,3])/(2*(a[2,1] + a[2,2] + a[2,3])) #0.4353147 real


b <- table(training_testing_young_male$CDIFF, training_testing_young_male$X4.74607055_T)
# 0   1   2
# ACDIFF 573 807 266
# CDIFF   43  86  28

(b[1,2] + 2*b[1,3])/(2*(b[1,1] + b[1,2] + b[1,3])) #0.4067436 rebl
(b[2,2] + 2*b[2,3])/(2*(b[2,1] + b[2,2] + b[2,3])) #0.4522293 rebl

c <- table(training_testing_old_female$CDIFF, training_testing_old_female$X4.74607055_T)
# 0    1    2
# ACDIFF 1380 1910  731
# CDIFF   117  196   78

(c[1,2] + 2*c[1,3])/(2*(c[1,1] + c[1,2] + c[1,3])) #0.4192987 recl
(c[2,2] + 2*c[2,3])/(2*(c[2,1] + c[2,2] + c[2,3])) #0.4501279 recl

d <- table(training_testing_old_male$CDIFF, training_testing_old_male$X4.74607055_T)
# 0    1    2
# ACDIFF 1231 1666  656
# CDIFF   113  150   59

(d[1,2] + 2*d[1,3])/(2*(d[1,1] + d[1,2] + d[1,3])) #0.4190825 redl
(d[2,2] + 2*d[2,3])/(2*(d[2,1] + d[2,2] + d[2,3])) #0.4161491 redl

# AGE_binary PT_SEXFemale CDIFF      n     MAF
# <dbl> <fct>        <fct>  <int>
# 1          0 0            ACDIFF 12757   0.4067436
# 2          0 0            CDIFF    527   0.4522293
# 3          0 1            ACDIFF 16072   0.4073661
# 4          0 1            CDIFF    688   0.4353147
# 5          1 0            ACDIFF 12066   0.4190825
# 6          1 0            CDIFF   1464   0.4161491
# 7          1 1            ACDIFF 14043   0.4192987
# 8          1 1            CDIFF   2076   0.4501279

table(training_testing_young$CDIFF, training_testing_young$X4.74607055_T)

# 0    1    2
# ACDIFF 2342 3110 1122
# CDIFF   129  237   77

#CC in ACDIFF: 
2342/(2342+3110+1122) #0.3562519
3110/(2342+3110+1122) #0.4730758
1122/(2342+3110+1122) #0.1706723
#CC in CDIFF:
129/(129+237+77) #0.2911964
237/(129+237+77) #0.5349887
77/(129+237+77)  #0.1738149

#%ACDIFF
(2342+3110+1122)/(2342+3110+1122+129+237+77) #0.9368676
#%CDIFF
(129+237+77)/(2342+3110+1122+129+237+77) #0.06313239
table(training_testing_old$CDIFF, training_testing_old$X4.74607055_T)

# 0    1    2
# ACDIFF 2611 3576 1387
# CDIFF   230  346  137

#CC in ACDIFF:
2611/(2611+3576+1387) #0.344732
3576/(2611+3576+1387) #0.4721415
1387/(2611+3576+1387) #0.1831265

#CC in CDIFF:
230/(230+346+137) #0.3225806
346/(230+346+137) #0.4852735
137/(230+346+137) #0.1921459

#%ACDIFF
(2611+3576+1387)/(2611+3576+1387+230+346+137) #0.9139616
#%CDIFF
(230+346+137)/(2611+3576+1387+230+346+137) #0.08603837


# simulation of genotypes based on the MAF estimated from MyCODE sample -------------------------------------------------



# simulation strategy I

length(merge.matched$CDIFF[merge.matched$CDIFF == 1])  
length(merge.matched$CDIFF[merge.matched$CDIFF == 0])  

Nc <- 4755  # number of cases
Nn <- 54938  # number of controls

4755+54938 #59693

status <- c(rep(1, Nc), rep(0, Nn))  # case/non-case status
genIL8 <- matrix(nrow=Nc+Nn, ncol=2)

#matched population controls
#X4.74607055_T(rs2227306)
#MAF was set at 0.177 in control and 0.185 in case
f0<-sqrt(0.177);f1<-sqrt(0.185);o1<-f1*(1-f0)/f0/(1-f1); o1 #1.039217

genIL8[1:Nc, 1] <- rbinom(Nc, 1, f1)+rbinom(Nc, 1, f1);

genIL8[(Nc+1):(Nc+Nn), 1] <- rbinom(Nn, 1, f0)+rbinom(Nn, 1, f0)

genIL8[, 2] <- status

colnames(genIL8) <- c("X4.74607055_T", "CDIFF1")

merge.matched.genotyped <- cbind(merge[order(merge$CDIFF, decreasing = TRUE),], genIL8)

save(merge.matched.genotyped, file = "merge.notmatched.genotypesimulated.RData", version = 2)

#Simulation strategy II
#need to consider age and sex

training$split <- 1
testing$split <- 2
training_testing <- rbind(training, testing) #59693
colnames(training_testing)

training_testing$AGE_binary <- NA
training_testing[training_testing$AGE >= 0, ]$AGE_binary <- 1
training_testing[training_testing$AGE < 0, ]$AGE_binary <- 0

summary <- training_testing %>% 
  group_by(AGE_binary, PT_SEXFemale, CDIFF) %>% 
  summarise(n = n())
summary
# AGE_binary PT_SEXFemale CDIFF      n     MAF(estimated from MyCode sample)
# <dbl> <fct>        <fct>  <int>
# 1          0 0            ACDIFF 12757   0.4067436
# 2          0 0            CDIFF    527   0.4522293
# 3          0 1            ACDIFF 16072   0.4073661
# 4          0 1            CDIFF    688   0.4353147
# 5          1 0            ACDIFF 12066   0.4190825
# 6          1 0            CDIFF   1464   0.4161491
# 7          1 1            ACDIFF 14043   0.4192987
# 8          1 1            CDIFF   2076   0.4501279

N <- summary$n
genIL8 <- matrix(nrow=Nc+Nn, ncol=2)
f <- c(0.4067436, 0.4522293, 0.4073661, 0.435314, 0.4190825, 0.4161491, 0.4192987, 0.4501279)
summary$MAF <- f
save(summary, file = "summary_simulation_conditional_table.RData", version = 2)

f1 <- summary$MAF[1]
genIL8[1:N[1], 1] <- rbinom(N[1], 1, f1)+rbinom(N[1], 1, f1);

f2 <- summary$MAF[2]
genIL8[(N[1]+1):(N[1]+N[2]), 1] <- rbinom(N[2], 1, f2)+rbinom(N[2], 1, f2)

f3 <- summary$MAF[3]
genIL8[(N[1]+N[2]+1):(N[1]+N[2]+N[3]), 1] <- rbinom(N[3], 1, f3)+rbinom(N[3], 1, f3)

f4 <- summary$MAF[4]
genIL8[(N[1]+N[2]+N[3]+1):(N[1]+N[2]+N[3]+N[4]), 1] <- rbinom(N[4], 1, f4)+rbinom(N[4], 1, f4)

f5 <- summary$MAF[5]
genIL8[(N[1]+N[2]+N[3]+N[4]+1):(N[1]+N[2]+N[3]+N[4]+N[5]), 1] <- rbinom(N[5], 1, f5)+rbinom(N[5], 1, f5)

f6 <- summary$MAF[6]
genIL8[(N[1]+N[2]+N[3]+N[4]+N[5]+1):(N[1]+N[2]+N[3]+N[4]+N[5]+N[6]), 1] <- rbinom(N[6], 1, f6)+rbinom(N[6], 1, f6)

f7 <- summary$MAF[7]
genIL8[(N[1]+N[2]+N[3]+N[4]+N[5]+N[6]+1):(N[1]+N[2]+N[3]+N[4]+N[5]+N[6]+N[7]), 1] <- rbinom(N[7], 1, f7)+rbinom(N[7], 1, f7)

f8 <- summary$MAF[8]
genIL8[(N[1]+N[2]+N[3]+N[4]+N[5]+N[6]+N[7]+1):(N[1]+N[2]+N[3]+N[4]+N[5]+N[6]+N[7]+N[8]), 1] <- rbinom(N[8], 1, f8)+rbinom(N[8], 1, f8)

sum(is.na(genIL8[, 1])) 

training_testing_new <- training_testing[with(training_testing, order(AGE_binary, PT_SEXFemale, CDIFF)),]
training_testing_new$X4.74607055_T <- genIL8[, 1]

training_testing_young <- training_testing[training_testing$AGE_binary == 0, ]
training_testing_young_male <- training_testing_young[training_testing_young$PT_SEXFemale == 0, ]
training_testing_young_male_ACDIFF <- training_testing_young_male[training_testing_young_male$CDIFF == "ACDIFF", ]
training_testing_young_male_CDIFF <- training_testing_young_male[training_testing_young_male$CDIFF == "CDIFF", ]
training_testing_young_female <- training_testing_young[training_testing_young$PT_SEXFemale == 1, ]
training_testing_young_female_ACDIFF <- training_testing_young_female[training_testing_young_female$CDIFF == "ACDIFF", ]
training_testing_young_female_CDIFF <- training_testing_young_female[training_testing_young_female$CDIFF == "CDIFF", ]

training_testing_old <- training_testing[training_testing$AGE_binary == 1, ]
training_testing_old_male <- training_testing_old[training_testing_old$PT_SEXFemale == 0, ]
training_testing_old_male_ACDIFF <- training_testing_old_male[training_testing_old_male$CDIFF == "ACDIFF", ]
training_testing_old_male_CDIFF <- training_testing_old_male[training_testing_old_male$CDIFF == "CDIFF", ]
training_testing_old_female <- training_testing_old[training_testing_old$PT_SEXFemale == 1, ]
training_testing_old_female_ACDIFF <- training_testing_old_female[training_testing_old_female$CDIFF == "ACDIFF", ]
training_testing_old_female_CDIFF <- training_testing_old_female[training_testing_old_female$CDIFF == "CDIFF", ]

training_testing_new <-do.call("rbind", list(training_testing_young_male_ACDIFF, training_testing_young_male_CDIFF, training_testing_young_female_ACDIFF, training_testing_young_female_CDIFF, training_testing_old_male_ACDIFF, training_testing_old_male_CDIFF, training_testing_old_female_ACDIFF, training_testing_old_female_CDIFF))
training_testing_new$X4.74607055_T <- genIL8[, 1]

summary_new <- training_testing_new %>% 
  group_by(AGE_binary, PT_SEXFemale, CDIFF, X4.74607055_T) %>% 
  summarise(n = n())
summary_new
# AGE_binary PT_SEXFemale CDIFF  X4.74607055_T     n
# <dbl> <fct>        <fct>          <int> <int>
#   1          0 0            ACDIFF             0  4370
# 2          0 0            ACDIFF             1  6285
# 3          0 0            ACDIFF             2  2102
# 4          0 0            CDIFF              0   147
# 5          0 0            CDIFF              1   271
# 6          0 0            CDIFF              2   109
# 7          0 1            ACDIFF             0  5641
# 8          0 1            ACDIFF             1  7690
# 9          0 1            ACDIFF             2  2741
# 10          0 1            CDIFF              0   228

E1 <- (summary_new$n[2]+2*summary_new$n[3])/(2*(summary_new$n[1] + summary_new$n[2] + summary_new$n[3])) 
E2 <- (summary_new$n[5]+2*summary_new$n[6])/(2*(summary_new$n[4] + summary_new$n[5] + summary_new$n[6])) 
E3 <- (summary_new$n[8]+2*summary_new$n[9])/(2*(summary_new$n[7] + summary_new$n[8] + summary_new$n[9])) 
E4 <- (summary_new$n[11]+2*summary_new$n[12])/(2*(summary_new$n[10] + summary_new$n[11] + summary_new$n[12]))
E5 <- (summary_new$n[14]+2*summary_new$n[15])/(2*(summary_new$n[13] + summary_new$n[14] + summary_new$n[15]))
E6 <- (summary_new$n[17]+2*summary_new$n[18])/(2*(summary_new$n[16] + summary_new$n[17] + summary_new$n[18]))
E7 <- (summary_new$n[20]+2*summary_new$n[21])/(2*(summary_new$n[19] + summary_new$n[20] + summary_new$n[21]))
E8 <- (summary_new$n[23]+2*summary_new$n[24])/(2*(summary_new$n[22] + summary_new$n[23] + summary_new$n[24]))

summary <- training_testing %>% 
  group_by(AGE_binary, PT_SEXFemale, CDIFF, X4.74607055_T) %>% 
  summarise(n = n())
summary

# AGE_binary PT_SEXFemale CDIFF  X4.74607055_T     n
# <dbl> <fct>        <fct>          <dbl> <int>
#   1          0 0            ACDIFF             0   573
# 2          0 0            ACDIFF             1   807
# 3          0 0            ACDIFF             2   266
# 4          0 0            CDIFF              0    43
# 5          0 0            CDIFF              1    86
# 6          0 0            CDIFF              2    28
# 7          0 1            ACDIFF             0  1769
# 8          0 1            ACDIFF             1  2303
# 9          0 1            ACDIFF             2   856
# 10          0 1            CDIFF              0    86

E1 <- (summary$n[2]+2*summary$n[3])/(2*(summary$n[1] + summary$n[2] + summary$n[3])) 
E2 <- (summary$n[5]+2*summary$n[6])/(2*(summary$n[4] + summary$n[5] + summary$n[6])) 
E3 <- (summary$n[8]+2*summary$n[9])/(2*(summary$n[7] + summary$n[8] + summary$n[9])) 
E4 <- (summary$n[11]+2*summary$n[12])/(2*(summary$n[10] + summary$n[11] + summary$n[12]))
E5 <- (summary$n[14]+2*summary$n[15])/(2*(summary$n[13] + summary$n[14] + summary$n[15]))
E6 <- (summary$n[17]+2*summary$n[18])/(2*(summary$n[16] + summary$n[17] + summary$n[18]))
E7 <- (summary$n[20]+2*summary$n[21])/(2*(summary$n[19] + summary$n[20] + summary$n[21]))
E8 <- (summary$n[23]+2*summary$n[24])/(2*(summary$n[22] + summary$n[23] + summary$n[24]))


E <- c(E1, E2, E3, E4, E5, E6, E7, E8)
#[1] 0.4067436 0.4522293 0.4073661 0.4353147 0.4190825 0.4161491 0.4192987 0.4501279 #genetic sample
#[1] 0.4111076 0.4639469 0.4097810 0.4316860 0.4174126 0.3924180 0.4233070 0.4381021 #nongenetic sample newmethod simulated

summary$MAF_estimated <- E

# chi-square tests --------------------------------------------------------



table(training_testing$CDIFF, training_testing$X4.74607055_T)  #oldmethod of simulation
# 0     1     2
# ACDIFF 18566 26635  9737
# CDIFF   1499  2350   906

chisq.test(training_testing$CDIFF, training_testing$X4.74607055_T)
#X-squared = 11.845, df = 2, p-value = 0.002679

chisq.test(training_testing[training_testing$AGE_binary == 0 & training_testing$PT_SEXFemale == 0, ]$CDIFF, training_testing[training_testing$AGE_binary == 0 & training_testing$PT_SEXFemale == 0, ]$X4.74607055_T)
#X-squared = 2.8171, df = 2, p-value = 0.2445

chisq.test(training_testing[training_testing$AGE_binary == 0 & training_testing$PT_SEXFemale == 1, ]$CDIFF, training_testing[training_testing$AGE_binary == 0 & training_testing$PT_SEXFemale == 1, ]$X4.74607055_T)
#X-squared = 2.7461, df = 2, p-value = 0.2533

chisq.test(training_testing[training_testing$AGE_binary == 1 & training_testing$PT_SEXFemale == 0, ]$CDIFF, training_testing[training_testing$AGE_binary == 1 & training_testing$PT_SEXFemale == 0, ]$X4.74607055_T)
#X-squared = 1.2516, df = 2, p-value = 0.5348

chisq.test(training_testing[training_testing$AGE_binary == 1 & training_testing$PT_SEXFemale == 1, ]$CDIFF, training_testing[training_testing$AGE_binary == 1 & training_testing$PT_SEXFemale == 1, ]$X4.74607055_T)
#X-squared = 8.8022, df = 2, p-value = 0.01226



# Oversampling strategies (ROSE vs. SMOTE) -------------------------------------------------


##SMOTE
smotetrain <- SMOTE(CDIFF ~ ., data = training, perc.over = 100, k = 5) # do not use training$CDIFF
summary(smotetrain)

library("ROSE")

##ROSE
### Calculate the total number of required cases in the over-sampled dataset
data <- training
print(table(data$CDIFF))
# ACDIFF  CDIFF 
# 9931    782
n_new <- 9931 / (1 - 0.50)

oversampling_result <- ovun.sample(formula = CDIFF ~ ., data = data,
                                   method = "over", N = n_new, seed = 2018)

# Verify the Class-balance of the over-sampled dataset
oversampled_credit <- oversampling_result$data
prop.table(table(oversampled_credit$CDIFF))

table(oversampling_result$data$CDIFF)
# ACDIFF  CDIFF 
# 9931   9931
rosetrain <- oversampling_result$data
save(rosetrain, file = "genetic_rosetrain.RData", version = 2)


# logistic regression analysis controlled for confounding factors  --------



training_testing_young <- training_testing[training_testing$INDEX_AGE < 0, ] #7017
training_testing_old <- training_testing[training_testing$INDEX_AGE >= 0, ] #8287

glm.fit <- glm(CDIFF ~ X4.74607055_T + as.factor(PT_SEXFemale), data = training_testing_young, family = binomial)
summary(glm.fit)

# Call:
#   glm(formula = CDIFF ~ X4.74607055_T + as.factor(PT_SEXFemale), 
#       family = binomial, data = training_testing_young)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -0.4604  -0.3628  -0.3393  -0.3173   2.4553  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)              -2.46704    0.10351 -23.834  < 2e-16 ***
#   X4.74607055_T             0.13811    0.06977   1.980   0.0478 *  
#   as.factor(PT_SEXFemale)1 -0.49684    0.10336  -4.807 1.53e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 3305.0  on 7016  degrees of freedom
# Residual deviance: 3279.1  on 7014  degrees of freedom
# AIC: 3285.1
# 
# Number of Fisher Scoring iterations: 5

glm.fit <- glm(CDIFF ~ X4.74607055_T + as.factor(PT_SEXFemale), data = training_testing_old, family = binomial)
summary(glm.fit)

# Call:
#   glm(formula = CDIFF ~ X4.74607055_T + as.factor(PT_SEXFemale), 
#       family = binomial, data = training_testing_old)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -0.4455  -0.4326  -0.4201  -0.4184   2.2521  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)              -2.45349    0.07519  -32.63   <2e-16 ***
#   X4.74607055_T             0.06171    0.05512    1.12    0.263    
# as.factor(PT_SEXFemale)1  0.07005    0.07870    0.89    0.373    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 4860.7  on 8286  degrees of freedom
# Residual deviance: 4858.7  on 8284  degrees of freedom
# AIC: 4864.7
# 
# Number of Fisher Scoring iterations: 5


training_testing_female <- training_testing[training_testing$PT_SEXFemale == 1, ] #9626
training_testing_male <- training_testing[training_testing$PT_SEXFemale == 0, ] #5678
glm.fit <- glm(CDIFF ~ X4.74607055_T + INDEX_AGE, data = training_testing_female, family = binomial)
summary(glm.fit)

# Call:
#   glm(formula = CDIFF ~ X4.74607055_T + INDEX_AGE, family = binomial, 
#       data = training_testing_female)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -0.4812  -0.4071  -0.3728  -0.3400   2.5661  
# 
# Coefficients:This 
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   -2.66515    0.06372 -41.827  < 2e-16 ***
#   X4.74607055_T  0.11905    0.05623   2.117   0.0342 *  
#   INDEX_AGE      0.20499    0.04053   5.058 4.25e-07 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 4899.5  on 9625  degrees of freedom
# Residual deviance: 4868.7  on 9623  degrees of freedom
# AIC: 4874.7
# 
# Number of Fisher Scoring iterations: 5


glm.fit <- glm(CDIFF ~ X4.74607055_T + INDEX_AGE, data = training_testing_male, family = binomial)
summary(glm.fit)


# Call:
#   glm(formula = CDIFF ~ X4.74607055_T + INDEX_AGE, family = binomial, 
#       data = training_testing_male)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -0.4822  -0.4278  -0.4155  -0.4042   2.2832  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   -2.40868    0.07580 -31.776   <2e-16 ***
#   X4.74607055_T  0.05386    0.06775   0.795    0.427    
# INDEX_AGE     -0.07783    0.05004  -1.555    0.120    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 3285.2  on 5677  degrees of freedom
# Residual deviance: 3282.2  on 5675  degrees of freedom
# AIC: 3288.2
# 
# Number of Fisher Scoring iterations: 5


# density plot after propensity score matching ----------------------------



head(training)
#after propensity score matching
library(ggplot2)

tiff("Rplot_densityplot_training_genetic_indexage_ratio5.tiff", units="in", width=9, height=8, res=300) 
p <- ggplot(training, aes(x = INDEX_AGE, colour = (CDIFF == "CDIFF")), alpha = 0.9)  # or AGE for data without genetic information
p + ggtitle("The training dataset of the cohort w/ genetic data w/propensity score match(1:5)") + geom_density() + theme(text = element_text(size=15)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            panel.background = element_blank(), axis.line = element_line(colour = "black")) # Your example
dev.off()


# create chi-square matrix for input features and create heatmaps -----------------------------


head(training)
training_select <- training[, c(1, 2, 14:24)]

training_select$INDEX_AGE_Binary <- ifelse(training_select$INDEX_AGE >= 0, 1, 0)
training_select <- training_select[, c(1:11, 13:14)]

#select this function for all cells without missingness
chisqmatrix <- function(x) {
  names = colnames(x);  num = length(names)
  m = matrix(nrow=num,ncol=num,dimnames=list(names,names))
  for (i in 1:(num-1)) {
    for (j in (i+1):num) {
      #browser()
      m[j,i] = chisq.test(x[, i, drop = TRUE],x[, j, drop = TRUE])$p.value
    }
  }
  return (m)
}

#select this function to get all cells filled to avoid error message
chisqmatrix <- function(x) {
  names = colnames(x);  num = length(names)
  m = matrix(nrow=num,ncol=num,dimnames=list(names,names))
  for (i in 1:num) {
    for (j in 1:num) {
      #browser()
      m[j,i] = chisq.test(x[, i, drop = TRUE],x[, j, drop = TRUE])$p.value
    }
  }
  return (m)
}


colnames(training_select)
# [1] "CDIFF"            "X4.74607055_T"    "AB_risk"          "chemo_risk"       "PPI_risk"         "steroid_risk"    
# [7] "TNF_risk"         "transplant_risk"  "HIV"              "T2DM"             "PT_SEXFemale"     "IBD_risk"        
# [13] "INDEX_AGE_Binary"

# [1] "CDIFF"            "transplant"       "AB"               "CHEMO"            "PPI"              "steroid"         
# [7] "TNF"              "HIV"              "T2D"              "IBD"              "PT_SEXFemale"     "X4.74607055_T"   
# [13] "INDEX_AGE_Binary"

#optional
colnames(training_select) <- c("CDIFF", "rs2227306(IL8)", "Antibiotics", "Chemotherapy", "Proton Pump Inhibitors", "Corticosteroid", "Anti-TNF Med", "Transplant Med", "HIV", "Type 2 Diabetic", "Sex(Female)", "Inflammatory Bowel Dis", "Age at index(binary)")
colnames(training_select) <- c("CDIFF", "Transplant Med", "Antibiotics", "Chemotherapy", "Proton Pump Inhibitors", "Corticosteroid", "Anti-TNF Med", "HIV", "Type 2 Diabetic", "Inflammatory Bowel Dis", "Sex(Female)", "rs2227306(IL8)", "Age at index(binary)")

mat <- chisqmatrix(training_select)
mat
#mat[-1, -ncol(mat)]

# filled the NA with 1
mat[is.na(mat)] <- 1
mat

mat_transformed <- function(x) {-log(x, base = 10)}
mat_final <- apply(mat, 2, mat_transformed)
mat_final[mat_final == "Inf"] <- 300
mat_final

mat_final[mat_final >= 25] <- 25
mat_final
colnames(mat_final)

heatmap(mat_final, scale = "none")

#heatmap2
library(gplots)
my_palette <- colorRampPalette(c("green", "blue"))(n = 100)
heatmap.2(as.matrix(mat_final), col=my_palette, 
          breaks=colors, density.info="none", trace="none", 
          dendrogram=c("col"), symm=F,symkey=F,symbreaks=T, scale="none")

tiff("Rplot_propensitymatching_training_chisquare_pvalue_heatmap.tiff", units="in", width=9, height=8, res=300)  
heatmap.2(as.matrix(mat_final),col=greenred(75),
          dendrogram=c("both"),trace="none",
          Rowv=TRUE,margins = c(8,16),
          cexRow=1.0,cexCol=1.0)
dev.off()

tiff("Rplot_matching_training_chisquare_pvalue_heatmap_genetic.tiff", units="in", width=9, height=8, res=300)  
heatmap.2(as.matrix(mat_final),col=greenred(75),
          dendrogram=c("none"),trace="none",
          Rowv=TRUE,margins = c(8,16),
          cexRow=1.0,cexCol=1.0, srtCol = 45,
          main= expression(paste("p value for 2X2 Chi-square tests \n (The cohort w/ genetic data)")),
          key=TRUE, 
          key.par = list(mar= c(8, 0, 0, 10)),
          symm=F,
          symkey=T,
          symbreaks=T)
dev.off()


# ML process with customized hyperparameter tunning -----------------------


training <- training[, c(1, 3:13)]
testing <- testing[, c(1, 3:13)]
smotetrain <- smotetrain[, c(1, 3:13)]

xgbGrid <- expand.grid(nrounds = c(150), 
                       max_depth = c(2, 3, 4),
                       eta = 0.40,
                       rate_drop = 0.10,
                       skip_drop = 0.10,
                       colsample_bytree = c(0.6, 0.7, 0.8),
                       min_child_weight = c(1, 2, 3),
                       subsample = 1,
                       gamma = c(0, 0.05, 0.1))

gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3, 4, 5),
                        n.trees = c(100, 150, 200), #(0:50)*50
                        shrinkage = seq(.005, .2,.005),
                        n.minobsinnode = 10) # you can also put something like c(5, 10, 15, 20)

#with repeats
ctrl <- trainControl(method="repeatedcv", 
                     number = 5,
                     classProbs = T, 
                     summaryFunction = twoClassSummary,
                     verboseIter = TRUE,
                     savePredictions = TRUE,
                     repeats = 10,
                     allowParallel = T)
#without repeats
ctrl <- trainControl(method="repeatedcv", 
                     number = 5,
                     classProbs = T, 
                     summaryFunction = twoClassSummary,
                     verboseIter = TRUE,
                     savePredictions = TRUE,
                     allowParallel = T)

set.seed(1234)
modeled <- train(CDIFF ~., data = smotetrain,
                 method = "gbm", 
                 trControl = ctrl,
                 #metric = "ROC", 
                 #preProc = c("center", "scale"),
                 tuneGrid = gbmGrid)

set.seed(1234)
modeled <- train(CDIFF ~., data = smotetrain,
                             method = "xgbDART", 
                             trControl = ctrl,
                             #metric = "ROC", 
                             #preProc = c("center", "scale"),
                             tuneGrid = xgbGrid)

registerDoMC(cores=4)

modeled <- train(CDIFF ~., data = smotetrain, method = "xgbDART", trControl=ctrl)
modeled <- train(CDIFF ~., data = smotetrain, method = "gbm", trControl=ctrl) 

str(smotetrain)
importance <- varImp(modeled, scale=FALSE)
# summarize importance
print(importance) 

plot(importance, top = 12)


TrainP <- predict(modeled, newdata=smotetrain)  #change from smotetrain to training
TestP <- predict(modeled, newdata=testing)
traincm <- confusionMatrix(data=TrainP, smotetrain$CDIFF, positive = "CDIFF")  #change from smotetrain$CDIFF to training$CDIFF
testcm <- confusionMatrix(data=TestP, testing$CDIFF, positive = "CDIFF")
traincm$overall

testcm$overall

traincm$byClass

testcm$byClass

x <- c(traincm$overall[1], traincm$overall[3], traincm$overall[4], traincm$byClass[1], traincm$byClass[2], 
       traincm$byClass[3], traincm$byClass[4], traincm$byClass[5], traincm$byClass[6], traincm$byClass[7],
       traincm$byClass[8], traincm$byClass[9], traincm$byClass[10], traincm$byClass[11],
       testcm$overall[1], testcm$overall[3], testcm$overall[4], testcm$byClass[1], testcm$byClass[2], 
       testcm$byClass[3], testcm$byClass[4], testcm$byClass[5], testcm$byClass[6], testcm$byClass[7],
       testcm$byClass[8], testcm$byClass[9], testcm$byClass[10], testcm$byClass[11])
x <- t(data.frame(x))
x

probs <- predict(modeled, testing, type = "prob")[,2]
pred1 <- prediction(probs, testing$CDIFF)
perf <- performance(pred1, "tpr", "fpr")
perf1 <- performance(pred1, "auc")

#perf1@x.values[[1]]
perf1@y.values
plot(perf@x.values[[1]], perf@y.values[[1]], type = 'l', ylab = perf@y.name, xlab = perf@x.name, col = '1', lwd = 2)
lines(perf@x.values[[1]], perf@y.values[[1]], col = '2', lty = 2)
abline(a = 0, b = 1)

# ML process with automatic hyperparameter tunning ------------------------


modeloutput = data.frame(matrix(nrow=0, ncol=10))
cname <- c("DeltaTime", "Time To Train", "Time to Predict", 
           "Train Accuracy", "Train AccuracyLower", "Train AccuracyUpper", 
           "Train Sensitivity", "Train Specificity", "Train Pos Pred Value", "Train Neg Pred Value",
           "Train Precision", "Train Recall", "Train F1", "Train Prevalence", "Train Detection Rate",  "Train Detection Prevalence", "Balanced Accuracy",
           "Test Accuracy", "Test AccuracyLower", "Test AccuracyUpper", 
           "Test Sensitivity", "Test Specificity", "Test Pos Pred Value", "Test Neg Pred Value",
           "Test Precision", "Test Recall", "Test F1", "Test Prevalence", "Test Detection Rate",  "Test Detection Prevalence", "Balanced Accuracy"
           )
ctrl <- trainControl(method="repeatedcv", 
                     number = 5,
                     classProbs = T, 
                     summaryFunction = twoClassSummary,
                     verboseIter = TRUE,
                     savePredictions = TRUE,
                     repeats = 10,
                     allowParallel = T)

#smotetrain or rosetrain depending upon ...

for (model in c("glm", "knn", "xgbDART", "nnet", "pcaNNet", "treebag", "nb", "gbm", "AdaBag", "LogitBoost", "C5.0", "svmRadialWeights")){
  cat(model)
  start=Sys.time()
  modeled <- train(CDIFF ~., data = smotetrain, method = model, trControl=ctrl)
  #modeled <- train(CDIFF ~., data = rosetrain, method = model, trControl=ctrl)
  cat("Trained")
  DeltaTimeModel <- Sys.time()-start
  TrainP <- predict(modeled, newdata=smotetrain)
  #TrainP <- predict(modeled, newdata=rosetrain)
  TestP <- predict(modeled, newdata=testing)
  traincm <- confusionMatrix(data=TrainP, smotetrain$CDIFF, positive = "CDIFF")
  #traincm <- confusionMatrix(data=TrainP, rosetrain$CDIFF, positive = "CDIFF")
  testcm <- confusionMatrix(data=TestP, testing$CDIFF, positive = "CDIFF")
  
  DeltaTime <- Sys.time()-start
  DeltaTimeCM <- DeltaTime-DeltaTimeModel
  
  x <- c(format(DeltaTime), format(DeltaTimeModel), format(DeltaTimeCM), 
         traincm$overall[1], traincm$overall[3], traincm$overall[4], traincm$byClass[1], traincm$byClass[2], 
         traincm$byClass[3], traincm$byClass[4], traincm$byClass[5], traincm$byClass[6], traincm$byClass[7],
         traincm$byClass[8], traincm$byClass[9], traincm$byClass[10], traincm$byClass[11],
         testcm$overall[1], testcm$overall[3], testcm$overall[4], testcm$byClass[1], testcm$byClass[2], 
         testcm$byClass[3], testcm$byClass[4], testcm$byClass[5], testcm$byClass[6], testcm$byClass[7],
         testcm$byClass[8], testcm$byClass[9], testcm$byClass[10], testcm$byClass[11])
  x <- t(data.frame(x))
  colnames(x) <- cname
  rownames(x) <- model
  modeloutput <- rbind(modeloutput, x)
  rm(DeltaTime, DeltaTimeCM, DeltaTimeModel, traincm, testcm)
  cat("Data logged")
  write.csv(modeloutput, "Modelvalues_trial_nongenetic_matching_ratio5_IL8simulated_newmethod_included.csv", sep = "|")
  saveRDS(modeled, file=paste0(model, "_nongenetic_matching_ratio5_IL8simulated_newmethod_included.rds"))
}


readRDS("nb.rds")

readRDS("gbm.rds")

modeloutput
# DeltaTime  Time To Train   Time to Predict    Train Accuracy Train Sensitivity Train Specificity
# nb                2.673587 mins  1.697155 mins    0.9764321 mins  0.56022827275394 0.999061308301555 0.194534076464261
# nb1               2.657776 mins  1.686607 mins    0.9711693 mins  0.56022827275394 0.999061308301555 0.194534076464261
# gbm               58.59025 secs  58.26617 secs    0.3240721 secs 0.866289767727139 0.869815195071869 0.863351911606532
# AdaBag            14.50065 mins  14.42407 mins   0.07657557 mins 0.790794421184565  0.80393077148724 0.779847462599003
# LogitBoost        19.00821 secs  14.26081 secs     4.747405 secs  0.84868929838129 0.863185684951599 0.836608976239366
# C5.0              9.609812 mins  9.578692 mins    0.0311197 mins 0.879063441691778  0.88019947198592 0.878116749779994
# svmRadialWeights 1.613348 hours 1.605147 hours 0.008201852 hours 0.854236112963012 0.842827808741566 0.863743033147551
# Test Accuracy  Test Sensitivity   Test Specificity
# nb               0.922132232279922 0.998737591094279 0.0102459016393443
# nb1              0.922132232279922 0.998737591094279 0.0102459016393443
# gbm              0.840082579005876 0.866758478223446  0.522540983606557
# AdaBag           0.783600656397226 0.804154472944282  0.538934426229508
# LogitBoost       0.826531152400614 0.857405175876513  0.459016393442623
# C5.0             0.845376105023556 0.873472198312963  0.510928961748634
# svmRadialWeights 0.815732359324546 0.840477420095255  0.521174863387978




# create ROC curve on testing dataset and calculate 95%CI -----------------


#genetic data without genetic features
modeled_gbm <- readRDS("gbm_genetic_without.rds")
modeled_xgb <- readRDS("xgbDART_genetic_without.rds")
modeled_glm <- readRDS("glm_genetic_without.rds")
modeled_nnet <- readRDS("nnet_genetic_without.rds")
modeled_treebag <- readRDS("treebag_genetic_without.rds")
modeled_c5 <- readRDS("C5.0_genetic_without.rds")
modeled_svm <- readRDS("svmRadialWeights_genetic_without.rds")
modeled_lb <- readRDS("LogitBoost_genetic_without.rds")

ctrl <- trainControl(method="repeatedcv", 
                     number = 5,
                     classProbs = T, 
                     summaryFunction = twoClassSummary,
                     verboseIter = TRUE,
                     savePredictions = TRUE,
                     repeats = 10,
                     allowParallel = T)


probs_gbm <- predict(modeled_gbm, testing, type = "prob")[,2]
pred1_gbm <- prediction(probs_gbm, testing$CDIFF)
perf_gbm <- performance(pred1_gbm, "tpr", "fpr")
perf1_gbm <- performance(pred1_gbm, "auc")

probs_xgb <- predict(modeled_xgb, testing, type = "prob")[,2]
pred1_xgb <- prediction(probs_xgb, testing$CDIFF)
perf_xgb <- performance(pred1_xgb, "tpr", "fpr")
perf1_xgb <- performance(pred1_xgb, "auc")

probs_glm <- predict(modeled_glm, testing, type = "prob")[,2]
pred1_glm <- prediction(probs_glm, testing$CDIFF)
perf_glm <- performance(pred1_glm, "tpr", "fpr")
perf1_glm <- performance(pred1_glm, "auc")

probs_nnet <- predict(modeled_nnet, testing, type = "prob")[,2]
pred1_nnet <- prediction(probs_nnet, testing$CDIFF)
perf_nnet <- performance(pred1_nnet, "tpr", "fpr")
perf1_nnet <- performance(pred1_nnet, "auc")

probs_treebag <- predict(modeled_treebag, testing, type = "prob")[,2]
pred1_treebag <- prediction(probs_treebag, testing$CDIFF)
perf_treebag <- performance(pred1_treebag, "tpr", "fpr")
perf1_treebag <- performance(pred1_treebag, "auc")

probs_c5 <- predict(modeled_c5, testing, type = "prob")[,2]
pred1_c5 <- prediction(probs_c5, testing$CDIFF)
perf_c5 <- performance(pred1_c5, "tpr", "fpr")
perf1_c5 <- performance(pred1_c5, "auc")

probs_svm <- predict(modeled_svm, testing, type = "prob")[,2]
pred1_svm <- prediction(probs_svm, testing$CDIFF)
perf_svm <- performance(pred1_svm, "tpr", "fpr")
perf1_svm <- performance(pred1_svm, "auc")
#roc.test(perf1_gbm, perf1_glm, reuse.auc=FALSE)

probs_lb <- predict(modeled_lb, testing, type = "prob")[,2]
pred1_lb <- prediction(probs_lb, testing$CDIFF)
perf_lb <- performance(pred1_lb, "tpr", "fpr")
perf1_lb <- performance(pred1_lb, "auc")

library(pROC)

gbm_CI <- ci.auc(testing$CDIFF, probs_gbm) 
gbm_CI <- paste(round(gbm_CI[2], 3), "[", round(gbm_CI[1], 3), "-", round(gbm_CI[3], 3), "]")

xgb_CI <- ci.auc(testing$CDIFF, probs_xgb) 
xgb_CI <- paste(round(xgb_CI[2], 3), "[", round(xgb_CI[1], 3), "-", round(xgb_CI[3], 3), "]")

glm_CI <- ci.auc(testing$CDIFF, probs_glm) 
glm_CI <- paste(round(glm_CI[2], 3), "[", round(glm_CI[1], 3), "-", round(glm_CI[3], 3), "]")

nnet_CI <- ci.auc(testing$CDIFF, probs_nnet) 
nnet_CI <- paste(round(nnet_CI[2], 3), "[", round(nnet_CI[1], 3), "-", round(nnet_CI[3], 3), "]")

treebag_CI <- ci.auc(testing$CDIFF, probs_treebag) 
treebag_CI <- paste(round(treebag_CI[2], 3), "[", round(treebag_CI[1], 3), "-", round(treebag_CI[3], 3), "]")

c5_CI <- ci.auc(testing$CDIFF, probs_c5) 
c5_CI <- paste(round(c5_CI[2], 3), "[", round(c5_CI[1], 3), "-", round(c5_CI[3], 3), "]")

svm_CI <- ci.auc(testing$CDIFF, probs_svm) 
svm_CI <- paste(round(svm_CI[2], 3), "[", round(svm_CI[1], 3), "-", round(svm_CI[3], 3), "]")

lb_CI <- ci.auc(testing$CDIFF, probs_lb) 
lb_CI <- paste(round(lb_CI[2], 3), "[", round(lb_CI[1], 3), "-", round(lb_CI[3], 3), "]")


getwd()
tiff("Rplot.tiff", units="in", width=8, height=8, res=300)  ##type = "s" in the following code, can be changed 

plot(perf_gbm@x.values[[1]], perf_gbm@y.values[[1]], type = 'l', ylab = perf_gbm@y.name, xlab = perf_gbm@x.name, col = '1', lwd = 2)
lines(perf_xgb@x.values[[1]], perf_xgb@y.values[[1]], col = '2', lty = 2)
lines(perf_glm@x.values[[1]], perf_glm@y.values[[1]], col = '3', lty = 2)
lines(perf_nnet@x.values[[1]], perf_nnet@y.values[[1]], col = '4', lty = 2)
lines(perf_treebag@x.values[[1]], perf_treebag@y.values[[1]], col = '5', lty = 2)
lines(perf_c5@x.values[[1]], perf_c5@y.values[[1]], col = '6', lty = 2)
lines(perf_svm@x.values[[1]], perf_svm@y.values[[1]], col = '7', lty = 2)
lines(perf_lb@x.values[[1]], perf_lb@y.values[[1]], col = '8', lty = 2)

abline(a = 0, b = 1)

model1.lab <- paste("gbm AUC", gbm_CI)
model2.lab <- paste('xgbDART AUC', xgb_CI)
model3.lab <- paste('glm AUC', glm_CI)
model4.lab <- paste('nnet AUC', nnet_CI)
model5.lab <- paste('treebag AUC', treebag_CI)
model6.lab <- paste('C5.0 AUC', c5_CI)
model7.lab <- paste('svmRadialWeights AUC', svm_CI)
model8.lab <- paste('LogitBoost AUC', lb_CI)

legend('bottomright', c(model1.lab, model2.lab, model3.lab, model4.lab, model5.lab, model6.lab, model7.lab, model8.lab), col = c(1:8), lwd = c(2, 1, 1, 1, 1, 1,1,1), lty = c(1, 2, 2, 2, 2, 2, 2,2), cex = 1, bty = 'n')
dev.off()

# Delong's test for AUC -----------------------------------------------------

probs_gbm <- predict(modeled_gbm, testing, type = "prob")[,2]
pred1_gbm <- prediction(probs_gbm, testing$CDIFF)
perf_gbm <- performance(pred1_gbm, "tpr", "fpr")
perf1_gbm <- performance(pred1_gbm, "auc")

result.predicted.prob.gbm <- predict(modeled_gbm, testing, type = "prob") # Prediction

result.roc.gbm <- roc(testing$CDIFF, result.predicted.prob.gbm$CDIFF, smooth = FALSE) # Draw ROC curve.
plot(result.roc.gbm, print.thres="best", print.thres.best.method="closest.topleft")

result.predicted.prob.glm <- predict(modeled_glm, testing, type = "prob") # Prediction

result.roc.glm <- roc(testing$CDIFF, result.predicted.prob.glm$CDIFF) # Draw ROC curve.
plot(result.roc.glm, print.thres="best", print.thres.best.method="closest.topleft")


result.coords <- coords(result.roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
print(result.coords)#to get threshold and accuracy


roc.test(result.roc.gbm, result.roc.glm, method="bootstrap")

roc.test(result.roc.gbm, result.roc.glm, paired=TRUE, plot = TRUE, smooth = FALSE, method="delong")

roc.test(result.roc.gbm, result.roc.glm, paired=TRUE, method="venkatraman") #very slow
roc.test(result.roc.gbm, result.roc.glm, paired=FALSE, method="specificity", specificity=0.9)


# Comparison on specificity and sensitivity
roc.test(result.roc.gbm, result.roc.glm, method="specificity", specificity=0.9)
roc.test(result.roc.gbm, result.roc.glm, method="sensitivity", sensitivity=0.9)

#multiple ROC curve comparison
result.predicted.prob.gbm <- predict(modeled_gbm, testing, type = "prob") # Prediction
result.roc.gbm <- roc(testing$CDIFF, result.predicted.prob.gbm$CDIFF, smooth = FALSE)
result.predicted.prob.xgb <- predict(modeled_xgb, testing, type = "prob") # Prediction
result.roc.xgb <- roc(testing$CDIFF, result.predicted.prob.xgb$CDIFF, smooth = FALSE)
result.predicted.prob.glm <- predict(modeled_glm, testing, type = "prob") # Prediction
result.roc.glm <- roc(testing$CDIFF, result.predicted.prob.glm$CDIFF, smooth = FALSE)
result.predicted.prob.nnet <- predict(modeled_nnet, testing, type = "prob") # Prediction
result.roc.nnet <- roc(testing$CDIFF, result.predicted.prob.nnet$CDIFF, smooth = FALSE)
result.predicted.prob.treebag <- predict(modeled_treebag, testing, type = "prob") # Prediction
result.roc.treebag <- roc(testing$CDIFF, result.predicted.prob.treebag$CDIFF, smooth = FALSE)
result.predicted.prob.c5 <- predict(modeled_c5, testing, type = "prob") # Prediction
result.roc.c5 <- roc(testing$CDIFF, result.predicted.prob.c5$CDIFF, smooth = FALSE)
result.predicted.prob.svm <- predict(modeled_svm, testing, type = "prob") # Prediction
result.roc.svm <- roc(testing$CDIFF, result.predicted.prob.svm$CDIFF, smooth = FALSE)
result.predicted.prob.lb <- predict(modeled_lb, testing, type = "prob") # Prediction
result.roc.lb <- roc(testing$CDIFF, result.predicted.prob.lb$CDIFF, smooth = FALSE)

all_tests <- combn(
  list(
    "gbm" = result.roc.gbm,
    "xgb" = result.roc.xgb,
    "glm" = result.roc.glm,
    "nnet" = result.roc.nnet,
    "treebag" = result.roc.treebag,
    "c5" = result.roc.c5,
    "svm" = result.roc.svm,
    "lb" = result.roc.lb
  ),
  FUN = function(x, ...) roc.test(x[[1]], x[[2]]),
  m = 2,
  simplify = FALSE, 
  reuse.auc = TRUE, 
  method = "delong", #or changed into "bootstrap", boot.n = 10000
  paired=TRUE,
  na.rm = TRUE
)

all_tests[[1]]

tests_names <- combn(
  list("gbm", "xgb", "glm", "nnet", "treebag", "c5", "svm", "lb"), 
  m = 2, 
  FUN = paste, 
  simplify = TRUE, 
  collapse = "_"
)
all_tests <- setNames(all_tests, tests_names)

all_tests$gbm_xgb$p.value
all_tests$gbm_xgb$estimate
all_tests$gbm_xgb$statistic
gbm_xgb <- c(all_tests$gbm_xgb$p.value, all_tests$gbm_xgb$estimate, all_tests$gbm_xgb$statistic)
all_tests

gbm_xgb <- c(all_tests$gbm_xgb$p.value, all_tests$gbm_xgb$estimate, all_tests$gbm_xgb$statistic)
gbm_glm <- c(all_tests$gbm_glm$p.value, all_tests$gbm_glm$estimate, all_tests$gbm_glm$statistic)
gbm_nnet <- c(all_tests$gbm_nnet$p.value, all_tests$gbm_nnet$estimate, all_tests$gbm_nnet$statistic)
gbm_treebag <- c(all_tests$gbm_treebag$p.value, all_tests$gbm_treebag$estimate, all_tests$gbm_treebag$statistic)
gbm_c5 <- c(all_tests$gbm_c5$p.value, all_tests$gbm_c5$estimate, all_tests$gbm_c5$statistic)
gbm_svm <- c(all_tests$gbm_svm$p.value, all_tests$gbm_svm$estimate, all_tests$gbm_svm$statistic)
gbm_lb <- c(all_tests$gbm_lb$p.value, all_tests$gbm_lb$estimate, all_tests$gbm_lb$statistic)

vectorList <- list(gbm_xgb, gbm_glm, gbm_nnet, gbm_treebag, gbm_c5, gbm_svm, gbm_lb)
result_final <- as.data.frame(do.call(rbind, vectorList))

colnames(result_final) <- c("p_value", "AUC_gbm", "AUC_others", "Z_Statistics")
result_final$algorithms <- c("xgb", "glm", "nnet", "treebag", "c5", "svm", "lb")
result_final

#        p_value   AUC_gbm AUC_others Z_Statistics algorithms
# 1 0.9582151730 0.7098534  0.7095361   0.05239348        xgb
# 2 0.0003090291 0.7098534  0.6828235   3.60761121        glm
# 3 0.2133498232 0.7098534  0.7008965   1.24440696       nnet
# 4 0.0011453021 0.7098534  0.6734246   3.25216327    treebag
# 5 0.1668406473 0.7098534  0.7012011   1.38242695         c5
# 6 0.0186298102 0.7098534  0.6899299   2.35285674        svm
# 7 0.0011733764 0.7098534  0.6785956   3.24527405         lb



# comparing model with or without genetic feature -------------------------


#for genetic sample without genetic features, match age and sex with ratio 5 and sampling up for training,  xgbDART_nongenetic_matching_ratio5.rds
modeled_gbm <- readRDS("gbm_nongenetic_matching_ratio10.rds")
modeled_xgb <- readRDS("xgbDART_nongenetic_matching_ratio10.rds")
modeled_glm <- readRDS("glm_nongenetic_matching_ratio10.rds")
modeled_nnet <- readRDS("nnet_nongenetic_matching_ratio10.rds")
modeled_treebag <- readRDS("treebag_nongenetic_matching_ratio10.rds")
modeled_c5 <- readRDS("C5.0_nongenetic_matching_ratio10.rds")
modeled_svm <- readRDS("svmRadialWeights_nongenetic_matching_ratio10.rds")
modeled_lb <- readRDS("LogitBoost_nongenetic_matching_ratio10.rds")

#for genetic sample with genetic features, match age and sex with ratio 5 and sampling up for training
modeled_gbm_X4 <- readRDS("gbm_genetic_matching_ratio10.rds")
modeled_xgb_X4 <- readRDS("xgbDART_genetic_matching_ratio10.rds")
modeled_glm_X4 <- readRDS("glm_genetic_matching_ratio10.rds")
modeled_nnet_X4 <- readRDS("nnet_genetic_matching_ratio10.rds")
modeled_treebag_X4 <- readRDS("treebag_genetic_matching_ratio10.rds")
modeled_c5_X4 <- readRDS("C5.0_genetic_matching_ratio10.rds")
modeled_svm_X4 <- readRDS("svmRadialWeights_genetic_matching_ratio10.rds")
modeled_lb_X4 <- readRDS("LogitBoost_genetic_matching_ratio10.rds")

#multiple ROC curve comparison
result.predicted.prob.gbm_X4 <- predict(modeled_gbm_X4, testing, type = "prob") # Prediction
result.roc.gbm_X4 <- roc(testing$CDIFF, result.predicted.prob.gbm_X4$CDIFF, smooth = FALSE)
result.predicted.prob.xgb_X4 <- predict(modeled_xgb_X4, testing, type = "prob") # Prediction
result.roc.xgb_X4 <- roc(testing$CDIFF, result.predicted.prob.xgb_X4$CDIFF, smooth = FALSE)
result.predicted.prob.glm_X4 <- predict(modeled_glm_X4, testing, type = "prob") # Prediction
result.roc.glm_X4 <- roc(testing$CDIFF, result.predicted.prob.glm_X4$CDIFF, smooth = FALSE)
result.predicted.prob.nnet_X4 <- predict(modeled_nnet_X4, testing, type = "prob") # Prediction
result.roc.nnet_X4 <- roc(testing$CDIFF, result.predicted.prob.nnet_X4$CDIFF, smooth = FALSE)
result.predicted.prob.treebag_X4 <- predict(modeled_treebag_X4, testing, type = "prob") # Prediction
result.roc.treebag_X4 <- roc(testing$CDIFF, result.predicted.prob.treebag_X4$CDIFF, smooth = FALSE)
result.predicted.prob.c5_X4 <- predict(modeled_c5_X4, testing, type = "prob") # Prediction
result.roc.c5_X4 <- roc(testing$CDIFF, result.predicted.prob.c5_X4$CDIFF, smooth = FALSE)
result.predicted.prob.svm_X4 <- predict(modeled_svm_X4, testing, type = "prob") # Prediction
result.roc.svm_X4 <- roc(testing$CDIFF, result.predicted.prob.svm_X4$CDIFF, smooth = FALSE)
result.predicted.prob.lb_X4 <- predict(modeled_lb_X4, testing, type = "prob") # Prediction
result.roc.lb_X4 <- roc(testing$CDIFF, result.predicted.prob.lb_X4$CDIFF, smooth = FALSE)

result.predicted.prob.gbm <- predict(modeled_gbm, testing, type = "prob") # Prediction
result.roc.gbm <- roc(testing$CDIFF, result.predicted.prob.gbm$CDIFF, smooth = FALSE)
result.predicted.prob.xgb <- predict(modeled_xgb, testing, type = "prob") # Prediction
result.roc.xgb <- roc(testing$CDIFF, result.predicted.prob.xgb$CDIFF, smooth = FALSE)
result.predicted.prob.glm <- predict(modeled_glm, testing, type = "prob") # Prediction
result.roc.glm <- roc(testing$CDIFF, result.predicted.prob.glm$CDIFF, smooth = FALSE)
result.predicted.prob.nnet <- predict(modeled_nnet, testing, type = "prob") # Prediction
result.roc.nnet <- roc(testing$CDIFF, result.predicted.prob.nnet$CDIFF, smooth = FALSE)
result.predicted.prob.treebag <- predict(modeled_treebag, testing, type = "prob") # Prediction
result.roc.treebag <- roc(testing$CDIFF, result.predicted.prob.treebag$CDIFF, smooth = FALSE)
result.predicted.prob.c5 <- predict(modeled_c5, testing, type = "prob") # Prediction
result.roc.c5 <- roc(testing$CDIFF, result.predicted.prob.c5$CDIFF, smooth = FALSE)
result.predicted.prob.svm <- predict(modeled_svm, testing, type = "prob") # Prediction
result.roc.svm <- roc(testing$CDIFF, result.predicted.prob.svm$CDIFF, smooth = FALSE)
result.predicted.prob.lb <- predict(modeled_lb, testing, type = "prob") # Prediction
result.roc.lb <- roc(testing$CDIFF, result.predicted.prob.lb$CDIFF, smooth = FALSE)

all_tests <- combn(
  list(
    "gbm" = result.roc.gbm,
    "xgb" = result.roc.xgb,
    "glm" = result.roc.glm,
    "nnet" = result.roc.nnet,
    "treebag" = result.roc.treebag,
    "c5" = result.roc.c5,
    "svm" = result.roc.svm,
    "lb" = result.roc.lb,
    "gbm_X4" = result.roc.gbm_X4,
    "xgb_X4" = result.roc.xgb_X4,
    "glm_X4" = result.roc.glm_X4,
    "nnet_X4" = result.roc.nnet_X4,
    "treebag_X4" = result.roc.treebag_X4,
    "c5_X4" = result.roc.c5_X4,
    "svm_X4" = result.roc.svm_X4,
    "lb_X4" = result.roc.lb_X4
  ),
  FUN = function(x, ...) roc.test(x[[1]], x[[2]]),
  m = 2,
  simplify = FALSE, 
  reuse.auc = TRUE, 
  method = "delong", 
  paired=TRUE,
  na.rm = TRUE
)

all_tests[[1]]

tests_names <- combn(
  list("gbm", "xgb", "glm", "nnet", "treebag", "c5", "svm", "lb", "gbm_X4X6", "xgb_X4X6", "glm_X4X6", "nnet_X4X6", "treebag_X4X6", "c5_X4X6", "svm_X4X6", "lb_X4X6"), 
  m = 2, 
  FUN = paste, 
  simplify = TRUE, 
  collapse = "_"
)
all_tests <- setNames(all_tests, tests_names)


all_tests$gbm_gbm_X4
gbm_gbm_X4 <- c(all_tests$gbm_gbm_X4$p.value, all_tests$gbm_gbm_X4$estimate, all_tests$gbm_gbm_X4$statistic)
xgb_xgb_X4 <- c(all_tests$xgb_xgb_X4$p.value, all_tests$xgb_xgb_X4$estimate, all_tests$xgb_xgb_X4$statistic)
glm_glm_X4 <- c(all_tests$glm_glm_X4$p.value, all_tests$glm_glm_X4$estimate, all_tests$glm_glm_X4$statistic)
nnet_nnet_X4 <- c(all_tests$nnet_nnet_X4$p.value, all_tests$nnet_nnet_X4$estimate, all_tests$nnet_nnet_X4$statistic)
treebag_treebag_X4 <- c(all_tests$treebag_treebag_X4$p.value, all_tests$treebag_treebag_X4$estimate, all_tests$treebag_treebag_X4$statistic)
c5_c5_X4 <- c(all_tests$c5_c5_X4$p.value, all_tests$c5_c5_X4$estimate, all_tests$c5_c5_X4$statistic)
svm_svm_X4 <- c(all_tests$svm_svm_X4$p.value, all_tests$svm_svm_X4$estimate, all_tests$svm_svm_X4$statistic)
lb_lb_X4 <- c(all_tests$lb_lb_X4$p.value, all_tests$lb_lb_X4$estimate, all_tests$lb_lb_X4$statistic)

vectorList <- list(gbm_gbm_X4, xgb_xgb_X4, glm_glm_X4, nnet_nnet_X4, treebag_treebag_X4, c5_c5_X4, svm_svm_X4, lb_lb_X4)
result_final <- as.data.frame(do.call(rbind, vectorList))

colnames(result_final) <- c("p_value", "AUC_nongenetics", "AUC_genetics", "Z_Statistics")
result_final$algorithms <- c("gbm", "xgb", "glm", "nnet", "treebag", "c5", "svm", "lb")
result_final


# create feature importance bar and raddar plots ---------------------------------------


## barplot
temp <- readRDS("gbm_genetic_without_new.rds")
#estimate variable importance
importance <- varImp(temp, scale = T)
print(importance)

tiff("Rplot_importanceplot_training_gbm_genetic_without_new.tiff", units="in", width=9, height=8, res=300) 
plot(importance, xlab = "importance(gbm)", main = "Feature importance for the cohort w/ genetic data but not included")
dev.off()



#raddar plot

importance <- varImp(temp, scale = T) #make sure scale is true
print(importance)

xgbDART <- importance$importance

xgbDART$rs2227306 <- "without" #change from with to without
xgbDART$source <- "genetics" #change from nongenetics to genetics
xgbDART$matching <- "ratio_5" #no, ratio_10, ratio_5
xgbDART$algorithms <- "xgbDART"
xgbDART$feature <- row.names(xgbDART)

gbm <- importance$importance

gbm$rs2227306 <- "without"  #change from with to without
gbm$source <- "genetics"
gbm$matching <- "ratio_5"
gbm$algorithms <- "gbm"
gbm$feature <- row.names(gbm)

glm <- importance$importance

glm$rs2227306 <- "without" #change from with to without
glm$source <- "genetics"
glm$matching <- "ratio_5"
glm$algorithms <- "glm"
glm$feature <- row.names(glm)

nnet <- importance$importance

nnet$rs2227306 <- "without" #change from with to without
nnet$source <- "genetics"
nnet$matching <- "ratio_5"
nnet$algorithms <- "nnet"
nnet$feature <- row.names(nnet)

df <- do.call("rbind", list(gbm, xgbDART, glm, nnet))

#library(data.table)
setDT(df)
#df_wide <- dcast(df, algorithms ~ feature, value.var = "Overall")
df_wide <- dcast(df,  ... ~ feature, value.var = "Overall") #keep other columns
colnames(df_wide)
#with genotype data
colnames(df_wide) <- c("rs2227306", "source", "matching", "algorithms", "Antibiotics", "Age at index", "Chemotherapy", "HIV", "Inflamatory Bowel Dis", "Proton Pump Inhibitors", "Sex", "Type 2 Diabetic", "Anti-TNF Med", "rs2227306/T", "Corticosteroid", "Transplant Med")

df_wide <- df_wide[, c("rs2227306", "source", "matching", "algorithms", "Antibiotics", "Age at index", "Chemotherapy", "HIV", "Inflamatory Bowel Dis", "Proton Pump Inhibitors", "Sex", "Type 2 Diabetic", "Anti-TNF Med", "Corticosteroid", "Transplant Med")]

df_wide_select <- df_wide[, 4:ncol(df_wide)]

library(ggradar)
tiff("Rplot_radar_importanceplot_training_nongenetic_nomatching_IL8simulated_newmethod_notincluded.tiff", units="in", width=12, height=10, res=300) 
ggradar(
  df_wide_select, 
  values.radar = c("0", "50", "100"),
  grid.min = 0, grid.mid = 50, grid.max = 100,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = c("#00AFBB", "#E7B800", "#FC4E07", "#EE82EE"),
  #group.colours = c("#E7B800", "#FC4E07", "#EE82EE"),
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom"
)
dev.off()

df_wide_select <- df_wide[df_wide$algorithms != "nnet", 4:ncol(df_wide)]
tiff("Rplot_radar3_importanceplot_training_nongenetic_nomatching_IL8simulated_newmethod_notincluded.tiff", units="in", width=12, height=10, res=300) 
ggradar(
  df_wide_select, 
  values.radar = c("0", "50", "100"),
  grid.min = 0, grid.mid = 50, grid.max = 100,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  #group.colours = c("#00AFBB", "#E7B800", "#FC4E07", "#EE82EE"),
  group.colours = c("#E7B800", "#FC4E07", "#EE82EE"),
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom"
)
dev.off()


table <- read.table("Table_AUC_comparison_newmethod.txt", header = T, stringsAsFactors = F)

# AUC/accuracy forest plot ---------------------------------------------------------


head(table)
# Index  p_value       AUC AUC_lower AUC_upper Z_Statistics algorithms     source matching rs2227306
# 1     1 6.24e-03 0.8177056     0.807     0.829   -2.7350757        gbm nongenetic       no   without
# 2     2 9.10e-01 0.8193095     0.808     0.831    0.1130072        xgb nongenetic       no   without
# 3     3 1.22e-01 0.7505024     0.737     0.764   -1.5478655        glm nongenetic       no   without
# 4     4 4.14e-07 0.8115372     0.800     0.823    5.0623630       nnet nongenetic       no   without
# 5     5 2.24e-02 0.7742220     0.761     0.787    2.2833437    treebag nongenetic       no   without
# 6     6 6.03e-01 0.8137551     0.802     0.825   -0.5201652         c5 nongenetic       no   without

table_select <- table[table$algorithms %in% c("glm", "gbm", "xgb"), ]

tiff("test1.tiff", units="in", width=12, height=8, res=600)
ggplot(table, aes(x = matching, y = AUC, col = rs2227306, shape = source)) + 
  geom_errorbar(aes(ymin = AUC_lower, ymax = AUC_upper), height = 0.1, size = 0.5, stat="identity", position=position_dodge(0.6)) + 
  geom_point(aes(size = -log10(p_value)), stat="identity", position = position_dodge(0.6))  +
  scale_shape_manual(values=c("\u25CF","\u25CB")) +  #circle
  scale_y_continuous(limits = c(0.55, 0.85)) +
  geom_hline(aes(yintercept = 0.75), linetype='dashed', col = 'black', width = 1) +
  facet_wrap(. ~ algorithms, scales="free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 12, colour = "black")) +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(face="bold", color="#993333", 
                                 size=12, angle=45),
        legend.text=element_text(size = 12)) + 
  ylab("ROCAUC & 95% CI") +
  xlab("Cohort w/ or w/o propensity score matching")
dev.off()

tiff("test.tiff", units="in", width=12, height=8, res=600)
ggplot(table_select, aes(x = matching, y = AUC, col = rs2227306, shape = source)) + 
  geom_errorbar(aes(ymin = AUC_lower, ymax = AUC_upper), height = 0.1, size = 0.5, stat="identity", position=position_dodge(0.6)) + 
  geom_point(aes(size = -log10(p_value)), stat="identity", position = position_dodge(0.6))  +
  scale_shape_manual(values=c("\u25CF","\u25CB")) +  #circle
  scale_y_continuous(limits = c(0.6, 0.85)) +
  geom_hline(aes(yintercept = 0.75), linetype='dashed', col = 'black', width = 1) +
  facet_wrap(. ~ algorithms, scales="free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 12, colour = "black")) +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(face="bold", color="#993333", 
                                 size=12, angle=45),
        legend.text=element_text(size = 12)) + 
  ylab("ROCAUC & 95% CI") +
  xlab("Cohort w/ or w/o propensity score matching")
dev.off()

#make the figure for Accuracy comparison
table2 <- read.csv("table2_update.csv", header = T, stringsAsFactors = F)

head(table2)
unique(table2$algorithms)
table2_select <- table2[table2$algorithms %in% c("glm", "gbm", "xgbDART"), ]
table2_select <- table2[table2$algorithms %in% c("glm", "gbm", "xgbDART", "nnet", "treebag", "C5.0", "LogitBoost", "svmRadialWeights"), ]
head(table2_select, 14)
colnames(table2_select)
tiff("test4.tiff", units="in", width=12, height=8, res=600)
ggplot(table2_select, aes(x = matching, y = Test_Accuracy, col = rs2227306, shape = source)) + 
  geom_errorbar(aes(ymin = Test_AccuracyLower, ymax = Test_AccuracyUpper), height = 0.1, size = 0.5, stat="identity", position=position_dodge(0.6)) + 
  geom_point(stat="identity", position = position_dodge(0.6))  +
  scale_shape_manual(values=c("\u25BC","\u25B2")) + #triangle
  scale_y_continuous(limits = c(0.6, 0.85)) +
  geom_hline(aes(yintercept = 0.75), linetype='dashed', col = 'black', width = 1) +
  facet_wrap(. ~ algorithms, scales="free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 12, colour = "black")) +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(face="bold", color="#993333", 
                                 size=12, angle=45),
        legend.text=element_text(size = 12)) + 
  ylab("Accuracy & 95% CI") +
  xlab("Cohort w/ or w/o propensity score matching")
dev.off()



# make figure for f1 score comparison across different oversampling --------



#f1 score for training datasets
tiff("test10.tiff", units="in", width=12, height=8, res=600)
ggplot(table2_select, aes(x = matching, y = Train_F1, col = rs2227306, fill = source)) + 
  geom_bar(stat="identity", width = 0.5, position = position_dodge(0.6))  +
  scale_y_continuous(limits = c(0.0, 1.0)) +
  scale_color_grey() +
  facet_wrap(. ~ algorithms, scales="free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 12, colour = "black")) +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(face="bold", color="#993333", 
                                 size=12, angle=45),
        legend.text=element_text(size = 12)) + 
  ylab("F1 score in the training dataset") +
  xlab("Cohort w/ or w/o propensity score matching")
dev.off()

head(table2_select)
#Recall for testing datasets
tiff("test10.tiff", units="in", width=12, height=8, res=600)
ggplot(table2_select, aes(x = matching, y = Test_Recall, col = rs2227306, fill = source)) + 
  geom_bar(stat="identity", width = 0.5, position = position_dodge(0.6))  +
  scale_y_continuous(limits = c(0.0, 0.8)) +
  scale_color_grey() +
  facet_wrap(. ~ algorithms, scales="free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 12, colour = "black")) +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(face="bold", color="#993333", 
                                 size=12, angle=45),
        legend.text=element_text(size = 12)) + 
  ylab("Recall in the testing dataset") +
  xlab("Cohort w/ or w/o propensity score matching")
dev.off()

#Recall for training datasets
tiff("test11.tiff", units="in", width=12, height=8, res=600)
ggplot(table2_select, aes(x = matching, y = Train_Recall, col = rs2227306, fill = source)) + 
  geom_bar(stat="identity", width = 0.5, position = position_dodge(0.6))  +
  scale_y_continuous(limits = c(0.0, 1.0)) +
  scale_color_grey() +
  facet_wrap(. ~ algorithms, scales="free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 12, colour = "black")) +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(face="bold", color="#993333", 
                                 size=12, angle=45),
        legend.text=element_text(size = 12)) + 
  ylab("Recall in the training dataset") +
  xlab("Cohort w/ or w/o propensity score matching")
dev.off()

#make the figure for F1 comparison between rose and smote
table3 <- read.csv("D:/Geisinger/X19981_backup_10072020/Desktop/C_DIFF_machine_learning/table3.csv", header = T, stringsAsFactors = F)
head(table3)
unique(table3$algorithms)
table3_select <- table3[table3$algorithms %in% c("glm", "gbm", "xgbDART"), ]
table3_select <- table3[table3$algorithms %in% c("glm", "gbm", "xgbDART", "nnet", "treebag", "C5.0", "LogitBoost", "svmRadialWeights"), ]
head(table3_select, 14)
colnames(table3_select)

head(table3_select)

#f1 score for testing datasets
tiff("test9.tiff", units="in", width=12, height=8, res=600)
ggplot(table3_select, aes(x = Sampling_strategy, y = Test.F1, fill = rs2227306)) + 
  geom_bar(stat="identity", width = 0.5, position = position_dodge(0.6))  +
  scale_y_continuous(limits = c(0.0, 0.80)) +
  scale_color_grey() +
  facet_wrap(. ~ algorithms, scales="free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 12, colour = "black")) +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(face="bold", color="#993333", 
                                 size=12, angle=45),
        legend.text=element_text(size = 12)) + 
  ylab("F1 score in the testing dataset") +
  xlab("Cohort w/ the different sampling strategyy")
dev.off()

#f1 score for training datasets
tiff("test10.tiff", units="in", width=12, height=8, res=600)
ggplot(table3_select, aes(x = Sampling_strategy, y = Train.F1, fill = rs2227306)) + 
  geom_bar(stat="identity", width = 0.5, position = position_dodge(0.6))  +
  scale_y_continuous(limits = c(0.0, 1.0)) +
  scale_color_grey() +
  facet_wrap(. ~ algorithms, scales="free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 12, colour = "black")) +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(face="bold", color="#993333", 
                                 size=12, angle=45),
        legend.text=element_text(size = 12)) + 
  ylab("F1 score in the training dataset") +
  xlab("Cohort w/ the different sampling strategy")
dev.off()


F1_score<-function(Recall, Precision)   {
  F1<-2*Recall*Precision/(Recall+Precision)
  return(F1)
}

recall_test <- sensitivity(preds1, test_data$target)
precision_test <- posPredValue(preds1, test_data$target)

F1_model_c5<-F1_score(recall_test,precision_test)

print(F1_model_c5)



