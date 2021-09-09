# Prediction-of-CDI

## Prediction of Clostridioides difficile infection using a genetic variant from IL8 and clinical risk factors from EHR
#### Jiang Li1, Vaibhav Sharma2, Durgesh Chaudhary3, Vishakha Sharma4, Venkatesh Avula1, Harshit S. Khara5  , Donna M. Wolk6, Ramin Zand4, Vida Abedi1,7
##### Geisinger Healthcare System

## Introduction: Clostridioides difficile(C. difficile), is a major cause of hospital-associated and community-acquired diarrhea. Our goal was to develop a prediction model of Clostridioides difficile infection (CDI) by integration of common clinical risk factors extracted from electronic health record (EHR) and genetic risk factor, in general as well as age and sex matched subpopulations. 

## Methods: We used patient-level data, in a retrospective design, from EHR (January 1, 2009 to December 31, 2017), eight machine learning algorithms (Logistic Regression, Gradient Boosted Classification, Extreme Gradient Boosting, Bagging for tree, Neural Network, C5.0, LogitBoost, and Support Vector Machine), inclusion or exclusion of genetic data – rs2227306(IL8), two simulation strategies for samples without genetic data, and two oversampling strategies to handle data imbalance, to model symptomatic CDI. The phenotyping algorithm to identify CDI was adapted from eMERGE. The data was split into training(70%) and testing(30%). Propensity score matching(PSM) was used to control confounding effects of age and sex. We evaluated the model performance using the area under the receiver operating characteristic (AUROC) and compared AUROCs from two models using the DeLong test.   

## Result: We included data from 5911 cases and 69086 controls. Inclusion of rs2227306(IL8) in the optimal models outperformed the corresponding base model (AUROCgbm=0.72[0.694-0.746] vs. 0.710[0.683-0.737], p = 0.006). The contribution of the genetic feature to the prediction became uncertain and algorithm dependent, particularly in elderly patients after PSM.
## Conclusion: The cost and benefit of including genetic feature into the prediction models should be thoroughly evaluated to determine its value in general as well as subpopulations. 
