Here are the explanation for every part of the framework we developed in this study:

1) Sample Generation: To generate RMPEA samples, you need to execute SampleGeneration.ipynb. Part A generates Samples.txt that are the 50,000 samples at T=850K and part B 
generates Samples_T.txt that includes those 50,000 samples for different tempeatures.

2) Phase Labeling: To label the samples generated in the earlier step with expected phases using Thermo-CalC, you need to execute PhaseLabeling.ipynb. This code generates Phases.txt. 
Note that we are not sharing Phases.txt and Phases_T.txt to comply with Thermo-Calc’s guidelines for sharing of calculated data.

3) Feature Labeling: To label the samples with the 51 input features, you need to execute FeatureLabeling.ipynb. Part A generates the primary features labels.txt and then part B 
generates the remaining ones and results labels_extra.txt that consists of both primary and secondary features. Note that to properly execute part A, you need to first create Phi.txt
using publicly available ASAP code.

4) Feature Engineering: In our study, we conduct two steps of feature engineering and you need to execute FeatureEngineering.ipynb to implement that. As for the first feature 
engineering, you need to execute Part A for one time following with executing Part B for multiple times untill stop reeiving the message "Number of deleted RMPEAs=". Then, you can 
execute Part C to generate the final dataset. This dataset further goes through the second feature engineering to remove those highly correlated features. To do the second feature 
engienering, you thus need to execute part A of it to see which features need to be removed and then execute Part B to complete the feature engineering process and generate the final 
dataset needed to start the training process as follows: xlo_1.npy,ylo.npy, xlo_test1.npy, ylo_test.npy. One should note that for the ease of following the name of the files, in the
remaining steps, while we use xlo_1.npy and xlo_test1.npy files, we re-name them as xlo.npy and xlo_test.npy

5) Training Process: To first do the Bayesian optimization, you need to execute MLP_Finetuning.ipynb. Following choosing the appropriate artitecture, MLP.ipynb is then needed to execute.

6) Performance Evaluation: Figures 2-4 and Table 1 can be obtained using PerformanceEvaluation.ipynb.

7) Spider Plot: Figure 5 of the manscuript is prepared via executing SpiderPlot.ipynb

8) Empirical Equation: Figure 6 of the mansucript is prepared via executing EmpiricalEquation.ipynb. Note that you need to first train the SVM using SVM.ipynb