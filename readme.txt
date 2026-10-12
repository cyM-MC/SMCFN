1 Titile: Tailored Low-Rank Matrix Factorization for Similarity Matrix Completion

2 Enviroment: All the algorithms are implemented in MATLAB R2022a. 

3 Detailed Description:
3.1 MDI folder contains all the missing data imputation (MDI) baseline methods, including Mean,LR,GR(original GROUSE),KFMC. 
3.2 MC folder contains all the matrix completion (MC) baseline methods, including BPMF, FGSR, NNFN, facNNFN.
3.3 SMC folder contains all the similarity matrix completion methods, inclduing CMC, DMC, FMC(original SMC). 
3.3 Ablation folder contains all ablation study, including SMC_F, SMC_NR, SMC_GD, SMC_ON, SMC_fON.
3.4 Ours folder contains our proposed SMCFN and SMCRN.


4 How to evaluate? 
4.1 Firstly, all hyperparameters setting, including weighted parameter $\lambda$, stepsize $\gamma$, rank $r$, and iterations number $T$, can be found in testMDI.m/testMC.m file. 


4.2 Test the effectivenss/efficiency for all comparison methods. 
testMDI.m: test code for comparison missing data imputation (MDI) methods.
testMC.m: test code for comparison matrix completion (MC) methods and similarity matrix completion (SMC) methods.
testablation.m: test code for comparison Ablation methods.

Given the number of search candidates 'noff' and query items 'non', missing ratio 'p_list', with the various settings of hyperparameter, we can test the average performance and standard devation of RMSE, Recall, nDCG of each baseline methods. 


For example, we first load the toydata as the entire database. Then, we randomly select 1,000 data samples as search candidates, 200 data samples as query items from the entire database. Missing ratio $\rho$= 0.7, rank $r$= 100, stepsize $\gamma$=0.001, weighted parameter $\lambda$=0.001. Then, by running testMC.m, the RMSE,Recall,nDCG,Time are shown as: 



                   CMC        DMC        FMC       NNFN      facNNFN     FGSR       BPMF       SMCFN      SMCRN 
                 _______    _______    _______    _______    _______    _______    _______    _______    _______

    RMSE          0.7184          1    0.73495     1.0574     3.8021     0.7624     0.7224    0.69668    0.69227
    Recall       0.74194    0.73944    0.65944     0.5825    0.51028    0.61326    0.61326    0.71493    0.68875
    nDCG         0.62146    0.61857    0.62934    0.62521     0.6769    0.73401    0.73401    0.82379    0.82053
    Time             441     415.31     421.11     115.78     118.75    110.625    114.925     111.12     118.12
    stdRMSE      0.02534          0    0.01194    0.03271    0.05932    0.03374    0.03374       0.01    0.00992
    stdRecall    0.00983    0.01276    0.00726    0.00724    0.01238    0.01124    0.01124    0.01035    0.01074
    stdnDCG       0.0038    0.00493    0.00389    0.00414    0.00531    0.01008    0.01008    0.00427    0.00433
    stdTime      0.04222    0.06666    0.08672    0.50847    0.11745    0.01398    0.32337          0          0