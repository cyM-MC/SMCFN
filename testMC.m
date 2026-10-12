%%
%clear all; clc;
warning off;
addpath("MC\"); 

%random toy data
load toydata.mat;
 
 
%% Data Normalization
X = double(X);
n = size(X,2);
% normalize to [-1, 1]
X = X - repmat(min(X,[],2), 1, n);
X = X ./ (repmat(max(X,[],2), 1, n) + eps);
X = (X - 0.5) * 2;

%% Parameter Setting
p_list = [0.7,0.8,0.9]; % vxalues of missing ratio p 
niter = 5;                % number of it erations
noff = 1000;               % number of search candidates
non = 200;                % number of query items
seed = 42;              % random seed  

% paras
LowR = 100;  lambda= 0.001; gamma = 0.001; 
alpha = 0.1; 
maxiters = 1000;
%recall 
topk = floor(0.2*(non+noff));

%mask M for BPMF,FGSR
M = ones(noff + non);
M(1:noff,1:noff) = 0; 

for k = 1 : length(p_list)
    p = p_list(k);
    rng(seed);

    %% Similarity Search with Missing Data
    for i = 1 : niter
        fprintf('\nIters = %1.0f: ', i);

        %% Data Construction - introduce random missing values into a dataset
        rp = randperm(n);
        Xoff = X(:, rp(1:noff));
        Xon = X(:, rp(noff+1 : noff+non));
        Xtrue = [Xoff, Xon];
        Strue  = similarity(Xtrue);

        Xon(rand(size(Xon)) < p) = NaN;
        Xmiss = [Xoff, Xon];
        Smiss = similarity(Xmiss, 'miss');
        Fnorm = norm(Smiss-Strue, 'fro')^2;
        FnormmSpp = norm(Smiss(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2;
 
        %% PSD  
        %% CMC 
        [hatSCMC,TimeCMC] = CMC(Smiss); 
        RMSE(i,1) = norm(hatSCMC(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        RecallC(i,1) = Acc(hatSCMC, Strue,topk);
        nDCGC(i,1) = nDCG(Strue, hatSCMC, topk);
        Time(i,1) = TimeCMC;

        %% DMC 
        [hatSDMC,TimeDMC] = DMC(Smiss); 
        RMSE(i,2) = norm(hatSDMC(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        RecallC(i,2) = Acc(hatSDMC, Strue, topk);
        nDCGC(i,2) = nDCG(Strue, hatSDMC, topk);
        Time(i,2) = TimeDMC;
       
        %% FMC 
        [hatSSMC, TimeFMC] = SMC(Smiss, noff, non);
        RMSE(i,3) = norm(hatSSMC(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        RecallC(i,3) = Acc(hatSSMC, Strue, topk);
        nDCGC(i,3) = nDCG(Strue, hatSSMC, topk);    
        Time(i,3) = TimeFMC;
       
        %% Low Rank
        %% WWW Alg1
        [hatSNNFN, TimeNNFN] = NNFN(Smiss, gamma, lambda, maxiters); 
        RMSE(i,4) = norm(hatSNNFN(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        RecallC(i,4) = Acc(hatSNNFN, Strue,topk);
        nDCGC(i,4) = nDCG(Strue, hatSNNFN, topk);
        Time(i,4) = TimeNNFN;
       
        %% WWW Alg2
        [hatSfacNNFN, TimefNNFN] = facNNFN(Smiss, LowR, gamma, lambda, maxiters); 
        RMSE(i,5) = norm(hatSfacNNFN(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        RecallC(i,5) = Acc(hatSfacNNFN, Strue,topk);
        nDCGC(i,5) = nDCG(Strue, hatSfacNNFN, topk);
        Time(i,5) = TimefNNFN;
       
        %% FGSR
        [hatFGSR, TimeFGSR] = FGSR(Smiss, M, LowR, 1, 0.001);  
        RMSE(i,6) = norm(hatFGSR(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        RecallC(i,6) = Acc(hatFGSR, Strue, topk);
        nDCGC(i,6) = nDCG(Strue, hatFGSR, topk);
        Time(i,6) = TimeFGSR;
       
        %% BPMF
        [hatSBPMF,TimeBPMF] = BPMF(Smiss, LowR, maxiters); 
        RMSE(i,7) = norm(hatSBPMF(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        RecallC(i,7) = Acc(hatSBPMF, Strue, topk);
        nDCGC(i,7) = nDCG(Strue, hatSBPMF, topk);
        Time(i,7) = TimeBPMF;

        %% Our PSD + Low Rank + Nuclear Norm
        %% SMCFN
        [hatSSSNN, TimeSMCFN] = SMCFN(Smiss, LowR, gamma, lambda, maxiters);
        RMSE(i,8) = norm(hatSSSNN(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        RecallC(i,8) = Acc(hatSSSNN, Strue, topk);
        nDCGC(i,8) = nDCG(Strue, hatSSSNN, topk);
        Time(i,8) = TimeSMCFN;
        %% SMCRN
        [hatSSMCRN, TimeSMCRN] = SMCRN(Smiss, LowR, gamma, lambda, alpha, maxiters); 
        RMSE(i,9) = norm(hatSSMCRN(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        RecallC(i,9) = Acc(hatSSMCRN, Strue, topk);
        nDCGC(i,9) = nDCG(Strue, hatSSMCRN, topk); 
        Time(i,9) = TimeSMCRN;
        
        fprintf('Finish.');
    end

    %%  fprintf(': noff=%1.0f, non=%1.0f, p=%1.2f, niter=%1.0f\n', noff, non, p, niter);
    fprintf(' r=%1.0f, search=%1.0f, query=%1.0f, p=%1.2f, niter=%1.0f', LowR, noff, non, p, niter);
    fprintf('\n'); 
    stat = [mean(RMSE);  mean(RecallC); mean(nDCGC); mean(Time); std(RMSE);std(RecallC);std(nDCGC);std(Time);];
    Stat = roundn(stat, -5);
    Table = table(Stat(:,1), Stat(:,2),Stat(:,3),Stat(:,4),Stat(:,5),Stat(:,6),Stat(:,7),Stat(:,8),Stat(:,9), ...
        'VariableNames',{'CMC','DMC','FMC','NNFN','facNNFN','FGSR','BPMF','SMCFN','SMCRN',},...
        'RowNames',{'RMSE';'Recall';'nDCG';'Time';'stdRMSE';'stdRecall';'stdnDCG';'stdTime'});
    disp(Table)
end
