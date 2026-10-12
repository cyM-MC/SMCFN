%% test for MC
% clear all; clc;
warning off;
addpath(genpath('MVI'));   
 
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
p_list = [0.5,0.6,0.7,0.8,0.9]; % vxalues of missing ratio p
niter = 5;                % number of it erations
noff = 1000;               % number of offline samples
non = 200;                % number of online samples
seed = 42;              % random seed for reproducing results

% paras
r = 100; gamma = 0.001; lambda = 0.001; maxiters = 1000;

% recall
topk = 0.2*(non+noff);

%mask M for BPMF,FGSR
M = zeros(noff+non);
M(1:noff,1:noff) = 1;

for k = 1 : length(p_list)
    p = p_list(k);
    rng(seed);

    %% Online Scenario
    for i = 1 : niter
        fprintf('\nIters = %1.0f: ', i);

        %% Data Construction -introduce random missing values into a dataset
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
 
        %% MVI
        %% MEAN
        [hatMEAN] = MEAN(Xoff, Xon, 'on'); 
        hatSMEAN = similarity(hatMEAN); 
        cliphatMEAN = min(hatSMEAN, 1);
        RMSE(i,1) = norm(cliphatMEAN(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        recall(i,1) = Acc(cliphatMEAN, Strue,topk);
        ndcg(i,1) = nDCG(Strue, cliphatMEAN, topk); 
        
        %% LR
        [hatLR] = LR(Xoff, Xon, 'on'); 
        hatSLR = similarity(hatLR);
        cliphatSLR = min(hatSLR, 1);
        RMSE(i,2) = norm(cliphatSLR(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        recall(i,2) = Acc(cliphatSLR, Strue,topk);
        ndcg(i,2) = nDCG(Strue, cliphatSLR, topk); 
        
        %% GR 
        [hatGR] = GR(Xoff, Xon, 'off'); 
        hatSGR = similarity(hatGR);
        cliphatSGR = min(hatSGR,1);
        RMSE(i,3) =  norm(cliphatSGR(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        recall(i,3) = Acc(cliphatSGR, Strue,topk); 
        ndcg(i,3) = nDCG(Strue, cliphatSGR, topk); 
        
        %% KFMC 
        [hatKFMC] = KFMC(Xmiss, Xtrue, 'on', 'rbf');  
        hatSKFMC = similarity(hatKFMC);
        cliphatSKFMC = min(hatSKFMC,1);
        RMSE(i,4) =  norm(cliphatSKFMC(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        recall(i,4) = Acc(cliphatSKFMC, Strue,topk);
        ndcg(i,4) = nDCG(Strue, cliphatSKFMC, topk);  

        fprintf('Finish.');
    end

    %%  fprintf(': noff=%1.0f, non=%1.0f, p=%1.2f, niter=%1.0f\n', noff, non, p, niter);
    fprintf(' noff=%1.0f, non=%1.0f, p=%1.2f, niter=%1.0f', noff, non, p, niter); 
    fprintf('\n');
    stat = [mean(RMSE);mean(recall); mean(ndcg); std(RMSE); std(recall); std(ndcg)];
    Stat = roundn(stat, -4);
    Table = table(Stat(:,1), Stat(:,2),Stat(:,3),Stat(:,4),...
        'VariableNames',{'MEAN','LR','GR','KFMC'}, ...
        'RowNames',{'RMSE';'Recall';'nDCG';'stdRMSE';'stdRecall';'stdnDCG'});
    disp(Table)
end




