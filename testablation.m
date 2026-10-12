%%
%clear all; clc;
warning off;

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
seed = 32;              % random seed  

% paras
LowR = 100;  lambda= 0.001; gamma = 0.001; maxiters = 1000;
%recall 
topk = 0.2*(non+noff);
% rand V(fixed)
Vinitial = randn(non+noff, LowR);

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

        %% Nuclear norm
        %% ab1
        [hatSv1] = ab1(Smiss, LowR, gamma, lambda, maxiters); 
        RMSE(i,1) = norm(hatSv1(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        recall(i,1) = Acc(hatSv1, Strue, topk);
        ndcg(i,1) = nDCG(Strue, hatSv1, topk);
        
        %% ab2
        [hatSv2] = ab2(Smiss, LowR, gamma, maxiters); 
        RMSE(i,2) = norm(hatSv2(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        recall(i,2) = Acc(hatSv2, Strue, topk);
        ndcg(i,2) = nDCG(Strue, hatSv2, topk);
        
        %% ab3
        [hatSv3] = ab3(Smiss, LowR, gamma, lambda, maxiters); 
        RMSE(i,3) = norm(hatSv3(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        recall(i,3) = Acc(hatSv3, Strue,topk);
        ndcg(i,3) = nDCG(Strue, hatSv3, topk);
        
        %% ab4
        [hatSv4] = ab4(Smiss, gamma, lambda, maxiters) ; 
        RMSE(i,4) = norm(hatSv4(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        recall(i,4) = Acc(hatSv4,Strue, topk); 
        ndcg(i,4) = nDCG(Strue, hatSv4, topk);
 
        %% ab5
        [hatSv5] = ab5(Smiss, LowR, gamma, lambda, maxiters) ; 
        RMSE(i,5) = norm(hatSv5(noff+1:noff+non,1:noff+non) - Strue(noff+1:noff+non,1:noff+non), 'fro')^2 / FnormmSpp;
        recall(i,5) = Acc(hatSv5, Strue,topk);
        ndcg(i,5) = nDCG(Strue, hatSv5, topk); 
      
        fprintf('Finish.');
    end

    %%  fprintf(': noff=%1.0f, non=%1.0f, p=%1.2f, niter=%1.0f\n', noff, non, p, niter);
    %fprintf(' LowR=%1.0f, search=%1.0f, query=%1.0f, p=%1.2f, niter=%1.0f', LowR, noff, non, p, niter);
    fprintf('\n'); 
 
    stat = [mean(RMSE); mean(recall); mean(ndcg); std(RMSE); std(recall); std(ndcg)];
    Stat = roundn(stat, -5);
    Table = table(Stat(:,1), Stat(:,2), Stat(:,3), Stat(:,4), Stat(:,5), ...
    'VariableNames', {'ab1', 'ab2', 'ab3', 'ab4', 'ab5'}, ...
    'RowNames', {'RMSE', 'Recall', 'nDCG', 'stdRMSE', 'stdRecall', 'stdnDCG'});

disp(Table)
end
