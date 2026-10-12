function [imputedX, time] = GR(Xoff, Xon, model) 

if (nargin < 3)
    model = 'on';
end
low = -1; high = 1;

n_off = size(Xoff, 2);
n_on = size(Xon, 2);
imputedX = [Xoff, Xon];
if mod(n_on, 100) == 0
    time = zeros(1, n_on/100);
end

max_rank = 10; 
step_size = 0.1;
max_Cycles = 50;

if strcmp(model, 'on')
    tic;
    for i = 1 : n_on
        Xmiss = [Xoff, Xon(:,i)];
        M = double(~isnan(Xmiss));
        [numr, numc] = size(Xmiss);
        [I, J, S] = miss_grouse(Xmiss, M);
        [U,V,~] = grouse(I,J,S,numr,numc,max_rank,step_size,max_Cycles);
        Ximp = U*V';
        imputedX(:, n_off+i) = Ximp(:, end);
        if mod(i, 100) == 0
            time(1, i/100) = toc;
        end
    end
elseif strcmp(model, 'seq')
    tic;
    for i = 1 : n_on
        Xmiss = [Xoff, Xon(:,i)];
        M = double(~isnan(Xmiss));
        [numr, numc] = size(Xmiss);
        [I, J, S] = miss_grouse(Xmiss, M);
        [U,V,~] = grouse(I,J,S,numr,numc,max_rank,step_size,max_Cycles);
        Ximp = U*V';
        imputedX(:, n_off+i) = Ximp(:, end);
        Xoff = [Xoff, Ximp(:, end)];
        if mod(i, 100) == 0
            time(1, i/100) = toc;
        end
    end    
elseif strcmp(model, 'off')
    tic;
    Xmiss = [Xoff, Xon];
    Xzero = impute_zero(Xmiss);
    M = double(~isnan(Xmiss));
    [numr, numc] = size(Xmiss);
    [I, J, S] = miss_grouse(Xmiss, M);
    [U,V,~] = grouse(I,J,S,numr,numc,max_rank,step_size,max_Cycles);
    imputedX = U*V';
    imputedX = Xzero + imputedX .* (1-M);
    time = toc;
end
imputedX(imputedX<low) = low;
imputedX(imputedX>high) = high;
end

%%
function [I, J, S] = miss_grouse(Xmiss, M)

I = []; J = []; S = []; kk = 0;
[numr, numc] = size(Xmiss);
for ii = 1:numc
    for jj = 1:numr
        if M(jj,ii) == 1
            kk = kk + 1;
            J(kk,1) = ii;
            I(kk,1) = jj;
            S(kk,1) = Xmiss(jj,ii);
        end
    end
end
end