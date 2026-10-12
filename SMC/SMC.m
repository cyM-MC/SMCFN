function [S, Time] = SMC(Smiss, noff, non)

% @param Smiss Initial similarity matrix
% @param noff Number of offline samples
% @param non Number of online samples
%
% @return S Calibrated similarity matrix
startTime = cputime; % Start the CPU timer
tolSMC = [];
 

Sobs = Smiss(1:noff, 1:noff);
Simp = Sobs;

tic;
for i = 1 : non
    v = Smiss(1:noff+i-1, noff+i);
    vonl = onlcal(Simp, v);
    Simp = [Simp, vonl; vonl', 1];
      
end

S = Simp;
elapsedTime = cputime - startTime; % Stop the timer and record the elapsed time
 Time = elapsedTime; % Stop the CPU timer and record the total CPU time
end

%% One-step OAMC
function [vonl] = onlcal(S0, v0)
%
% @param S0 Pairwise similarity matrix
% @param v0 Intial similarity vector
% @return vonl Calibrated similarity vector
tolSMC=[];
n = size(S0, 1);
tol = 1e-4;

[U, S, ~] = svd(S0);
s = diag(S);
C = U * diag(sqrt(s));
Cinv = diag(1./s) * C';

y0 = Cinv * v0;
if norm(y0) <= 1
    ycal = y0;
else
 
    lambda_min = max(sqrt(s'.^2 * y0.^2) - max(s), 0);
    lambda_max = sqrt(s'.^2 * y0.^2) - min(s);
    lambda = lambda_min;
    ycal = (s ./ (s+lambda)) .* y0;
    ylen = norm(ycal);
    iters = 0;
    %while (ylen > 1) || (ylen < 1-tol)
    ycalprev =ylen;
    while iters < 50 %maxiters
        lambda = 0.5*(lambda_min + lambda_max);
        ycal = (s ./ (s+lambda)) .* y0;
        ylen = norm(ycal);
        % use bisection search to find optimal lambda
        if ylen > 1
            lambda_min = lambda;
        elseif ylen < 1-tol
            lambda_max = lambda;
        end
        %         if norm(ylen-ycal,'inf')<=tol
        %             break;
        %         end

        iters = iters + 1;
    tolSMCt=norm(ycalprev-ycal,'fro');%/ norm(ycalprev,'fro');
    tolSMC = [tolSMC, tolSMCt]; % Store tolerance value

    ycalprev = ycal;
    end
    %     iter
    save('tolSMC.mat', 'tolSMC');
end
vonl = C * ycal;
end