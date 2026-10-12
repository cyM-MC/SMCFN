function [S, Time] = CMC(S0, maxiter)
    if nargin < 2
        maxiter = 20;
    end

    r = 10; % number of partitions
    n = size(S0, 1);
    np = cell(r, 1);

    % Convergence
    S_prev = S0;
    tolCMC = [];
    
    startTime = cputime; % Start the CPU timer

    for t = 1 : maxiter
        tic; % Start the timer

        p = reshape(mod(randperm(2 * n), n) + 1, 2 * n / r, r);
        for i = 1 : r
            np{i} = unique(p(:, i));
        end
        St = cyclic_cal(S0, np, 2);
        
        tolCMCt = norm(St - S_prev, 'fro') / norm(St,'fro');
        tolCMC = [tolCMC, tolCMCt];
        S_prev = S0;
        S0 = St;
        
        elapsedTime = toc; % Stop the timer and record the elapsed time
    end

    save('tolCMC.mat', 'tolCMC'); % Save the tolerance values
    S = S0;
    Time = cputime - startTime; % Stop the CPU timer and record the total CPU time
end

%%
function [X, iter] = nearpsd(A, maxits, low, high) 
tolCMC = [];
if  ~isequal(A,A')
    A = (A + A') / 2;
end
if nargin < 4
    high = 1;
end
if nargin < 3
    low = -1; 
end
if nargin < 2
    maxits = 20;
end

% threshold for convergence & eigs
tolconv = 1.0e-4;
toleigs = 1.0e-5;

n = size(A,1);

U = zeros(n);
Y = A;

[V, D] = eig(Y);
d = diag(D);

iter = 0;
while 1
    T = Y - U;

    % project onto psd matrices
    [Q, D] = eig(T);
    d = diag(D);
    p = d > toleigs*d(n);
    X = Q(:,p) * D(p,p) * Q(:,p)';

    % update correction
    U = X - T;

    % maximum iteration & convergence test
    iter = iter + 1;
    if iter == maxits
        %fprintf('Max iterations reached. ');
        break; 
    end 

    if norm(Y-X,'inf')/norm(Y,'inf') <= tolconv 
    	break;
    end
     
    % problem-dependent here
    Y = X;
    Y(1:n+1:n*n) = 1;
    Y(Y<low) = low;
    Y(Y>high) = high;
end

Y(1:n+1:n*n) = 1;
Y(Y<low) = low;
Y(Y>high) = high;
 

end

%%
function [Jt] = cyclic_cal(J0, p, iter) 

if (nargin < 3)
    iter = 2;
end
r = length(p);
nc = zeros(r, 1);
It = cell(r,1);
for i = 1 : r
    nc(i) = length(p{i});
    It{i} = zeros(nc(i));
end

Jt = J0;
for t = 1 : iter 
    for i = 1 : r
            P = Jt(p{i}, p{i}) - It{i};
            Q = nearpsd(P, 15);
            Jt(p{i},p{i}) = Q;
            It{i} = Q - P;
    end
end
end

