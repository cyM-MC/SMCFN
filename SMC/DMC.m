function [S, Time] = DMC(Smiss, maxiter)
    if nargin < 2
        maxiter = 20;
    end

    low = -1;
    high = 1;
    S = Smiss; % Initialize S with S0 
    % Convergence
    S_prev = S;
   
    tolDMC = []; % Initialize the array to store tolerance values

    startTime = cputime; % Start the CPU timer

    for iter = 1:maxiter
        St = nearpsd(S, maxiter, low, high); % Update S in each iteration

        tolDMCt = norm(St - S_prev, 'fro') / norm(St,"fro"); % Calculate tolerance
        tolDMC = [tolDMC, tolDMCt]; % Store tolerance value 
         
        S_prev = S;

        elapsedTime = cputime - startTime; % Stop the timer and record the elapsed time
    end
    save('tolDMC.mat', 'tolDMC'); % Save the tolerance values
    S = Smiss;
    Time = elapsedTime; % Stop the CPU timer and record the total CPU time
end

%%
function [Y, iter] = nearpsd(A, maxits, low, high)

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
    if  norm(Y-X,'inf')/norm(Y,'inf') <= tolconv 
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
