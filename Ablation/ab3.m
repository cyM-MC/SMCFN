%ab3 SMCGD
function [hatS, Time] = ab3(Smiss, r, gamma, lambda, maxiters)
% Initial V0 randomly
tolhatVab3 = [];

tolhatSab3 = [];

[n, ~] = size(Smiss);
V0 = randn(n, r);

% Paras initialization
Vt = V0;
V_prev = Vt;
tol = Inf;
iters = 0;
hatSt = Vt * Vt';
% Initialize time variables
Time = zeros(1, maxiters);
startTime = cputime; 


% GD process
while iters < maxiters
    gradient_F = CalculteGrad_F(Vt, Smiss, lambda); % Calculate the gradient of F

    V_prev = Vt; % Update V
    hatSprev = hatSt; 
    Vt = Vt - gamma * gradient_F;

    tol = norm(Vt - V_prev, 'fro'); % Calculate the norm difference of Vt vs. Vt-1
    ab3 = tol / norm(Vt, 'fro');
    tolhatVab3 = [tolhatVab3, ab3];

    hatSt =  Vt * Vt'; 

    [A, D] = eig(hatSt);
    d = diag(D);
    p = d>1e-5*d(n);
    hatSt_psd = A(:,p) * D(p,p) * A(:,p)';

    Vt = cholcov(hatSt_psd,'lower');
    Vt = Vt';
    Vt = Vt(:,1:r);
    
    hatSt = Vt * Vt';
    
    tolhatSt = norm(hatSt - hatSprev, 'fro');
    tolhatSab3 = [tolhatSab3, tolhatSt];

    iters = iters + 1;

end

save('ab3.mat', 'tolhatVab3');
save('tolab3.mat', 'tolhatSab3');

hatV = Vt; % Output the estimated hatV and hatS
hatS = hatV * hatV';

Time = cputime - startTime; % Record the total CPU time
 

end
 

% Gradient Calculation Function
function grad_F = CalculteGrad_F(V, S, lambda)  
    grad_F = 2 * (V*V' - S) * V + lambda * V * ((V'*V)^(-0.5));
end
