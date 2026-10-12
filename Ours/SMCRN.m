%Alg2: SMCRN
function [hatS, rankS, Time] = SMCRN(Smiss, r, gamma, lambda, alpha, maxiters)
    % Initial V0 randomly
    tolgeSMC = [];

    tolhatgeSMC = [];
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
    
   % Measure memory before the operation 
    % GD process
    while iters < maxiters
        grad_F = CalculteGrad_F(Vt, Smiss, lambda, alpha); % Calculate the gradient of F
        
        V_prev = Vt; % Update V

        hatSprev = hatSt;
        Vt = Vt - gamma * grad_F;
        
        tol = norm(Vt - V_prev, 'fro'); % Calculate the norm difference of Vt vs. Vt-1
        tolgeSMCt = tol / norm(Vt, 'fro');
        tolgeSMC = [tolgeSMC, tolgeSMCt];
        
        hatSt =  Vt * Vt';
        hatSt = 2 * (hatSt - min(hatSt, [], 'all')) / range(hatSt, 'all') - 1;  
   
        tolhatSt = norm(hatSt - hatSprev, 'fro')/ norm(hatSt,'fro');
        tolhatgeSMC = [tolhatgeSMC, tolhatSt];

        rankS = rank(hatSt);
        iters = iters + 1; 
        
    end
     
    save('tolgeSMC.mat', 'tolgeSMC');
    
    save('tolhatgeSMC.mat', 'tolhatgeSMC');
    hatV = Vt; % Output the estimated hatV and hatS
    hatS = hatV * hatV'; 
    Time = cputime - startTime; % Record the total CPU time

    info = whos('hatV');  
end
% 

% Gradient Calculation Function 
function grad_F = CalculteGrad_F(V, S, lambda, alpha) 
    n = size(V, 1); 
    grad1 = 2 * (V * V' - S) * V; 
    grad2 = zeros(size(V));
    for i = 1:n
        Vi = V(i, :);
        norm_Vi = norm(Vi);
        if norm_Vi > alpha
            grad2(i,: ) = 4 * lambda * (norm_Vi)^3 * Vi / norm_Vi;
        end
    end 
    grad_F = grad1 + grad2;
end
 