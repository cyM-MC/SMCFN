%SMCFN
function [hatS, rankS, Time] = SMCFN(Smiss, r, gamma, lambda, maxiters)
    % Initial V0 randomly
    tolSSNN = []; 
    tolhatSSNN = [];
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
        gradient_F = CalculteGrad_F(Vt, Smiss, lambda); % Calculate the gradient of F
        
        V_prev = Vt; % Update V

        hatSprev = hatSt;
        Vt = Vt - gamma * gradient_F;
        
        tol = norm(Vt - V_prev, 'fro'); % Calculate the norm difference of Vt vs. Vt-1
        tolSSNNt = tol / norm(Vt, 'fro');
        tolSSNN = [tolSSNN, tolSSNNt];
        
        hatSt =  Vt * Vt';
        hatSt = 2 * (hatSt - min(hatSt, [], 'all')) / range(hatSt, 'all') - 1;
        tolhatSt = norm(hatSt - hatSprev, 'fro')/ norm(hatSt,'fro');
        tolhatSSNN = [tolhatSSNN, tolhatSt];
        
        rankS = rank(hatSt);
        iters = iters + 1; 
        
    end
    
    save('tolSSNN.mat', 'tolSSNN');
    
    save('tolhatSSNN.mat', 'tolhatSSNN');
    hatV = Vt; % Output the estimated hatV and hatS
    hatS = hatV * hatV'; 
    
    Time = cputime - startTime; % Record the total CPU time

    info = whos('hatV'); 
    mem = info.bytes;
end
 

% Gradient Calculation Function
function grad_F = CalculteGrad_F(V, S, lambda)
    grad_F = 4 * (V * V' - S) * V + 2 * lambda * V ;
end