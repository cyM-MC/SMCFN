% variation 2 SMC_NR
% S = VV^T (no norm)
function [hatS, hatV, totalItes, clockTime, cpuTime] = ab2(Smiss, r, gamma, maxiters)
    % Initial V0 randomly
    [n, ~] = size(Smiss);
    V0 = randn(n, r); 
    
    % Paras initialization
    Vt = V0;
    V_prev = Vt;
    tol = Inf;
    iters = 0;
    
    % Time
    clockStart = tic; % Clock time
    [~, cpuStart] = system('echo $SECONDS'); % CPU time
    
    % GD process
%    while tol > 1e-4 && iters < 100
       %Vt = normMatrix(Vt);

     while iters < maxiters
        gradient_F = CalculteGrad_F(Vt, Smiss); % Calculate the gradient of F
        
        V_prev = Vt; % Update V

        Vt = Vt  - gamma * gradient_F;  
        
        tol = norm(Vt - V_prev, 'fro'); %;/norm(Vt,'fro'); % Calculate the norm difference of Vt vs. Vt-1
        
        iters = iters + 1;
    end
    %Vt = normMatrix(Vt);
    hatV = Vt; % Output the estimated hatV and hatS
    hatS = hatV * hatV';
    
    clockTime = toc(clockStart); % Clock time
   
    [~, cpuEnd] = system('echo $SECONDS'); % CPU time
    cpuTime = str2double(cpuEnd) - str2double(cpuStart);
    totalItes = iters;
end

% Gradient Calculation Function
function grad_F = CalculteGrad_F(V, S)
    grad_F = 2 * (V * V' - S) * V ;
end