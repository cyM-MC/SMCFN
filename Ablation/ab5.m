%ab5 SMCfON
%min S = lambda/2 (W_F+H_F)
function [hatS, Time] = ab5(S0, r, gamma, lambda, maxiters)
    % Initial V0 randomly
    [n, ~] = size(S0);
    W0 = randn(n, r); 
    H0 = randn(n, r);

    % Paras initialization
    Wt = W0; W_prev = Wt;
    Ht = H0; H_prev = Ht;
    St = S0; S_prev = St;

    tolW = Inf;
    tolH = Inf; 
    iters = 0;

    tolab5 = [];
    
    % Initialize time variables
    Time = 0;
    startTime = cputime;

    % GD process
    while iters < maxiters
        S_prev = St;
        gradient_FW = CalculateGrad_FW(Wt, Ht, S0, lambda); % Calculate the gradient of F
        gradient_FH = CalculateGrad_FH(Wt, Ht, S0, lambda);

        W_prev = Wt; % Update V
        Wt = Wt - gamma * gradient_FW;

        H_prev = Ht; % Update H
        Ht = Ht - gamma * gradient_FH;

        tolW = norm(Wt - W_prev, 'fro'); % Calculate the norm difference of Wt vs. Wt-1
        tolH = norm(Ht - H_prev, 'fro');
        St = Wt * Ht'; 
        tolSt = norm(St - S_prev, 'fro') / norm(St, 'fro');
        iters = iters + 1;
        tolab5 = [tolab5, tolSt];
        
        % Calculate time taken for the iteration
        Time = cputime - startTime;
    end

    hatW = Wt; % Output the estimated hatW and hatH
    hatH = Ht;
    hatS = hatW * hatH';  

    save('tolab5.mat', 'tolab5');
end

% Gradient Calculation Function
function grad_FW = CalculateGrad_FW(W, H, S, lambda)
    grad_FW = lambda* (W * H' - 2 * S) * H ;
end

function grad_FH = CalculateGrad_FH(W, H, S, lambda)
    grad_FH = lambda* (W * H' - 2 * S)' * W ;
end