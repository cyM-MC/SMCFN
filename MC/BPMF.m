function [hatR, Time] = BPMF(Smiss, r, maxIter)
    % Smiss: Input matrix with missing values (NaNs)
    % r: Number of latent factors
    % maxIter: Number of iterations for the algorithm

    % Initialize parameters
    [M, N] = size(Smiss);
    U = abs(rand(M, r));
    V = abs(rand(N, r));
    
    hatRt = U * V';
    % Hyperparameters
    alpha = 1;
    beta = 0.01;
    
    % Initialize time variables
    Time = zeros(1, maxIter);
    startTime = cputime;
    
    % Main loop
    for iter = 1:maxIter
        hatR_prev = hatRt;
        
        % Update U
        U = U .* (((Smiss * V) + (alpha * U)) ./ ((U * (V' * V)) + (alpha * U)));
        
        % Update V
        V = V .* (((Smiss' * U) + (alpha * V)) ./ ((V * (U' * U)) + (beta * V)));
        hatR = U * V';
        hatRt = hatR;  
        
        % Calculate reconstruction error
        tolBPMFt = norm(hatR - hatR_prev, 'fro') / norm(hatR, 'fro');
        tolBPMF(iter) = tolBPMFt; % Store tolBPMFt in the tolBPMF array
        % Display iteration information
        %disp(['Iteration: ' num2str(iter) ', Reconstruction error: ' num2str(tolBPMFt)]);
    end 
    
    save('tolBPMF.mat', 'tolBPMF');
    Time = cputime - startTime;

    hatR = hatRt;
    hatR = mapminmax(hatR, 0, 1); 

end