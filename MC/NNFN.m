function [hatS, Time] = NNFN(Smiss, gamma, eta, maxiters)

    St = Smiss;
    tolNNFN = [];

    % initial S
    [n, ~] = size(Smiss);
    hatS = randn(n, n);
    iters = 0;

    startTime = cputime; % Start the CPU timer

    while iters < maxiters
        S_prev = St;
        Z = hatS - gamma * (2 * (hatS - Smiss));

        [U, Sigma, V] = svd(Z, 'econ');

        Sigma_prox = max(diag(Sigma) - eta, 0);

        hatS = U * diag(Sigma_prox) * V';
        iters = iters + 1;

        St = hatS;
        tolNNFNt = norm(St - S_prev, 'fro') / norm(St, 'fro');

        tolNNFN = [tolNNFN, tolNNFNt];
        
        save('tolNNFN.mat', 'tolNNFN');
    end

Time = cputime - startTime; % Stop the CPU timer and record the total CPU time
end
 