function S = similarity(X, choice)
% Approximate a similarity matrix for samples with NaN values.
%
% @param X       d*n, each column is a sample
% @param choice  Default "true" (calculate the true similarity matrix)
% @return S      n*n

if (nargin < 2)
    choice = 'true';
end

if strcmp(choice, 'miss')
    [~, n] = size(X);
    Idx = isnan(X);
    Xzero = X; Xzero(Idx) = 0;
    XX = Xzero' * Xzero;
    Norm = zeros(n);
    Norm_ini = sqrt(sum(Xzero.^2, 1));
    for i = 1 : n
        idx = Idx(:, i);
        if sum(idx) == 0
            Norm(i, :) = Norm_ini;
        else
            idx = repmat(idx, 1, n);
            XI = Xzero; XI(idx) = 0;
            Norm(i, :) = sqrt(sum(XI.^2, 1));
        end
    end
    S = XX ./ (Norm .* Norm');
    S(isnan(S)) = 0;
elseif strcmp(choice, 'true')
    Norm = sqrt(sum(X.^2, 1));
    S = (X' * X) ./ (Norm' * Norm);
elseif strcmp(choice, 'basic')
    [d, n] = size(X);
    O =~isnan(X);
    S = zeros(n);
    for i = 1 : n
        for j = i + 1 : n
            k = O(:,i) & O(:,j);
            S(i,j) = X(k,i)' * X(k,j) / (norm(X(k,i))*norm(X(k,j))) ;
        end
    end
    S = S + S' + eye(n);
    S(isnan(S)) = 0;
end
end