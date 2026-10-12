function ndcg = nDCG(Strue, hatS, topk)
% nDCG@topk 

[numRows, ~] = size(Strue);
DCG = zeros(numRows, 1);
IDCG = zeros(numRows, 1);

for i = 1:numRows 
    % Calculate estimated ranking
    [~, predictedIndices] = sort(hatS(i,:), 'descend');
    predictedRelevance = Strue(i, predictedIndices);
    
    %  DCG
    for k = 1:topk
        DCG(i) = DCG(i) + (2^predictedRelevance(k) - 1) / log2(k+1);
    end
    
    % Ideal ranking
    [~, idealIndices] = sort(Strue(i,:), 'descend');
    idealRelevance = Strue(i, idealIndices);
    
    % Calculate IDCG
    for k = 1:topk
        IDCG(i) = IDCG(i) + (2^idealRelevance(k) - 1) / log2(k+1);
    end
    
    % avoid zero
    if IDCG(i) == 0
        DCG(i) = 0;
    else
        DCG(i) = DCG(i) / IDCG(i);
    end
end

% mean nDCG
ndcg = mean(DCG);
end
