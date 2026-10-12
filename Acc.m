function [simX, simY, acc] = Acc(X, Y, topk)
%Recall@topk
[numRows, ~] = size(X); % Assuming both matrices have the same size
overlapCounts = zeros(numRows, 1);

% Sort and get indices of topk indx for each matrix
[~,idxX] = sort(X,2,'descend');

[~,idxY] = sort(Y,2,'descend');

[m,~] = size(X);
simX = idxX(:,1:topk);
simY = idxY(:,1:topk);

for i = 1:length(idxX(:,1))
    overlapCounts(i) = length(intersect(simX(i,:),simY(i,:)));
end
count = sum(overlapCounts);
acc = count/(m*topk);
end