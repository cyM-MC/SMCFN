function [imputedX, time] = LR(Xoff, Xon, model) 
if (nargin < 3)
    model = 'on';
end
low = -1; high = 1;

n_off = size(Xoff, 2);
n_on = size(Xon, 2);
imputedX = [Xoff, Xon];
if mod(n_on, 100) == 0
    time = zeros(1, n_on/100);
end

if strcmp(model, 'on')
    tic;
    for i = 1 : n_on
        y = Xon(:,i);
        idx = isnan(y);
        y_obs = y(~idx);
        X_mis = Xoff(idx,:);
        X_obs = Xoff(~idx,:);
        beta = regress(y_obs, [ones(sum(~idx),1), X_obs]);
        y_mis = sum(beta' .* [ones(sum(idx),1), X_mis], 2);
        y(idx) = y_mis;
        imputedX(:, n_off+i) = y;
        if mod(i, 100) == 0
            time(1, i/100) = toc;
        end
    end
elseif strcmp(model, 'seq')
    tic;
    for i = 1 : n_on
        y = Xon(:,i);
        idx = isnan(y);
        y_obs = y(~idx);
        X_mis = Xoff(idx,:);
        X_obs = Xoff(~idx,:);
        beta = regress(y_obs, [ones(sum(~idx),1), X_obs]);
        y_mis = sum(beta' .* [ones(sum(idx),1), X_mis], 2);
        y(idx) = y_mis;
        imputedX(:, n_off+i) = y;
        Xoff = [Xoff, y];
        if mod(i, 100) == 0
            time(1, i/100) = toc;
        end
    end    
elseif strcmp(model, 'off')
    fprintf('No Offline Version!');
end
imputedX(imputedX<low) = low;
imputedX(imputedX>high) = high;
end