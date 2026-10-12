function [imputedX, time] = MEAN(Xoff, Xon, model) 
if (nargin < 3)
    model = 'on';
end
low = -1; high = 1;

n_off = size(Xoff, 2);
n_on = size(Xon, 2);
if mod(n_on, 100) == 0
    time = zeros(1, n_on/100);
end

if strcmp(model, 'on')
    tic;
    imputedX = [Xoff, Xon];
    rowmean = mean(Xoff, 2);
    for i = 1 : n_on
        x = Xon(:, i);
        idx = isnan(x);
        imputedX(idx, n_off+i) = rowmean(idx);
        if mod(i, 100) == 0
            time(1, i/100) = toc;
        end
    end
elseif strcmp(model, 'seq')
    tic;
    for i = 1 : n_on
        rowmean = mean(Xoff, 2);
        x = Xon(:, i);
        idx = isnan(x);
        x(idx) = rowmean(idx);
        Xoff = [Xoff, x];
        if mod(i, 100) == 0
            time(1, i/100) = toc;
        end
    end
    imputedX = Xoff;
elseif strcmp(model, 'off')
    tic;
    Ximp = [Xoff, Xon];
    idx = isnan(Ximp);
    imputedX = Ximp;
    imputedX(idx) = 0;
    imputedX = imputedX + nanmean(Ximp,2).*idx;
    time = toc;
end
imputedX(imputedX<low) = low;
imputedX(imputedX>high) = high;
end