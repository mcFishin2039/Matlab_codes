

function [smoothedEstimate] = movingAverageSmooth(x,y,xo,winsize)
    
    % initialize
    x = x(:);
    y = y(:);
    xo = xo(:);

    smoothedEstimate = zeros(size(xo)); % empty matrix
    
    % perform average smoothing 
    for n=1:length(xo)
        Ix=find(abs(x-xo(n)) <= winsize); % locate points within given window size
        smoothedEstimate(n) = mean(y(Ix)); % estimate the model at this location
    end


end 