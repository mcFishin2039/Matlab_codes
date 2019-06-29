% calculate the root mean square error for a 1-D dataset 

% input: number of data points (N), data vector of true model (y), 
%        data vector of observed data (x)
% note that the length of vector y and x must be the same

function [RMSE] = rmse(N, x, y)
    
    % initialize 
    squareDiff = zeros(1,N);

    % calculate the square of the difference 
    for ii = 1:N
        squareDiff(ii) = ( y(ii) - x(ii) )^2;
    end 

    % calculate mean square error 
    mse = (1/N)*sum(squareDiff);
    
    % calculate root mean square error
    RMSE = sqrt(mse);
    
end 

