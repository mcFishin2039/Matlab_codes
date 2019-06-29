% script to understand parametric and non-parametric statistical models, 
% cross-validation, and application to a 1-D ice-flow model
% 
% 
% Geostatistics and Data Analysis (GEOPH 522)
% Homework #3
% 
% Author: Marlon D. Ramos

clear; close all; clc

load ice_velocity.txt

z = ice_velocity(:,1); % ice depth, [meters]
v = ice_velocity(:,2); % velocity, [meters/second]
L = length(v); % number of data points

order = 0:1:4; % polynomial order vector 
l = length(order); % number of order values
pColors = {'g', 'r', 'y', 'c', 'b'}; % color plotting 

fig = figure(1); clf 

% initialize 
rms_error = zeros(1,l); 

%--------------------------------------------------------------
% PARAMETRIC APPROACH
% fitting the data to n-order polynomial fit
for ii = 1:l
    
    % generate the polymonial coefficients of a specified order
    p = polyfit(z, v, order(ii) ); 
    
    % store the polynomial coefficients using a cell array
    C{ii} = p; 
    
    % evaluate the polynomial data at points of the indepenent variable
    Vm = polyval(p,z);  
                      
    % calculate the root mean square error
    rms_error(ii) = rmse(L,v,Vm);
                       
    % plot the polynomial fit               
    plot(z, Vm, 'Color', pColors{ii} , 'LineWidth', 3); hold on 
    legendInfo{ii} = ['RMSE = ', num2str(rms_error(ii))];
    
end 

% append legend
legend(legendInfo, 'Location', 'SouthWest')

% plot the original ice velocity data as a function of depth 
plot(z, v, 'o', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k'); grid on 
title('parametric parameter fit '); xlabel('velocity [meters/year]'); ylabel('depth [meters]');
fig.Color = 'w'; fig.Position = [100 400 700 500];
set( findall( fig, '-property', 'FontSize' ), 'FontSize', 18 );

% assessing uncertainty of the model parameters ---------------------------

d_90 = round(L*90/100.0, 0); % 90 percent of the data samples
nIter = 1000;

% initialize arrays
RMS = zeros(nIter,l); % root mean square

% for every iteration
for jj = 1:nIter

    % for every polynomial order
    for kk = 1:l
     
    % randomly sample a 90% of the data
    random_z = randsample(z, d_90); % depths
    random_v = randsample(v, d_90); % velocities
    p = polyfit(random_z, random_v, order(kk) );  % n-order polymial fit
    polyCoeff{jj,kk} = p; % store
    Vmm = polyval(p, random_z); % model estimate

    end 

end 

% extract polynomial coefficients and assess their mean 
% and standard deviation
a0 = [polyCoeff{1:end,1}]; 
a1 = [polyCoeff{1:end,2}]; 
a2 = [polyCoeff{1:end,3}]; 
a3 = [polyCoeff{1:end,4}]; 
a4 = [polyCoeff{1:end,5}]; 

% average values for each parameter coefficient
m1 = mean(a0); m2 = mean(a1); m3 = mean(a2); m4 = mean(a3); m5 = mean(a4);
means = [m1 m2 m3 m4 m5];

% standard deviations
s1 = std(a0); s2 = std(a1); s3 = std(a2); s4 = std(a3); s5 = std(a4);
standard_dev = [s1 s2 s3 s4 s5];

%%
% cross validation 

windowSize = 10;

RMSE = zeros(nIter,l); 

for jj = 1:nIter % for each iteration
    for n = 1:l % for each polynomial fit

        Ix2 = randsample(L, d_90); % randomly sample 90% of the data

        % training set
        zt = z(Ix2); % temporary depth vector
        vt = v(Ix2); % temporary velocity vector

        % testing set
        Ix3 = ones(size(z)); 
        Ix3(Ix2) = 0; 
        Ix3 = find(Ix3 == 1); % test set - the 10% leftover
        ztest = z(Ix3); 
        vtest = v(Ix3);

        % we want to sample at values where we did not already train with
        vmod = nonparametric_smooth(zt,vt,ztest, windowSize); 
        RMSE(jj,n) = sqrt(mean((vmod-vtest).^2));

    end 
    
end 

rmse1 = RMSE(:,1); % n = 0
rmse2 = RMSE(:,2); % n = 1
rmse3 = RMSE(:,3); % n = 2
rmse4 = RMSE(:,4); % n = 3
rmse5 = RMSE(:,5); % n = 4


% plot the distribution of rmse errors
fig2 = figure(2); clf 
numBars = 12;
subplot(511)
hist(rmse1, numBars); xlabel('rmse values for n = 0'); ylabel('frequency');

subplot(512) 
hist(rmse2, numBars); xlabel('rmse values for n = 1'); ylabel('frequency');

subplot(513)
hist(rmse3, numBars); xlabel('rmse values for n = 2'); ylabel('frequency'); 

subplot(514)
hist(rmse4, numBars); xlabel('rmse values for n = 3'); ylabel('frequency');

subplot(515)
hist(rmse5, numBars); xlabel('rmse values for n = 4');ylabel('frequency');
fig2.Color = 'w'; fig2.Position = [500 90 800 900];

%%
%--------------------------------------------------------------
% NON-PARAMETRIC APPROACH
% generate a moving window of some sample size
% we remove the assumption that our data fits a curve 
% at each depth range, you are providing a best estimate for the 
% velocity. As seen below, we use the mean as the best estimate. 
% a more elegant approach to this problem would be to have weight the data 
% closer to the middle of the window more than the data further away from
% the window

windows = [3 10 50]; % window lengths [meters]

% intialize
Y = zeros(length(windows), length(z) );

for index = 1:3
    Y(index,:) = movingAverageSmooth(z,v,z,windows(index)); % smoothed average estimate 
end 

fig3 = figure(3); clf
subplot(211)

p1 = plot(z,Y(1,:), 'LineWidth', 2); hold on 
p2 = plot(z,Y(2,:), 'LineWidth', 2);
p3 = plot(z,Y(3,:), 'LineWidth', 2);
original = plot(z,v,'o'); % original data 

grid on
title('moving window average of ice velocity'); 
xlabel('velocity [meters/year]'); ylabel('depth [meters]');
legend([p1, p2, p3, original], 'window size = 3 m', 'window size = 10 m', 'window size = 50 m',...
                                    'data', 'Location', 'Southwest')

% applying a WEIGHTED MOVING WINDOW over the same window lengths    
Y2 = zeros(length(windows), length(z) ); % initialize

for index = 1:3
    Y2(index,:) = nonparametric_smooth(z,v,z,windows(index)); % smoothed average estimate 
end 

subplot(212)

pp1 = plot(z,Y2(1,:), 'LineWidth', 2); hold on 
pp2 = plot(z,Y2(2,:), 'LineWidth', 2);
pp3 = plot(z,Y2(3,:), 'LineWidth', 2);
original = plot(z,v,'o'); % original data 

grid on
title('moving, weigted window average of ice velocity'); 
xlabel('velocity [meters/year]'); ylabel('depth [meters]');
legend([pp1, pp2, pp3, original], 'window size = 3 m', 'window size = 10 m', 'window size = 50 m',...
                                    'data', 'Location', 'Southwest')

fig3.Color = 'w'; fig3.Position = [500 90 700 700];
%% use a brute-force grid search technique to find the model parameters that minimize 
%%  the rms misfit
%%
% surface velocity u(z=0)
u_surf = v(1);

% choose an ensemble of potential paramter values for A 
A = 1e-18:0.1e-18:10e-18;  
A_length = length(A); 

% and for n 
nn = 2:0.01:4;
n_length = length(nn); 

% initialize
ICE_RMS = zeros(A_length, n_length);

% iterate over the values of A 
for indexA = 1:A_length

    % iterate over the values of power n
    for indexn = 1:n_length
    
        % calculate the expected ice velocity and the root mean square
        % error
        ICE_RMS(indexA, indexn) = icemodel_rms(v, z, A(indexA), nn(indexn), u_surf);
    
    % exit n loop 
    end 
% exit A loop
end 

% visualize the RMS error
clc
fig4 = figure(4); clf
colormap('jet')
imagesc(A, nn, ICE_RMS, [0 500]); 
cc = colorbar; 
ylabel(cc, 'root mean squared error')
ylabel('n values'); xlabel('A values');
title('RMSE(depth) - grid search'); 
fig4.Color = 'w'; fig4.Position = [500 90 700 600];
%%

% repeat applying the fminsearch function
% choose model parameters A and n contained in vector P 
fh=@(P)icemodel_rmse(v,z,u_surf,P); % create function handle 

% initialize
rmse_output = zeros(A_length, n_length);

tic 
% iterate over the values of A 
for indexA = 1:A_length

    % iterate over the values of power n
    for indexn = 1:n_length
        
        point = [A(indexA) nn(indexn) ];
        [Pbest,fval] = fminsearch(fh, point); % look for the minimum of the function evaluated at a starting point
        
        rmse_output(indexA, indexn) = fval;
        
    end
end 
toc

fig5 = figure(5); clf 
colormap('jet')
imagesc(A, nn, rmse_output); 
cbar = colorbar; ylabel(cbar, 'root mean squared error'); 
title('RMSE(depth) - gradient method'); 
ylabel('n values'); xlabel('A values');
fig5.Color = 'w'; fig5.Position = [400 90 700 600];

%%
%--------------------------------------------------
% randomly sample 90% of the dataset and find the optimum value of A using
% the gradient search method repeated nIter = 1000 times. 

constant_n = 3; % from theory

% initialize
rmse_output2 = zeros(nIter, 1);
bestA = []; % optimum values of A 

% pull out 90% of the data at random
random_z2 = randsample(z, d_90); % depths
random_v2 = randsample(v, d_90); % velocities [m/yr]

% modify function handle
fhh=@(P)icemodel_rmse(random_v2, random_v2, u_surf, P);


for jx = 1:nIter
           
        Astart = randsample(A,1); % choose a random value
        P = [ Astart constant_n ]; % starting point to evaluate
        [Pbest2,fval2] = fminsearch(fhh, P); % find values that achieve the minimum value
        rmse_output2(jx) = fval2; % rmse 
        
        % store
        bestA = [bestA; Pbest2(1) ];
    
end 

% plot the pdf (relative density histogram) of the distribution of A and
% the rms error
%%
figure(6); clf 
[Hist1, Xpos1] = hist(bestA); 
[Hist2, Xpos2] = hist(rmse_output2);

subplot(211)
hist(bestA) 
title('histogram of A')

subplot(212)
hist(rmse_output2)
title('histogram of RMS error')

%%
% 












