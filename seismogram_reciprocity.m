% script to demonstrate GF reciprocity 
% Author: Marlon Ramos 
% modified from python notebook for computational seismology
clear; close all; clc; 
roma=load('~/ScientificColourMaps5/roma/roma.mat'); 

nt    = 300;                                   % number of time steps
c0    = 1.0;                                   % acoustic velocity, [m/s]
eps   = 0.5;                                   % stability limit
isnap = 10  ;                                  % plot frequency
nx    = 200  ;                                 % number of grid points in x
nz    = nx    ;                                % number of grid points in y

% Initialize pressure field 
p     = zeros(nz, nx);                             % p at time n (now)
pold  = zeros(nz, nx);                             % p at time n-1 (past)
pnew  = zeros(nz, nx);                             % p at time n+1 (present)
d2px  = zeros(nz, nx);                            % 2nd space derivative of p in x-direction
d2pz  = zeros(nz, nx);                             % 2nd space derivative of p in z-direction
c     = zeros(nz, nx);                             % velocity field
c     = c + c0;                                    % velocity model (here homogeneous)
cmax  = max(max(c)); 

% Receivers
nr    = 20;                                              % number of receievers
rec   = zeros(2, nr);              % empty receivers
seis1  =zeros(nt, nr);                           % empty seismograms (at circular receivers)

phi   = linspace (0, 2 * pi, nr+1); 
for i=1:nr
    rec(1,i) = 100 +  floor(50 * cos(phi(i)));
    rec(2,i) = 100 +  floor(50 * sin(phi(i)));
end 


isx   = 10;                                          % source location in x-direction
isz   = 20;                                          % source location in z-direction

% Grid initialization
dx    = 1. / (nx-1);                                    % space increment in x-direction
dz    = dx;                                             % space increment in z-direction

x     = (0:nx)*dx;                        % initialize space coordinates in x-direction
z     = (0:nz)*dz;                         % initialize space coordinates in z-direction
dt    = eps * dx / (cmax);                           % calculate time step from stability criterion 

% Source time function
f0    = 1. / (10. * dt);                                % dominant frequency
t     = linspace(0 * dt, nt * dt, nt);               % initialize time axis
t0    = 5. / f0;                                        % shifting of source time function

src   = zeros(nt+1);                                 
% src   = exp(-1.0* (f0.^2) .* (t-t0).^2) .* (2* sin(t-t0).^2);              % Gaussian
src   = exp(-1.0* (f0.^2) .* (t-t0).^2); 
src   = diff(src) / dt;                              % first derivative of Gaussian
src   = [src 0.0];
src   = diff(src) / dt;                              % second derivative of Gaussian
src   = [src 0.0];

%
% Plotting the source-time function
f = figure(1); clf 

plot(src); grid on; 

% label
title('Ricker Wavelet Source Time Function')
xlabel('time (s)'); ylabel('amplitude')
set( findall( f, '-property', 'FontSize' ), 'FontSize', 18 );
set( findall( f, '-property', 'FontWeight' ), 'FontWeight', 'Bold' );
set( findall( f, '-property', 'Linewidth' ),'Linewidth', 3 );
%% 
% EXTRAPOLATION SCHEME AND PLOTS
f=figure(2); clf 
% Initialize Plot
subplot(121)
imagesc(pnew); axis xy; colormap(roma.roma); hold on

% Plot Receivers
plot(rec(1,:), rec(2,:), 'k<', 'MarkerFaceColor', 'k')
colorbar
axis([0 nx 0 nz])
xlabel('x', 'fontsize', 20); ylabel('z', 'fontsize', 20)
f.Position = [64 179 1339 617];
%% propagate wavefield 
p     = zeros(nz, nx);                             % p at time n (now)
pold  = zeros(nz, nx);                             % p at time n-1 (past)
pnew  = zeros(nz, nx);                             % p at time n+1 (present)
d2px  = zeros(nz, nx);                             % 2nd space derivative of p in x-direction
d2pz  = zeros(nz, nx);                             % 2nd space derivative of p in z-direction
c     = zeros(nz, nx);                             % velocity field
c     = c + c0;                                    % velocity model (here homogeneous)
cmax  = max(max(c)); 

tic 

subplot(122)
for it =1:nt                        % 5 point operator FD scheme
   hold off 
    % Space derivative in x-direction
    for i= 3: nx - 2
        d2px(i, :) = (- 1. / 12 * p(i + 2, :) + 4. / 3  * p(i + 1, :) - 5. / 2 * p(i, :) ...
                       + 4. / 3  * p(i - 1, :) - 1. / 12 * p(i - 2, :)) / dx^2;
    end 
    
    % Space derivative in z-direction
    for j=3: nz - 2
        d2pz(:,j) = (- 1. / 12 * p(:, j + 2)+ 4. / 3  * p(:, j + 1) - 5. / 2 * p(:, j) ...
                      + 4. / 3  * p(:, j - 1) - 1. / 12 * p(:, j - 2)) / dz^2;
    end 
    
    % Time Extrapolation
    pnew = (2*p) - pold + (dt ^ 2 * cmax ^ 2 * (d2px + d2pz));
    
    % inject source term at isx, isz
    pnew(isz, isx) = pnew(isz, isx) + src(it) / (dx * dz) * (dt ^ 2); 
    
    % Remap Time Levels
    pold = p; 
    p = pnew; 
    
    % Save Seismograms 
    for i = 1: nr
        seis1(it,i) = p(rec (1,i), rec(2, i) ); 
    end 
    
    
    % see my stf that I am injecting 
    amp = max(max(pnew)); 
    cLim = [min(min(pnew))/amp  max(max(pnew))/amp];
    % imagesc(pnew, cLim); 
    imagesc(pnew./(max(max(pnew)))); 
    title(['nt = ', num2str(it)], 'fontsize', 18); 
    xlabel('x', 'fontsize', 20); ylabel('z', 'fontsize', 20)
    axis tight; hold on; 
    
    % receivers 
    plot(rec(1,:), rec(2,:), 'k<', 'MarkerFaceColor', 'k'); colormap(roma.roma); colorbar
    pause(0.0001)
    
end 
toc 
%% plot recorded seismograms 
f2=figure(4); clf 
counter = 0; 
da = 0.15; 

yyaxis left 
for kk = 1:nr
    
    p1 = plot(t, seis1(:,kk)+counter, 'b-', 'linewidth', 2);  hold on; grid on; 
    counter = counter + da; 
end 
title('recorded seismograms', 'fontsize', 20)
xlabel('time (s)', 'fontsize', 20); % ylabel('amplitude', 'fontsize', 20)
xlim([min(t) max(t)])
f.Position=[117 152 801 626]; 
f.Color='w'; 

%% back-propagate the recorded wavefield 
% re-initialize pressure and spatial derivative fields 
p     = zeros(nz, nx);                             % p at time n (now)
pold  = zeros(nz, nx);                             % p at time n-1 (past)
pnew  = zeros(nz, nx);                             % p at time n+1 (present)
d2px  = zeros(nz, nx);                             % 2nd space derivative of p in x-direction
d2pz  = zeros(nz, nx);                             % 2nd space derivative of p in z-direction
c     = zeros(nz, nx);                             % velocity field
c     = c + c0;                                    % velocity model (here homogeneous)
cmax  = max(max(c)); 

f = figure(5); clf 
for it =1:nt                        % 5 point operator FD scheme
   hold off 
    % Space derivative in x-direction
    for i= 3: nx - 2
        d2px(i, :) = (- 1. / 12 * p(i + 2, :) + 4. / 3  * p(i + 1, :) - 5. / 2 * p(i, :) ...
                       + 4. / 3  * p(i - 1, :) - 1. / 12 * p(i - 2, :)) / dx^2;
    end 
    
    % Space derivative in z-direction
    for j=3: nz - 2
        d2pz(:,j) = (- 1. / 12 * p(:, j + 2)+ 4. / 3  * p(:, j + 1) - 5. / 2 * p(:, j) ...
                      + 4. / 3  * p(:, j - 1) - 1. / 12 * p(:, j - 2)) / dz^2;
    end 
    
    % Time Extrapolation
    pnew = (2*p) - pold + (dt ^ 2 * cmax ^ 2 * (d2px + d2pz));
    
    % inject recorded wavefield as new source term
    for ii = 1:nr
        pnew(rec(2, ii), rec(1,ii)) = pnew(rec(2, ii), rec(1, ii)) + (seis1(nt-it+1,ii)*dt^2); 
    end 
    
    % reinitialize pressure fields 
    pold = p; p = pnew; 
    
    % Save Seismograms
    for i = 1:nr
        seis1(it,i) = p(rec (1,i), rec(2, i) ); 
    end 
    
    % see my stf that I am injecting 
    amp = max(max(pnew)); 
    cLim = [min(min(pnew))/amp  max(max(pnew))/amp];
    % imagesc(pnew, cLim); 
    imagesc(pnew); 
    title(['nt = ', num2str(it)], 'fontsize', 18); 
    xlabel('x', 'fontsize', 18); ylabel('z', 'fontsize', 18)
    axis tight; hold on; 
    
    % receivers 
    plot(rec(1,:), rec(2,:), 'k<', 'MarkerFaceColor', 'k'); colormap('parula'); % colorbar
    pause(0.001)
    
end 


%% 
counter = 0; 
da = 3.0E-7; 
figure(4);
yyaxis right
for kk = 1:nr
    
    p2= plot(t, seis1(:,kk)+counter, '-', 'linewidth', 1);  hold on; grid on; 
    counter = counter + da; 
end 
legend([p1, p2], 'recorded seismograms', 're-injected seismograms')

xlabel('time (s)', 'fontsize', 18); ylabel('amplitude', 'fontsize', 18)
xlim([min(t) max(t)])


