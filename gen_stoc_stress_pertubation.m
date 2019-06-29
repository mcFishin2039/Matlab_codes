% author: Marlon Ramos 

clear; close all; clc
ax = 5;% characteristic correlation length in x, [km]
az = 5;   % characteristic correlation length in z, [km]

% distance arrays
dx = .1; dz = dx; 
x = 0.1:dx:30;    % [m]
z = 0.1:dz:30;  % [m]

kx = (2*pi)./(x);    % wave number x [rad/m]
kz = (2*pi)./(z);    % wave number z [rad/m]


[X, Z] = meshgrid(x, z);  
[KX, KZ] = meshgrid(kx, kz);  


R = sqrt( (X.^2./ax^2) + (Z.^2./az^2));       % compute distance matrix 
K = sqrt((ax^2 * KX.^2) + (az^2 * KZ.^2));    % compute wave number matrix 

%-----------------Gaussian model----------------
G = exp(-R.^2);                        % autocorrelation function
PSD_G = (ax*az/2) .* exp(-0.25 * K.^2);   % power spectral density

%-----------------Exponential model----------------
EX = exp(-R);                               % autocorrelation function
PSD_E = (ax*az) ./ ( (1 + K.^2).^(3/2) );   % power spectral density

%-----------------Von Karman model----------------
H = 1;                    % Hurst exponent
VC = besseli(H, R);       % autocorrelation function
PSD_VC = (ax*az) ./ ( (1 + K.^2).^(H) );   % power spectral density

%% plot gaussian acf and psd
f = figure(1); clf 
subplot(121)
imagesc(x, z, G); axis xy; ylabel('z'); xlabel('x'); title('gaussian acf')
colorbar; 

subplot(122)
imagesc(kx, kz, PSD_G); axis xy; ylabel('k_{z}'); xlabel('k_{x}'); title('gaussian psd')
colorbar;
f.Color = 'white'; 
set( findall( f, '-property', 'FontSize' ), 'FontSize', 18 );
%% plot exponential acf and psd
f = figure(2); clf 
subplot(121)
imagesc(x, z, EX); axis xy; ylabel('z'); xlabel('x'); title('exponential acf')
colorbar; 

subplot(122)
imagesc(kx, kz, PSD_E); axis xy; ylabel('k_{z}'); xlabel('k_{x}'); title('exponential psd')
colorbar;
f.Color = 'white'; 
set( findall( f, '-property', 'FontSize' ), 'FontSize', 18 );

%% plot von karman acf and psd
f = figure(3); clf 
subplot(121)
imagesc(x, z, VC); axis xy; ylabel('z'); xlabel('x'); title('von karman acf')
colorbar; 

subplot(122)
imagesc(kx, kz, PSD_VC); axis xy; ylabel('k_{z}'); xlabel('k_{x}'); title('von karman psd')
colorbar;
f.Color = 'white'; 
set( findall( f, '-property', 'FontSize' ), 'FontSize', 18 );




