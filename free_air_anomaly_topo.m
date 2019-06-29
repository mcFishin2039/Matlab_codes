% matab code to calculate the free-air gravity anomaly over any topographic
% profile 

clear; close all; clc; 
%% 
rw = 1000;   % water density, [kg/m^3]
rc = 2800;   % crustal density 
rm = 3300;   % mantle density 
G = 6.67E-11;  % gravitational constant 
g = 9.81;    % gravitational acceleration [m/s^2]
E = 6.5E10;  % Young's modulus [Pa]
nu = 0.25;   % Poisson's ratio 
s = 5000;    % mean sea level depth [m]
d = 6.0E3;   % mean crustla chickness [m]
ho= 5.0E3;   % height of gaussian seamount, [m]
sig = 2.0E4; % width of gaussian seamount

Te0 = 0.0E4;         % elastic plate thickness [m]
Do = E*Te0^3/(12*(1 - nu^2));

Te30 = 3.0E4; Do30 = E*Te30^3/(12*(1 - nu^2));

N = 2048; 
L = 4.0E6; 
dx = L/N; 
x = dx*(1:N) - (L/4); 
k = (-N/2):(N/2 -1); 
k = 2*pi*k./L;     % wave number, [radians/m]
ks = ifftshift(k); 

k4 = ks.^4; 
Tk0 = 2 * pi * G * (rc - rw ) * exp(-abs(ks*s)) .* (1.0 -exp(-abs(ks*d)))./(1.0 + (Do * k4 /(g*(rm - rc))));

Tk30 = 2 * pi * G * (rc - rw ) * exp(-abs(ks*s)) .* (1.0 - exp(-abs(ks*d)))./(1.0 + (Do30 * k4 /(g*(rm - rc))));

% topography 
topo = ho * exp(-x.*x/(2.0*sig^2));
ctopo = fft(topo); 

% corresponding gravity anomalies
grav0=real(ifft(Tk0.*ctopo));
grav30=real(ifft(Tk30.*ctopo));
%% 
f = figure(1); clf 
subplot(311)
semilogx(abs(ks), Tk0, 'b-'); hold on; 
semilogx(abs(ks), Tk30, 'r-'); grid on; 
xlabel('k (rad/m)', 'interpreter', 'latex'); 
ylabel('$$ \Delta g_{fa}(k)/H(k) $$ ', 'interpreter', 'latex'); 

subplot(312)
plot(x, topo); grid on; 
xlabel('x (m)', 'interpreter', 'latex'); 
title('topography')
ylabel('$$ h(x) $$ ', 'interpreter', 'latex'); 

subplot(313)
plot(x, grav0, 'k-'); hold on;  
plot(x, grav30, 'k--'); grid on;  
xlabel('x (m)', 'interpreter', 'latex'); 
ylabel('$$ g_{fa}(x) $$ ', 'interpreter', 'latex'); 

set( findall( f, '-property', 'FontSize' ), 'FontSize', 18 );
set( findall( f, '-property', 'LineWidth' ), 'LineWidth', 1.5 );
