% script to demonstrate some key multiple filter analysis and 
% frequency time analysis concepts 
% MFA - Dziewonski et al., (1969) BSSA
% FTAN - Levshin et al., (1972)
clear; close all; clc; 
T=0.001:0.01:2.5;    % period, (s)
f=1./T;           % frequency, (1/s)
omega=2*pi*f;    % angular frequency, (rad/s)

alpha=200;      % smoothing parameter

% ---window function----
omega_n=50;     % some center frequency 
Hn1=exp(-alpha*((omega-omega_n)./omega_n).^2);  % frequency-domain

Hn2=(sqrt(pi)*omega_n/(2*alpha)) * exp(-(omega_n^2*T.^2)/(4*alpha)) .* cos(omega_n*T);

%--- plot ---
figHandle=figure(1); clf
subplot(211)
plot(omega, Hn1, 'b-'); grid on; 
xlabel('\omega'); title('H_{n}(\omega)')
subplot(212)
plot(T, Hn2, 'k-'); grid on; 
xlabel('t'); title('h_{n}(t)')
set( findall( figHandle, '-property', 'FontSize' ), 'FontSize', 14 );




