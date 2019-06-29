% elastic dislocation models for dip-slip faults 
% Author: Marlon 
clear; close all; clc
% ----- Freund & Barnett (1976) ------
delta_u=1.0;       % slip (m)
theta=20;          % dip angle (degrees)
D=15;              % max fault depth (m)
x2=0:0.1:100;   % along-surface distance (m)
xD=D/tand(theta); 
xP=D/(cosd(theta)*sind(theta)); 

% vertical displacement
u3=delta_u*sind(theta)/pi * ( atand(D./(x2-D)) - (x2*D)./((x2-xD).^2+D^2) ...
                - pi/2*(1 - sign(x2)));
% horizontal displacement 
u2=delta_u*cosd(theta)/pi * ( atand((x2-D)./D) - ((x2-xP)*D)./((x2-xD).^2+D^2) ...
                - pi/2*(sign(x2)));
% plot               
figHandle = figure('color', 'w'); 
subplot(211)
plot(x2./D, u3, 'linewidth', 2); grid on; 
xlabel('x_{2}'); ylabel('displacement'); title('u_{3}')
subplot(212)
plot(x2, u2, 'b-', 'linewidth', 2); grid on; 
xlabel('x_{2}'); ylabel('displacement'); title('u_{2}')
set( findall( figHandle, '-property', 'FontSize' ), 'FontSize', 14 );

%------ Rani & Sing (199X) ------






