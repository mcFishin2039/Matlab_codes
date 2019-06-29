% script to make 1-D static finite element wave propagation 
% based off computational seismology text 
% Author: Marlon D. Ramos 

clear; close all; clc; 

% ----initialization ----
nx=50;            % number of points 
u=zeros(1,nx);    % displacement 
f=zeros(1,nx);    % source vector 
mu=1.0;           % shear modulus, [Pa]
x=linspace(0, 1, nx);
h=x(2) - x(1);    % element size
K=zeros(nx, nx);  % stiffness matrix, K

% populate K matrix 
for ii=1:nx
    for jj=1:nx
        if ii==jj
            K(ii,jj)=2*mu/h;
        elseif ii==jj+1
            K(ii,jj)=-mu/h;
        elseif jj==ii+1
            K(ii,jj)=-mu/h;
        else 
            K(ii,jj)=0; 
        end 
    end 
end 

% inject source term 
f(round(3*nx/4))=1.0;

% boundary conditions 
u(1) = 0.15; 
f(2)=u(1)/h;

% finite element solution 
u(2:end)=inv(K(2:end, 2:end)) * f(2:end)'; 

%----- plot------ 
figHandle=figure('color', 'w'); clf 
plot(x, u, 'r-', 'linewidth', 2); grid on; 
fSize=17;
title('1D static elasticity solution')
legend('finite elements')
ylabel('$$ u(x) $$ ', 'interpreter', 'latex'); 
xlabel('$$ x $$ ', 'interpreter', 'latex')
axis([0, 1, 0.04, .28])
set( findall( figHandle, '-property', 'FontSize' ), 'FontSize', 18 );



