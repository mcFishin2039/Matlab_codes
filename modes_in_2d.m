clear
close all
clc

% addpath to cool new scientific colormaps 
roma=load('~/ScientificColourMaps5/roma/roma.mat'); 

%% 
% This script computes the normal modes of a plate

Ly = 1.0;
Lx = 1.0;

nx = 101;
ny = 100;
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);

% make a vector of modes here in case you want to sum over modes
m = 1 : 20;
n = 1 : 20;

% allocate the solution matrix
modes = zeros(ny,nx);

% compute the term due to first mode
m_idx = 4;
n_idx = 6;
for ii = 1 : ny
    for jj = 1 : nx
        modes(ii,jj) = sin( m(m_idx) * pi * x(jj)/Lx ) * sin( n(n_idx) * pi * y(ii)/Ly );
    end
end

% compute the term due to second mode
m_idx = 1;
n_idx = 2;
for ii = 1 : ny
    for jj = 1 : nx
        modes(ii,jj) = modes(ii,jj) + sin( m(m_idx) * pi * x(jj)/Lx ) * sin( n(n_idx) * pi * y(ii)/Ly );
    end
end
% If you want a single mode, just comment the loop above

%% 
% plot the mode combination
h = figure(1);
imagesc(x,y,modes); 
colormap(roma.roma); 
c = colorbar;
title('modes 53 + 21');
xlabel('X-coord'); ylabel('Y-coord');
axis('square');
set( findall( h, '-property', 'FontSize' ), 'FontSize', 18 );
set( findall( h, '-property', 'FontWeight' ), 'FontWeight', 'Bold' );
set( findall( h, '-property', 'Linewidth' ), 'Linewidth', 2 );
% print(h,'-dpng', '-r500' , 'mode12_highres.png')
%% 
h = figure('color','w');
contourf(x,y,modes); c = colorbar;
title('modes 53 + 21');
xlabel('X-coord'); ylabel('Y-coord');
axis('square');
print(h,'-dpng' , 'mode12-lowres.png')

%% how about a plate that moves through time

nt = 100;
t = linspace(0,1,nt);
vel = 1; % [m/s] phase velocity

modes = zeros( ny, nx );

m_idx =1;
n_idx = 2;

% compute the eigen frequency
omega = (pi * vel) * sqrt( ( m(m_idx) / Lx )^2 + ( n(n_idx) / Ly )^2 );

h = figure('Color','w');
for kk = 1 : nt
    for ii = 1 : ny
        for jj = 1 : nx
            modes(ii,jj) = cos( omega * t(kk) ) * sin( m(m_idx) * pi * x(jj)/Lx ) * sin( n(n_idx) * pi * y(ii)/Ly );
        end
    end
    
    clf
    surf(x,y,modes); shading('interp'); c = colorbar; colormap(roma.roma)
    title(['m=' num2str(m(m_idx)) ', n=' num2str(n(n_idx)) ': t=' num2str(t(kk),'%0.2f') ' [s]']);
    xlabel('X-coord'); ylabel('Y-coord'); zlabel('Amplitude (a.u)');
    axis([0 Lx 0 Ly -1 1]);
    caxis([-1 1]);
    pause(0.001);

end

