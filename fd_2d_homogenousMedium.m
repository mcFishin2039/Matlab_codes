% script to understand 2D acoustic wave propagation in a homogeneous medium
% note that the grid is uniform.
% this script demonstrates solution accuracy with 3 and 5 point finite difference operators
% Author: Marlon
clear; close all; clc; 
%% 
%--------------- set up domain and initialize parameters --------------
nx = 500;    % grid points in x-direction
nz= nx; 

dx = 1;      % sample spacing in x
dz = dx;     % sample spacing in z
c0 = 580.0;  % acoustic wave velocity, [m/s]

isx =nx/2;  % x,z location to inject stf
isz = isx; 

irx = 330;   % receiver location in x, z
irz = irx;  

nt = 502;    % num time steps 
dt = 0.001;  % time step, (s)

eps = c0* dt/ dx;    % CFL staility criterion 

%% 
%--------------- source time function --------------
fdominant = 40.0;     % dominant source time function frequency, [Hz]
tshift = 4.0/fdominant;   % source time shift, [s]


t1 = 0*dt; t2= nt*dt; 
time = linspace(t1, t2, nt);    % time vector, [s]
src = -2 *(time - tshift).*fdominant^2 .* exp(-1.0 * (fdominant^2) * (time - tshift).^ 2);   % stf (time domain)

srf_spectrum = fft(src); 

figure(1); clf 
subplot(121); 
plot(time, src, 'linewidth', 3); grid on;  
xlim([0 max(time)])
% label
fSize = 16; 
title('Source Time Function', 'fontsize', fSize)
xlabel('time (s)', 'fontsize', fSize); ylabel('amplitude', 'fontsize', fSize)

subplot(122)
df = dt/4; 
freq = fftfreq(nt, dt/4.0); 
plot(freq, abs(srf_spectrum), 'r-', 'linewidth', 3); grid on;
xlim([0 200])
title('frequency response', 'fontsize', fSize); xlabel('frequency (Hz)', 'fontsize', fSize);

%% 
%--------------- initialize pressure field, velocity field, computation domain --------------
p=zeros(nz, nx);                 % p(i) 
pold=zeros(nz, nx);              % p(i-1)
pnew=zeros(nz, nx);              % p(i+1) 
d2px=zeros(nz, nx);              % 2nd spactial derivative of p in x
d2pz=zeros(nz, nx);              % 2nd spatial derivative of p in z

c= zeros(nz, nx);                % velocity field
c= c + c0;                       % velocity model (here homogeneous)

x = 1:nx; x = x*dx; 
z = 1:nz; z = z*dz; 

seis=zeros(1, nt);               % seismometers 
%% 
% --------------2D acoustic green's function ----------
G    = zeros(size(time)); 
r    = sqrt((x(isx) - x(irx)).^2 + (z(isz) - z(irz)).^ 2); 

for it=1:nt
    if time(it)- (abs(x(irx)-x(isx)))/c0  >= 0 
        % construct green's function
        G(it) = (1/(2 * pi * c0^2)).*(1/sqrt((time(it)^ 2) - (r.^ 2 /(c0^ 2)))); 
    end 
end 

Gc = conv(abs(G), src*dt);
figure(1); clf 
plot(Gc); 
%% 
Gc = Gc(1:nt); 
lim = max(Gc); 

%% 
figure(2); clf 
% normalize both stf and green's function
plot(src./max(src)); hold on;  
plot(Gc./lim)







