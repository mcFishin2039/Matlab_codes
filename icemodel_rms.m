% function to compute the rms error between successive model estimations of
% depth and the true depth
% 
% input: u - ix component of ice flow in the z direction
%        z - depth of measurements
%        a - A parameter 
%        n - n parameters
%        us - surface velocity u(z=0)
% 
% output: root mean square error between observed and predicted
% 

function rms = icemodel_rms(u,z,a,n,us)

rho = 917; % density
g = 9.81;  % gravitational acceleration
theta = 10*pi/180; % slope angle [degrees] --> [radians]
T1 = rho*g*sin(theta); % constant value

um = us - a.*(T1.^n).*z.^(n+1); % theoretical model
rms = sqrt(sum(u - um).^2 ); % root mean square error


end 
