% example with partial differential equaiton solver 
% 
clear; close all; clc; 

% wave equation coeffiecients (homogeneous wave equation)
c = 1;
a = 0;
f = 0;     % force term = 0
m = 1;

%% 
numberOfPDE = 1;
model = createpde(numberOfPDE);
geometryFromEdges(model,@squareg);
pdegplot(model,'EdgeLabels','on');
ylim([-1.1 1.1])
axis equal
title('Geometry With Edge Labels Displayed');
xlabel('x'); ylabel('y')

%% 
specifyCoefficients(model,'m',m,'d',0,'c',c,'a',a,'f',f);

%% Boundary Conditions
applyBoundaryCondition(model,'dirichlet','Edge',[2,4],'u',0);
applyBoundaryCondition(model,'neumann','Edge',([1 3]),'g',0);

%% 
generateMesh(model);
figure
pdemesh(model);
ylim([-1.1 1.1]);
axis equal
xlabel x
ylabel y

