% script to pre-process .SAC binary files for cross-correlated waveforms 
% the negative and positive lag correlations are separated and stacked 
% 
% output appropriate for multiple frequency analysis 
% (i.e., Computer Programs in Seismology, $ do_mft sac_file_name.sac)
% Author: Marlon D. Ramos

clear; close all; clc; 

% establish path to sac files 
working_dir='./'; 

% read in .sac files 
fileIn='egf.IC.KMI.IC.LSA.ZZ.sym.SAC';
[data]=rdsac(fileIn); 
% [data, npts, stat, delta] = sac2mat(fileIn); 

%% 
figure(1); clf 

% cut waveforms at 0 time 
nmax=npts/2; 
data_neg=data(0:nmax); 
subplot(412)
plot(data_neg); title('neg-xcoor waveform')
data_neg=flipud(data_neg); 
subplot(413); 
plot(data_neg); title('neg-xcoor waveform  FLIPPED')
data_pos=data(nmax:end); 
subplot(414); 
plot(data_pos); title('pos-xcoor waveform')

data_stack=data_pos + data_neg; 

subplot(411); 
plot(data);  hold on; 
plot(data_stack);
legend('original waveform','stacked x-corr waveform' )
% sum negative lag to positive lag

%% 
figure(2); clf 
plot(data_neg); hold on; 
plot(data_pos, 'r-'); 
plot(data_stack, 'k-');grid on
legend('neg-xorr','pos-xcoor', 'stacked x-corr waveform' )
xlim([0 2000])
%-----write out new sac file----

% convert to .sac file 
fileOut=strcat(fileIn,'.sum'); 
mksac(file_name,seismogram,now,'DELTA',deltat,'KSTNM','TEST', 'USER1', 4e+0, 'USER2', 2.00e+2)

