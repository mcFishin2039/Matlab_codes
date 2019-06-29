% script to read and write a sac data file 
% Author: Marlon Ramos 

data=load('eq_marlon'); 

%% 
seismogram=data.x; 
samp_rate=data.sr; 

N=length(seismogram); 
t=(1:N)/samp_rate;    % seconds
deltat=t(2)-t(1);     % sample spacing, (s)
% plot(t, seismogram) % q/c

% convert to .sac file 
file_name='TA.marlon.BHZ.sac';   
mksac(file_name,seismogram,now,'DELTA',deltat,'KSTNM','TEST', 'USER1', 4e+0, 'USER2', 2.00e+2)

