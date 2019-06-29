% plots a filtered matrix of data
% filtering method is Butterworth with np poles
% displays the data from bottom to the top 
% Input: DATA matrix, sample rate, lowest frequency (Hz), highest frequency (Hz),
% number of poles for the filter, option to filter 1 or not 0
% Output DATAf = filtered data matrix
%% IMPORTANT: the matrix is plotted as column 3; column 2 and column 1 from top to bottom : ie Z, N, E

function [DATAf] = plot_filtered_matrix(DATA, samprate, flo, fhi, np, filter_option,plot_option)

% Butterworth filter parameters
[pd,qd]=size(DATA);
if length(samprate(:,1))==1
    'same sample rate'
    samprate =ones(length(DATA(1,:)),1)*samprate; 
end 
DATA=detrend(DATA,'constant'); 
DATAf=[]; 

% 0 phase butterworth filter
for kd=1:qd
    % [b,a]=butter(np, [flo/samprate(kd)*2 fhi/samprate(kd)*2]);
    [b,a]=butter(np, [flo/samprate(kd) fhi/samprate(kd)]);
    if filter_option == 1
        DATAf(:,kd)=filtfilt(b,a,DATA(:,kd));   % zero phase filtered waveforms
    else
        DATAf(:,kd)=DATA(:,kd);
    end
end

if plot_option==1
    ld=length(DATAf(1,:));
    MA= max(abs(DATAf(1,1:ld-1)));
    lm = length(MA);
    %ld=length(DATAf(:,1));
    plot((1:ld-1)/samprate(1),DATAf(1:ld-1,qd)); hold on; % zoom on
    % text(1,0,[num2str(qd)])
    for km = 1: qd-1
        plot( (1:ld-1)/samprate(km),DATAf(1,qd-km)+sum(MA(qd-km:qd))); % zoom on;
        text(1,sum(MA(qd-km:qd)),[num2str(qd-km)])
    end
    xlabel('time (sec) from the origin time')
    title('bp filtered surface waveform  ')
end % if plot_option

end 