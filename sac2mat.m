%GET_SAC
%
%	get_sac
%
%	MATLAB M-file to read data directly from a SAC file.
%	Once the file is selected, it brings the following
%	data in the workspace:
%		npts:	number of data points
%		delta:	sampling interval
%		date:	year month and day of the begining of file
%		hour:		time of
%		minu:		begining
%		sec:		of file
%		stat:   station name
%		<station name>:	trace data
%
%	Interactive form.  See also fget_sac for programable form.
%
%	

%{
Guy Tytgat, 17 feb 95
[file,path] = uigetfile('*','Select a SAC file',300,300);

if isempty(file)
  disp('Error selecting file')
  clear file path
  return
elseif file == 0
  disp('No file selected')
  clear file path
  return
end
%%filename = sprintf('%s%s',path,file);
%}
function [data, npts, stat, delta, station, KCMPNM,KNETWK, KUSER0, KUSER1] = sac2mat(file)


filename = sprintf('%s%s',file);

[fid,message] = fopen(filename);

if fid ~= -1
  [fheader,count] = fread(fid,70,'float');
  [nheader,count] = fread(fid,35,'long');
  [lheader,count] = fread(fid,5,'long');
  [kheader,count] = fread(fid,[8,24],'char');
  station = setstr(kheader(:,1)'); % station name 
  KCMPNM=setstr(kheader(:,21)');   % station component
  KNETWK=setstr(kheader(:,22)');   % network
  KUSER0=setstr(kheader(:,18)'); 
  KUSER1=setstr(kheader(:,19)'); 
%%  station = station(f_nblank(station):length(station));	% remove leading blank
%%  station = station(1:f_blank(station)-1);	% remove all from next blank on
  stat = station;
%%  if (station(1) < 65 | station(1) > 90) & (station(1) < 97 | ...
%%     station(1) > 122)
%%    station = ['s',station];
%%  end

depmin =fheader(2); depmax = fheader(3);  odelta = fheader(5); 
b = fheader(6); e = fheader(7);  parriv = fheader(9); 
sarriv = fheader(11); pparriv = fheader(12);  ssarriv = fheader(13); 
scsarriv = fheader(14); pcparriv = fheader(15); 

  npts = nheader(10);
  delta = fheader(1);
  date = nheader(1:2);
  hour = nheader(3);
  minu = nheader(4);
  sec = nheader(5) + (nheader(6)/1000);

  data = zeros(npts,1);
  [data,count] = fread(fid,npts,'float');

%%  eval([station ' = data;'])

  close_stat = fclose(fid);
%  clear fheader nheader lheader kheader data station
%  clear file path filename fid message count close_stat
else
  disp(message)
  clear file path filename fid message
end

end 


