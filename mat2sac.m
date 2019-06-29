% Actions: write sac files from a matrix in memory
% July 25 2008
% INPUT file_name,DATA, deltaf, fheader, nheader, lheader,kheader

% saves data in a file filename.sac 

function [count1, count2, count3, count4,count5, count6]=mat2sac(file_name,data, deltaf, fheader, nheader, lheader,kheader)
%%%%%%%%%%% reads the names of the files from dummym
   'I am writing'
		wfilename = [file_name '.sac']
        [fid,message] = fopen(wfilename,'a');
        count1 = fwrite(fid,fheader,'float');
        [count2] = fwrite(fid,nheader,'long');
        [count3] = fwrite(fid,lheader,'long');
        [count4] = fwrite(fid,kheader,'char'); 
        [count5] = fwrite(fid,data,'float');
        [count6] = fwrite(fid,deltaf, 'float');
 	    close_stat = fclose(fid);
       
end 
  
    


