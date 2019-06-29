% wells and coppersmith regression equations 

function [M] = surface_rupture_length(a, b, SRL)
   % a, b=regression coefficients (a=5; b=1.22; reverse faults)
   % L= surface rupture length (km)
   % M=magnitude 
   
   M = a + (b*log10(SRL)); 

end 

function [N] = test(d, e)

end 