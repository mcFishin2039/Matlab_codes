% input: x - independent variable 
%        y - dependent variable 
%        xmod - places along the indepdent variable you want an estimate 
%        winsize - size of window (same units as x)

% output:  ymod - non parametric model density estimate

function ymod = nonparametric_smooth(x,y,xmod,winsize)

x = x(:); 
y = y(:); 
xmod = xmod(:);

ymod = zeros(size(xmod));

   for i = 1:length(xmod)

	dist = sqrt( (x-xmod(i)).^2 );
	ival = find( dist < winsize  );
	ival = ival(isfinite(  y(ival) ) );
	
	if isempty(ival)
		ymod(i) = NaN; 
	else
		weights = 15/16*(1 - (dist(ival)/winsize).^2).^2;
		ymod(i) = sum(weights.*y(ival))./sum(weights);

	end   


   end 

end 
