function r=icemodel_rmse(v,z,us,P)
% calculates rmse for ice flow model
% INPUT: v=observed velocity
%               z = observed depth
%              P = [A,n] flow law parameters
% OUTPUT: r = RMSE for model

    rho=917;
    g=9.8;
    theta=10*pi/180;
    T1=rho*g*sin(theta);
    vm=us-P(1).*T1.^P(2).*z.^(P(2)+1);
    r=sqrt(mean((v-vm).^2));

end 