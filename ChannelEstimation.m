function [ Est_H, Err_H ] = ChannelEstimation( M, Nm, tau, U, Bar_Xi, real_H, Beta, Ptr, sigma, isTranspose )
%CHANNELESTIMATION Summary of this function goes here
%   Detailed explanation goes here


Y_tr = cell(1,M);

for m = 1:1:M
    
    Z = sigma*sqrt(1/2)*(randn(tau,Nm) + 1i*randn(tau,Nm));
%     size(real_H([(m-1)*Nm+1:m*Nm],:))
%     size(Bar_Xi)
    Y_tr{m} = sqrt(tau*Ptr)*Bar_Xi*real_H(:,[(m-1)*Nm+1:m*Nm]) + Z;
    
end

Est_H = [];
Err_H = zeros(M,U);

for m = 1:1:M
       
    Est_h = [];
    for l = 1:1:U
        Var = sigma^2;
        for ll = 1:1:U
            Var = Var + tau*Ptr*Beta(m,ll)*(Bar_Xi(:,l)'*Bar_Xi(:,ll))^2;
        end
        Est_h = [Est_h; sqrt(tau*Ptr)*Beta(m,l)/(Var)*Bar_Xi(:,l)'*Y_tr{m}];
        Err_H(m,l) = Beta(m,l)*Nm*(1 - tau*Ptr*Beta(m,l)/(Var))/sigma^2;
    end
    
    if (isTranspose)
        Est_H = [Est_H; Est_h'];
    else
        Est_H = [Est_H, Est_h];
    end
    
    
end


end

