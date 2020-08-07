function [ StartingPoint ] = GetInitialization( M, K, L, Nm, sigma, H_d, H_u, G_SI, G_CCI, nu, Pd, Pu, eta, C_SE, C_EE, Pa, Ps, PAP, Pl, SysModel, APConn, Delta, ErrorEst )
%GETINITIALIZATION Summary of this function goes here
%   Detailed explanation goes here

if (nargin<nargin('GetInitialization')-1)
    Err_H_DL = zeros(M,K);
    Err_H_UL = zeros(M,L);
    Err_GCCI = zeros(K,L);
    Delta = 0.99;
elseif (nargin<nargin('GetInitialization'))
    Err_H_DL = zeros(M,K);
    Err_H_UL = zeros(M,L);
    Err_GCCI = zeros(K,L);
else
    Err_H_DL = ErrorEst{1};
    Err_H_UL = ErrorEst{2};
    Err_GCCI = ErrorEst{3};
end

Method = SysModel(1)
DLULSchemes = SysModel(2)

APs_DLs = APConn{1};
APs_ULs = APConn{2};

epspow = 10^(-3)*(PAP/nu(1)+Pd+Pu)

temp_H_d = H_d;
if (Method==2 && DLULSchemes==1)
    GG = G_SI'*G_SI;
    [U, S, V] = svd(GG);
    s = diag(S);
    for NoEig = 1:1:length(s)-1
        if (sum(s(1:NoEig))/sum(s)<Delta)
            continue
        else
            break
        end
    end
    Infpercent = sum(s(1:NoEig))/sum(s)
    NoEig
    U1 = U(:,1:NoEig);
%     U1 = U(:,1:end-32);
    
    
    P = eye(size(U1,1)) - U1*U1';
    temp_H_d = H_d*P;
end

[Q, R] = qr(temp_H_d',0);
T = R';
Q = Q';


% TranfMat = [];
% for m = 1:1:M
%     TranfMat = [TranfMat; eye(K)];
% end
 
% omega = 1*rand(1,K*M);
omega = 1*rand(1,K);%0.01*PAP/K

if (Method>=5)
    omega = 100*rand(1,K);
end

T_tilde = zeros(K,K);
for i = 1:1:K
    for j = 1:1:K
        if (i<j)
            T_tilde(i,j) = 0;
        elseif (i==j)
            T_tilde(i,j) = 1;
        else
            T_tilde(i,j) = -1/T(i,i)*(T(i,[j:(i-1)])*T_tilde([j:(i-1)],j));
        end
    end
end

if (Method==4 || Method==7)
    
    QT = zeros(M*Nm, K);
%     size(QT)
    
    for m=1:1:M
        idx_k = find(APs_DLs(m,:)>0);
        if (length(idx_k)==1)
            QT([(m-1)*Nm+1:m*Nm], idx_k) = H_d(idx_k,[(m-1)*Nm+1:m*Nm])';
        elseif (length(idx_k)>1)
%             size(QT([(m-1)*Nm+1:m*Nm], idx_k))
%             size(H_d(idx_k,[(m-1)*Nm+1:m*Nm]))
            QT([(m-1)*Nm+1:m*Nm], idx_k) = H_d(idx_k,[(m-1)*Nm+1:m*Nm])'*inv(H_d(idx_k,[(m-1)*Nm+1:m*Nm])*H_d(idx_k,[(m-1)*Nm+1:m*Nm])');
        end
    end
    
    for k = 1:1:K
        QT(:,k) = QT(:,k)/norm(QT(:,k));
    end
    
else

    if (DLULSchemes==1)

        QT = Q'*T_tilde; 

    elseif (DLULSchemes==2)

        QT = H_d'*inv(H_d*(H_d)');

    else

        QT = H_d';

    end

    if (DLULSchemes~=1)
        for k = 1:1:K
            QT(:,k) = QT(:,k)/norm(QT(:,k));
        end
    end
end

% QT_Blk = [];
% for m = 1:1:M
%     QT_Blk = blkdiag(QT_Blk,QT([(m-1)*Nm+1:m*Nm],:));
% end


W = QT*diag(omega.^(0.5));

DLMask = QT;

if (Method==2 && DLULSchemes==1)
    W = P*W;
    DLMask = P*QT;
end

% W = QT*diag(omega.^(0.5));

p = 1*rand(1,L);%0.01*Pl

if (Method>=5)
    p = 10*rand(1,L);
end

mu = ones(1,M);%rand(1, M);%%0.5*

psi_d = zeros(1,K);

for k = 1:1:K
    
    sum_hw_k = 0;
    for kp = 1:1:K
        if (kp==k)
            continue
        end
%         size(W)
        sum_hw_k = sum_hw_k + norm(H_d(k,:)*W(:,kp))^2;
    end
    
    
%     psi_d(k) = 1/(norm(H_d(k,:)*QT(:,k)))^2*(sum_hw_k+norm(G_CCI(k,:).* (p.^(0.5)) )^2 + sigma^2 );
    psi_d(k) = 1/(norm(H_d(k,:)*DLMask(:,k)))^2*(sum_hw_k+norm(G_CCI(k,:).* (p.^(0.5)) )^2 + sigma^2 + ...
                                                    norm(kron((Err_H_DL(:,k)'/Nm).^(0.5),ones(1,Nm))*W)^2 + ...
                                                    norm(Err_GCCI(k,:).* (p.^(0.5)))^2);
    
end


lambda_d = omega./psi_d;

psi_u = zeros(1,L);
ULMask = [];

if (Method==4 || Method==7)
    
    A_SIC = zeros(L, M*Nm);
    
    for m=1:1:M
        idx_l = find(APs_ULs(m,:)>0);
        if (length(idx_l)==1)
            A_SIC(idx_l, [(m-1)*Nm+1:m*Nm]) = H_u([(m-1)*Nm+1:m*Nm], idx_l)';
        elseif (length(idx_l)>1)
            A_SIC(idx_l, [(m-1)*Nm+1:m*Nm]) = inv(H_u([(m-1)*Nm+1:m*Nm], idx_l)'*H_u([(m-1)*Nm+1:m*Nm], idx_l))*H_u([(m-1)*Nm+1:m*Nm], idx_l)';
        end
    end
    
    for l = 1:1:L
        psi_u(l) = 1/(A_SIC(l,:)*H_u(:,l))^2*( norm(A_SIC(l,:)*H_u(:,[1:(l-1) (l+1):L])*diag(p([1:(l-1) (l+1):L]).^(0.5)))^2 + ...
                                                   norm(A_SIC(l,:)*G_SI*W)^2 + sigma^2*norm(A_SIC(l,:))^2 + ...
                                                    norm(A_SIC(l,:)*kron((Err_H_UL/Nm).^(0.5),ones(Nm,1))*diag(p.^(0.5)))^2);
        ULMask = [ULMask; A_SIC(l,:)./norm(A_SIC(l,:))];
    
    end
else

if (DLULSchemes==2 || DLULSchemes==3)
    if (DLULSchemes==2)
        A_SIC = inv(H_u'*H_u)*H_u';
    else
        A_SIC = H_u';
    end
    for l = 1:1:L
        psi_u(l) = 1/(A_SIC(l,:)*H_u(:,l))^2*( norm(A_SIC(l,:)*H_u(:,[1:(l-1) (l+1):L])*diag(p([1:(l-1) (l+1):L]).^(0.5)))^2 + ...
                                                   norm(A_SIC(l,:)*G_SI*W)^2 + sigma^2*norm(A_SIC(l,:))^2 + ...
                                                    norm(Err_H_UL*(p'.^(0.5)))^2);
        ULMask = [ULMask; A_SIC(l,:)./norm(A_SIC(l,:))];
    
    end
else
    for l = 1:1:L
        A_SIC = inv(H_u(:,[l:L])'*H_u(:,[l:L]))*H_u(:,[l:L])';
        if (l==L)
            psi_u(l) = 1/(A_SIC(1,:)*H_u(:,l))^2*( norm(A_SIC(1,:)*G_SI*W)^2 + sigma^2*norm(A_SIC(1,:))^2 + ...
                                                    norm(Err_H_UL*(p'.^(0.5)))^2 );
        else
            psi_u(l) = 1/(A_SIC(1,:)*H_u(:,l))^2*( norm(A_SIC(1,:)*H_u(:,[(l+1):L])*diag(p([(l+1):L]).^(0.5)))^2 + ...
                                                   norm(A_SIC(1,:)*G_SI*W)^2 + sigma^2*norm(A_SIC(1,:))^2 + ...
                                                    norm(Err_H_UL*(p'.^(0.5)))^2);
        end
        ULMask = [ULMask; A_SIC(1,:)./norm(A_SIC(1,:))];
    %     A_ZFSIC = [A_ZFSIC; A_SIC(1,:)];
    end 
end
end

lambda_u = p./psi_u;

xi = zeros(1,M);
sum_chi_xi = 0;

for m = 1:1:M
    sum_norm_w = 0;
    for k = 1:1:K
        w_km = W([((m-1)*Nm+1):(m*Nm)],k);
        sum_norm_w = sum_norm_w + norm(w_km)^2*(1/nu(1)+Pd/(norm(w_km)^2+epspow));
    end
    xi(m) = sum_norm_w; %+Pa+Pu
    sum_chi_xi = sum_chi_xi + (xi(m)/(2*mu(m))*mu(m)^2+mu(m)/(2*xi(m))*xi(m)^2);
    
%     DLEffPow_current(m) = sqrt(sum_norm_w);
%     ULEffPow_current(m) = sqrt(sum(p.*real(diag(ULMask(:,[(m-1)*Nm+1:m*Nm])*H_u([(m-1)*Nm+1:m*Nm],:))'.^2)));
end

chi = sum_chi_xi + sum(p) + sum(mu*Pu) + sum(mu*Pa) + sum((1-mu)*Ps);

phi = 1/(eta*C_SE+(1-eta)*C_EE/chi);

phi_bar = 1/(eta*C_SE*chi + (1-eta)*C_EE);

% phi = phi_bar*chi;

DLSel_current = ones(K,M);
ULSel_current = ones(L,M);


% StartingPoint = {DLMask, ULMask, omega,  p, mu, lambda_d, lambda_u, psi_d, psi_u, xi, {phi, phi_bar}, chi, {DLEffPow_current, ULEffPow_current} };
StartingPoint = {DLMask, ULMask, omega,  p, mu, lambda_d, lambda_u, psi_d, psi_u, xi, {phi, phi_bar}, chi, {DLSel_current, ULSel_current} };

end





