function [ OptimalValue, DLRate_PerUser_next, ULRate_PerUser_next, Solution_next, Status, PowScale ] = Get_optSolutionPerIteration( AllParameters, Solution_current, Iter, Init, prePowScale )
%GET_OPTSOLUTIONPERITERATION Summary of this function goes here
%   Detailed explanation goes here

penpar = 10^2;
Iter = Iter + 1;

pThres = 10^(-3);

M = AllParameters{1};
K = AllParameters{2};
L = AllParameters{3};
Nm = AllParameters{4};
sigma = AllParameters{5};
H_d = AllParameters{6};
H_u = AllParameters{7};
G_SI = AllParameters{8};
G_CCI = AllParameters{9};
nu = AllParameters{10};
P_AP = AllParameters{11};
P_UL = AllParameters{12};
Pd = AllParameters{13};
Pu = AllParameters{14};
eta = AllParameters{15};
C_SE = AllParameters{16};
C_EE = AllParameters{17};
if (length(AllParameters{18})==1)
    DLRate_Threshold = AllParameters{18};
    ULRate_Threshold = AllParameters{18};
else
    DLRate_Threshold = AllParameters{18}(1);
    ULRate_Threshold = AllParameters{18}(2);
end
Pa = AllParameters{19}(1);
PcirAP = AllParameters{19}(2);
PcirUE = AllParameters{19}(3);
Ps = AllParameters{20};
Pbh = AllParameters{21};
DLULSchemes = AllParameters{22};

if (length(AllParameters)<23)
    Err_H_DL = zeros(M,K);
    Err_H_UL = zeros(M,L);
    Err_GCCI = zeros(K,L);
else
    Err_H_DL = AllParameters{23}{1,1};
    Err_H_UL = AllParameters{23}{1,2};
    Err_GCCI = AllParameters{23}{1,3};
end

if (norm(G_SI)<10^(-20) && norm(G_CCI)<10^(-20))
    HD = 1;
else
    HD = 0;
end

% P_AP
Init;

epspow = 10^(-2)*(P_AP/nu(1)+Pd); % *3^(-Iter)+Pu+Pa-Ps
epspowDL = epspow; %10^(-3)*(P_AP)
epspowUL = 10^(-2)*P_UL(1);

if (prePowScale.keepSolution)
   epspow = epspow; %* 2^(prePowScale.epsScale);
end

epspow;

DLMask = Solution_current{1};
ULMask = Solution_current{2};
omega_current = Solution_current{3};
p_current = Solution_current{4};
mu_current = Solution_current{5};
lambda_d_current = Solution_current{6};
lambda_u_current = Solution_current{7};
psi_d_current = Solution_current{8};
psi_u_current = Solution_current{9};
xi_current = Solution_current{10};
phi_current = Solution_current{11}{1};
phi_bar_current = Solution_current{11}{2};
chi_current = Solution_current{12};
DLSel_current = Solution_current{13}{1};
ULSel_current = Solution_current{13}{2};


omega_next = sdpvar(1, K, 'full','real');
p_next = sdpvar(1, L, 'full','real'); 
% mu_next = ones(1,M);%sdpvar(1, M, 'full','real');%
lambda_d_next = sdpvar(1, K, 'full','real');
lambda_u_next = sdpvar(1, L, 'full','real');
psi_d_next = sdpvar(1, K, 'full','real');
psi_u_next = sdpvar(1, L, 'full','real');
xi_next = sdpvar(1, M, 'full','real');
phi_next = sdpvar(1, 1, 'full','real');
chi_next = sdpvar(1, 1, 'full','real');

W_current = DLMask*diag(omega_current.^(0.5));
DLPowPercent = zeros(K,M);
    
A = H_u*diag(p_current.^(0.5));
ULPowPercent = zeros(L,M);

for m = 1:1:M
    
    for k = 1:1:K
        Wmk = H_d(k,(m-1)*Nm+1:m*Nm)*W_current((m-1)*Nm+1:m*Nm,k);
        DLPowPercent(k,m) = (norm(Wmk)^2)/(norm(H_d(k,:)*W_current(:,k))^2);
    end

    for l = 1:1:L
        Aml = ULMask(l,[(m-1)*Nm+1:m*Nm])*A([(m-1)*Nm+1:m*Nm],l);
        ULPowPercent(l,m) = (norm(Aml)^2)/(norm(ULMask(l,:)*A(:,l))^2);%trace(Am'*Am);
    end

end

maxMu = max([DLPowPercent' ULPowPercent'],[],2)';

if (~(HD && P_UL(1)<10^(-3)))
    DLMask = diag(kron(mu_current,ones(1,Nm)))*DLMask;
end


if (~(HD && P_AP<10^(-3)))
    ULMask = ULMask*diag(kron(mu_current,ones(1,Nm)));
end


% W_next = DLMask*diag(omega_next);
% mu_estimate = zeros(1,M);
% for m = 1:1:M
%     mu_estimate(m) = norm(vec(W_current([(m-1)*Nm+1:m*Nm],:)))^2/(norm(vec(W_current([(m-1)*Nm+1:m*Nm],:)))^2 + epspow);
% %     mu_estimate(m) = xi_current(m)/(xi_current(m)+epspow);
% end

cons = [];

%% Cons. (14c)

for l = 1:1:L
    cons = [cons, p_next(l)>=0];
    cons = [cons, p_next(l)<=P_UL(l)];
end


%% Cons. (22) (34e)

% pen_mu = sdpvar(1, M, 'full','real');
for m = 1:1:M
    
    if (HD && P_AP<10^(-3))
        continue
    end
    
    Sm = zeros(1,M*Nm);
    Sm( ((m-1)*Nm+1):(m*Nm) ) = ones(1,Nm);
    Sm = diag(Sm);
%     cons = [cons, trace(DLMask'*Sm*DLMask*diag(omega_next))<=mu_next(m)*P_AP];
%     cons = [cons, sum(real(diag(DLMask'*Sm*DLMask))'.*omega_next)<=mu_next(m)*P_AP];

    
    cons = [cons, sum(real(diag(DLMask'*Sm*DLMask))'.*omega_next)<=(mu_current(m)+10^(-5))*P_AP]; % mu_current(m)*
%     cons = [cons, sum(real(diag(DLMask'*Sm*DLMask))'.*omega_next)<=xi_current(m)/(xi_current(m)+epspow)*P_AP];
%     cons = [cons, sum(real(diag(DLMask'*Sm*DLMask))'.*omega_next)<=mu_estimate(m)*P_AP];

%     cons = [cons, mu_next(m)>=0];
%     cons = [cons, mu_next(m)<=1];
    
%     cons = [cons, pen_mu(m)<=0];
%     cons = [cons, 2*mu_current(m)*mu_next(m)-mu_current(m)^2-mu_next(m)>=pen_mu(m)];
%     cons = [cons, mu_next(m)>=mu_next(m)^2+pen_mu(m)];
end

%% Cons. (40a) (41a)

omega_bar = sdpvar(1, K, 'full','real');
for k = 1:1:K
    
    if (HD && P_AP<10^(-3))
        cons = [cons, lambda_d_next(k) >=0];
        cons = [cons, lambda_d_next(k) <=10^(-4)];
        cons = [cons, omega_next(k) >=0];
        continue
    end
    
    cons = [cons, psi_d_next(k) >= 1/(norm(H_d(k,:)*DLMask(:,k))^2)* ...
                                    ( sum( abs(H_d(k,:)*DLMask(:,[1:(k-1) (k+1):K])).^2.*omega_next([1:(k-1) (k+1):K])) + ...
                                      sum(abs(G_CCI(k,:)).^2.*p_next) + sigma^2 + ...
                                      sum( abs(kron((Err_H_DL(:,k)'/Nm).^(0.5),ones(1,Nm))*DLMask).^2.*omega_next) + ... 
                                      sum(abs(Err_GCCI(k,:)).*p_next) )];
                                  
%     cons = [cons, psi_d_next(k) >=0];
    cons = [cons, omega_next(k) >=0];
%     cons = [cons, rcone(omega_bar(k), omega_next(k), 1/2)];
    cons = [cons, cone([omega_bar(k), 0.5*(omega_next(k)-1)], 0.5*(omega_next(k)+1) )];
    cons = [cons, 2*sqrt(omega_current(k))/psi_d_current(k)*omega_bar(k)-omega_current(k)/(psi_d_current(k)^2)*psi_d_next(k)>=lambda_d_next(k)];
    cons = [cons, omega_bar(k)>=0];
    cons = [cons, lambda_d_next(k) >=0];
    
    
    cons = [cons, 2*sqrt(omega_current(k))/psi_d_current(k)*omega_bar(k)-omega_current(k)/(psi_d_current(k)^2)*psi_d_next(k)>=0];
end


%% Cons. (40b) (41b)

p_bar = sdpvar(1, L, 'full','real');
ULMask_temp = ULMask;
for l = 1:1:L
    
    if (HD && P_UL(l)<10^(-3))
        cons = [cons, lambda_u_next(l) >=0];
        cons = [cons, lambda_u_next(l) <=10^(-4)];
        cons = [cons, psi_u_next(l)>=0];
        continue
    end
    
%     if (~(HD && P_AP<10^(-3)))
% %         ULMask = ULMask*diag(kron(xi_current./(xi_current+epspow).*mu_current,ones(1,Nm)));
%         ULMask = ULMask*diag(kron(mu_current,ones(1,Nm)));
%     end
    
    if (DLULSchemes==1)
    
        if (l==L)
            cons = [cons, psi_u_next(l) >= 1/(norm(ULMask(l,:)*H_u(:,l))^2)* ...
                                        ( sum( abs(ULMask(l,:)*G_SI*DLMask).^2.*omega_next) + sigma^2*norm(ULMask(l,:))^2 + ...
                                          sum( abs(ULMask(l,:)*kron((Err_H_UL/Nm).^(0.5),ones(Nm,1))).^2.*p_next) )];
        else

%             GSIDLMask = norm(G_SI*DLMask);
%             UMaskHu_SI = abs(ULMask(l,:)*G_SI*DLMask).^2;
%             UMaskHu_noise = sigma^2*norm(ULMask(l,:))^2;
            cons = [cons, psi_u_next(l) >= 1/(norm(ULMask(l,:)*H_u(:,l))^2)* ...
                                        ( sum( abs( ULMask(l,:)*H_u(:,(l+1):L) ).^2.*p_next((l+1):L) ) + ...
                                          sum( abs(ULMask(l,:)*G_SI*DLMask).^2.*omega_next) + sigma^2*norm(ULMask(l,:))^2 + ...
                                          sum( abs(ULMask(l,:)*kron((Err_H_UL/Nm).^(0.5),ones(Nm,1))).^2.*p_next) )];
        end
    else
        cons = [cons, psi_u_next(l) >= 1/(norm(ULMask(l,:)*H_u(:,l))^2)* ...
                                        ( sum( abs( ULMask(l,:)*H_u(:,[1:(l-1) (l+1):L]) ).^2.*p_next([1:(l-1) (l+1):L]) ) + ...
                                          sum( abs(ULMask(l,:)*G_SI*DLMask).^2.*omega_next) + sigma^2*norm(ULMask(l,:))^2 + ...
                                          sum( abs(ULMask(l,:)*kron((Err_H_UL/Nm).^(0.5),ones(Nm,1))).^2.*p_next) )];
    end
    
%     cons = [cons, psi_u_next(l)>=0];
    
%     cons = [cons, rcone(p_bar(l), p_next(l), 1/2)];
    cons = [cons, cone([p_bar(l), 0.5*(p_next(l)-1)], 0.5*(p_next(l)+1) )];
    cons = [cons, 2*sqrt(p_current(l))/psi_u_current(l)*p_bar(l)-p_current(l)/(psi_u_current(l)^2)*psi_u_next(l)>=lambda_u_next(l)];
    cons = [cons, p_bar(l)>=0];
    cons = [cons, lambda_u_next(l)>=0];
    
    
    
end

ULMask = ULMask_temp;
%% Cons. (45)



Coeff = zeros(K,M);
CoeffUL = zeros(L,M);
DLEff = sdpvar(1, M, 'full','real');
ULEff = sdpvar(1, M, 'full','real');
DLEffPow_sum = sdpvar(1, M, 'full','real');
ULEffPow_sum = sdpvar(1, M, 'full','real');
DLEffPow_next = sdpvar(1, M, 'full','real');
ULEffPow_next = sdpvar(1, M, 'full','real');

DLrateperAP = zeros(K,M);
ULrateperAP = zeros(L,M);
for m = 1:1:M
    
    if (HD && P_AP<10^(-3))
        continue
    end
%     Coeff{m} = zeros(1, K);
    maxNorm = 0;
    for k = 1:1:K
%         Coeff(k,m) = 1/nu + Pd/(norm(W_current([(m-1)*Nm+1:m*Nm],k))^2 + epspowDL);
%         Coeff(k,m) = Coeff(k,m)*norm(DLMask([(m-1)*Nm+1:m*Nm],k))^2;

%         maxNorm = max(norm(W_current([(m-1)*Nm+1:m*Nm],k))^2,maxNorm);
        DLrateperAP(k,m) = norm(H_d(k,:)*W_current(:,k))^2;
        Coeff(k,m) = norm(H_d(k,[(m-1)*Nm+1:m*Nm])*DLMask([(m-1)*Nm+1:m*Nm],k))^2;
    end
    
    for k = 1:1:K
%         Coeff(k,m) = Coeff(k,m)*(1/nu + Pd/(maxNorm + epspowDL));
%         Coeff(k,m) = Coeff(k,m)*(1/nu(1) + Pd/(DLrateperAP(k,m) + 0.01/M));
        Coeff(k,m) = Coeff(k,m)*(Pd/(DLrateperAP(k,m) + 0.01/M)) + 1/nu(1)*norm(DLMask([(m-1)*Nm+1:m*Nm],k))^2;
    end
    

    
%     if (~Init)
% %     cons = [cons, real(sum(Coeff{m}.*omega_current))<=xi_next(m)];
%         
%         cons = [cons, real(sum(Coeff(:,m)'.*omega_next))<=xi_next(m)+10^(-3)];
%         cons = [cons, real(sum(Coeff(:,m)'.*omega_next))>=xi_next(m)-10^(-3)];
% 
% 
%         DLEffPow_sum(m) = real(sum(Coeff(:,m)'.*omega_next));
%         
%         ULEffPow_sum(m) = sum(p_next.*real(diag(ULMask(:,[(m-1)*Nm+1:m*Nm])*H_u([(m-1)*Nm+1:m*Nm],:))'.^2));
%         
% 
%         DLEff(m) = DLEffPow_sum(m)/real(sum(Coeff(:,m)'.*omega_current)+epspow);
%         ULEff(m) = ULEffPow_sum(m)/(sum(p_current.*real(diag(ULMask(:,[(m-1)*Nm+1:m*Nm])*H_u([(m-1)*Nm+1:m*Nm],:))'.^2))+epspowUL); 
%         
% 
%     else
%         xi_next(m)=real(sum(Coeff(:,m)'.*omega_current));
        xi_next(m)=real(sum(Coeff(:,m)'.*omega_next));
        
%     end
    
%     cons = [cons, xi_next(m)<=P_AP/nu+Pd];
    maxNormUL = 0;
    for l = 1:1:L
%         maxNormUL = max(norm(ULMask(l,[(m-1)*Nm+1:m*Nm])*H_u([(m-1)*Nm+1:m*Nm],l))^2,maxNormUL);
        ULrateperAP(l,m) = norm(ULMask(l,:)*H_u(:,l))^2*p_current(l);
        CoeffUL(l,m) = norm(ULMask(l,[(m-1)*Nm+1:m*Nm])*H_u([(m-1)*Nm+1:m*Nm],l))^2;
    end
    
    for l = 1:1:L
%         CoeffUL(l,m) = CoeffUL(l,m)*(Pu/(maxNormUL + epspowUL));
        CoeffUL(l,m) = CoeffUL(l,m)*(Pu/(ULrateperAP(l,m) + 0.01/M));
    end
end


%% Cons. (46)
if (eta==0)
if (~(HD && P_AP<10^(-3)))
    
DLThput = sdpvar(1, K, 'full','real');
ULThput = sdpvar(1, L, 'full','real');

DLappx = log(1+lambda_d_current) - lambda_d_current./(1+lambda_d_current); 
for k = 1:1:K
     DLThput(k) = real(DLappx(k) + lambda_d_next(k)/(1+lambda_d_current(k)));
%     cons = [cons, DLappx(k) + lambda_d_next(k)/(1+lambda_d_current(k)) <= DLThput(k) ];
%     cons = [cons, DLThput(k) >= 0];
end

ULappx = log(1+lambda_u_current) - lambda_u_current./(1+lambda_u_current); 
for l = 1:1:L
    ULThput(l) = real(ULappx(l) + lambda_u_next(l)/(1+lambda_u_current(l)));
%     cons = [cons, ULappx(l) + lambda_u_next(l)/(1+lambda_u_current(l)) <= ULThput(l) ];
%     cons = [cons, ULThput(l)>=0];
end

SumThput = sum(DLThput) + sum(ULThput);

if (~Init)


%     cons = [cons, sum(xi_next.^2./(xi_current+epspow).*mu_current) + sum(p_next) + ...
%                   sum(xi_next./(xi_current+epspow).*mu_current*(Pa+Pu)) + sum((1-xi_next./(xi_current+epspow).*mu_current)*Ps) <= chi_next];

%     cons = [cons, sum(real(Coeff(:,m)'*omega_next').*mu_current') + sum(p_next) + ...
%                   sum(real(Coeff(:,m)'*omega_next').*mu_current'*(Pa+Pu)) + sum((1-real(Coeff(:,m)'*omega_next').*mu_current')*Ps) + SumThput*Pbh <= chi_next];
    
%     real(Coeff(:,m)'*omega_next')

%     cons = [cons, sum(real(Coeff'*omega_next').*mu_current') + sum(p_next) + sum(real(CoeffUL'*p_next').*mu_current') + ...
%                   sum(0.5*real(Coeff'*omega_next'+CoeffUL'*p_next').*mu_current'*Pa) + ...
%                   sum((1-0.5*real(Coeff'*omega_next'+CoeffUL'*p_next').*mu_current')*Ps) + SumThput*Pbh <= chi_next];

%     cons = [cons, sum(max([DLPowPercent' ULPowPercent'],[],2).*real(Coeff'*omega_next').*mu_current') + ...
%                   sum(p_next/nu(2)) + sum(max([DLPowPercent' ULPowPercent'],[],2).*real(CoeffUL'*p_next').*mu_current') + ...
%                   sum(max([DLPowPercent' ULPowPercent'],[],2).*mu_current'*Pa) + ...
%                   sum((1-max([DLPowPercent' ULPowPercent'],[],2).*mu_current')*Ps) + SumThput*Pbh/log(2) + ...
%                   M*PcirAP + (K+L)*PcirUE <= chi_next];

    cons = [cons, sum(max([DLPowPercent' ULPowPercent'],[],2).*real(Coeff'*omega_next').*mu_current') + ...
                  sum(p_next/nu(2)) + sum(max([DLPowPercent' ULPowPercent'],[],2).*mu_current'*Pu) + ...
                  sum(max([DLPowPercent' ULPowPercent'],[],2).*mu_current'*Pa) + ...
                  sum((1-max([DLPowPercent' ULPowPercent'],[],2).*mu_current')*Ps) + SumThput*Pbh/log(2) + ...
                  M*PcirAP + (K+L)*PcirUE <= chi_next];

    cons = [cons, chi_next>=0];

%     cons = [cons, chi_next<= 10^(10)*(M*P_AP + M*K*Pd + sum(P_UL) + M*Pu + M*Pa + M*Ps)];          
              
else


%     chi_next = sum(xi_next.^2./(xi_current+epspow).*mu_current)  +  sum(p_next) + ...
%                sum(xi_next./(xi_current+epspow).*mu_current*(Pa+Pu)) + sum((1-xi_next./(xi_current+epspow).*mu_current)*Ps);

%     chi_next = sum(real(Coeff'*omega_next').*mu_current') + sum(p_next) + ...
%                   sum(real(Coeff'*omega_next').*mu_current'*(Pa+Pu)) + sum((1-real(Coeff'*omega_next').*mu_current')*Ps) + SumThput*Pbh;

    chi_next = sum(max([DLPowPercent' ULPowPercent'],[],2).*real(Coeff'*omega_next').*mu_current') + ...
               sum(p_next/nu(2)) + sum(max([DLPowPercent' ULPowPercent'],[],2).*real(CoeffUL'*p_next').*mu_current') + ...
               sum(max([DLPowPercent' ULPowPercent'],[],2).*mu_current'*Pa) + ...
               sum((1-max([DLPowPercent' ULPowPercent'],[],2).*mu_current')*Ps) + SumThput*Pbh/log(2) + M*PcirAP + (K+L)*PcirUE;

end



end
end
            
% cons = [cons, chi_next<= M*P_AP + M*K*Pd + sum(P_UL) + M*Pu + M*Pa + M*Ps];
            
%% Cons. (48)


phi_bar = sdpvar(1, 1, 'full','real');
% cons = [cons, cone([1, 0.5*(phi_next-phi_bar)],0.5*(phi_next+phi_bar))];
% cons = [cons, eta*C_SE + (1-eta)*C_EE*(2/chi_current-chi_next/(chi_current^2))>=phi_bar];
% cons = [cons, phi_next>=0];
% cons = [cons, phi_bar>=0];
if (~(HD && P_AP<10^(-3)))
    if (eta==1)
        phi_next = 1;
    else
        if (~Init)
%         chi1 = eta*C_SE*chi_next + (1-eta)*C_EE;
% 
%         cons = [cons, chi_next-10^(-5) <= 2*phi_bar*phi_bar_current-phi_bar_current^2];
%         cons = [cons, cone([phi_bar, 0.5*(chi1-phi_next)],0.5*(chi1+phi_next))];
%         % cons = [cons, cone([sqrt(chi_current/(2*phi_bar_current))*phi_bar, sqrt(phi_bar_current/(2*chi_current))*chi_next, 0.5*(phi_next-1)], ...
%         %                     0.5*(phi_next+1)) ];
% 
%         cons = [cons, chi1>=0];                
%         cons = [cons, phi_bar>=0];                
        cons = [cons, phi_next>=0];
        cons = [cons, chi_next<=phi_next];
        else
%             chi1 = eta*C_SE*chi_next + (1-eta)*C_EE;
%             phi_bar = sqrt(real(chi_next));
%             phi_next = phi_bar^2/chi1;
            phi_next=chi_next;
        end
    end
end

%% Cons. (50)
% Rate_Threshold = 0;
pen_d = sdpvar(1, K, 'full','real');
% varrho = sdpvar(1, 1, 'full','real');
for k = 1:1:K
%     if (HD && P_AP<10^(-3))
%         pen_d(k)=0;
%         continue
%     end
%     cons = [cons, lambda_d_next(k)-exp(Rate_Threshold)+1>=0 ];

%     cons = [cons, lambda_d_next(k)>=varrho ];
%     cons = [cons, (eta*C_SE + (1-eta)*C_EE)*(lambda_d_next(k)-exp(Rate_Threshold)+1)>=pen_d(k) ];

    cons = [cons, (lambda_d_next(k)-exp(DLRate_Threshold)+1)>=pen_d(k) ];
    cons = [cons, pen_d(k)<=0];
end

pen_u = sdpvar(1, L, 'full','real');
for l = 1:1:L
%     if (HD && P_UL(l)<10^(-3))
%         pen_u(l)=0;
%         continue
%     end
%     cons = [cons, lambda_u_next(l)-exp(Rate_Threshold)+1>=0 ];

%     cons = [cons, lambda_u_next(l)>= varrho];

%     cons = [cons, (eta*C_SE + (1-eta)*C_EE)*(lambda_u_next(l)-exp(Rate_Threshold)+1)>= pen_u(l)];
    cons = [cons, (lambda_u_next(l)-exp(ULRate_Threshold)+1)>= pen_u(l)];
    cons = [cons, pen_u(l)<=0];
end


%% Objective func.
t_current = real((log(det( eye(K)+diag(lambda_d_current) )) + log(det( eye(L)+diag(lambda_u_current) )))/phi_current);

% obj = 0;

if (~Init)
    if (HD && P_AP<10^(-3))
        obj = logdet( eye(L)+diag(real(lambda_u_next)) );
    elseif (HD && P_UL(1)<10^(-3))
        obj = logdet( eye(K)+diag(real(lambda_d_next)) );
    else
        obj = logdet( eye(K)+diag(real(lambda_d_next)) ) + logdet( eye(L)+diag(real(lambda_u_next)) ) - t_current*real(phi_next);
    end
    
    obj = obj + penpar^(Iter)*(sum(pen_d) + sum(pen_u) );% + 2^(Iter)*sum(pen_mu);
    
% obj = logdet( diag(real(lambda_d_next)) );% + logdet( eye(L)+diag(real(lambda_u_next)) ) - t_current*real(phi_next);
% cons = [cons, logdet( eye(K)+diag(real(lambda_d_next)) ) + logdet( eye(L)+diag(real(lambda_u_next)) ) - t_current*real(phi_next)>=0];
else
%     if (HD && P_AP<10^(-3))
%         obj = 100*penpar*(sum(pen_u));% + 2^(Iter)*sum(pen_mu);
%     elseif (HD && P_UL(1)<10^(-3))
%         obj = 100*penpar*(sum(pen_d));% + 2^(Iter)*sum(pen_mu);
%     else
        obj = 100*penpar*(sum(pen_d) + sum(pen_u));% + 2^(Iter)*sum(pen_mu);
%     end
end



%% Optimization

%     cons = [cons, varrho>=10^(-5)];
%     myops = sdpsettings('solver','mosek');
    myops = sdpsettings('solver','sdpt3','verbose',0);
%     myops = sdpsettings('solver','sedumi','verbose',0);
    
    diagnotics = solvesdp(cons, -obj, myops);
%     diagnotics = solvesdp(cons, -varrho, myops)

%     diagnotics = solvesdp(cons, -(varrho+Sum_alpha+Sum_beta), myops)


%% Output
    Status.convergence = 0;
%     mu_next = ones(1,M);%double(xi_next)./(xi_current+10^(-5))

    obj = double(obj);
%     obj = double(varrho);
    DLRate_PerUser_next = log( ones(1,K)+double(lambda_d_next));
    DLRate_Sum = log(det( eye(K)+diag(double(lambda_d_next)) ));
    
    ULRate_PerUser_next = log( ones(1,L)+double(lambda_u_next));
    ULRate_Sum = log(det( eye(L)+diag(double(lambda_u_next)) ));
    
    if (HD)
        OptimalValue = DLRate_Sum+ULRate_Sum;
    else
        OptimalValue = (DLRate_Sum+ULRate_Sum)/double(phi_next);  %(eta + (1-eta)*10^4)*
    end
    omega = double(omega_next);
    p = double(p_next);
    
    W = DLMask*diag(omega.^(0.5));
    WW = zeros(1,M);
    DLPowPercent = zeros(K,M);
    
    A = H_u*diag(p.^(0.5));
    ULPowPercent = zeros(L,M);

    for m = 1:1:M
        Wm = W((m-1)*Nm+1:m*Nm,:);
        WW(m) = trace(Wm'*Wm);
        for k = 1:1:K
            Wmk = H_d(k,(m-1)*Nm+1:m*Nm)*W((m-1)*Nm+1:m*Nm,k);
            DLPowPercent(k,m) = 100*(norm(Wmk)^2)/(norm(H_d(k,:)*W(:,k))^2);
        end
        
        for l = 1:1:L
            Aml = ULMask(l,[(m-1)*Nm+1:m*Nm])*A([(m-1)*Nm+1:m*Nm],l);
            ULPowPercent(l,m) = 100*(norm(Aml)^2)/(norm(ULMask(l,:)*A(:,l))^2);%trace(Am'*Am);
        end

    end
    WW;
    DLPowPercent;
    ULPowPercent;
    DL = norm(H_d(k,:)*W(:,k))^2;
    UL = norm(ULMask(l,:)*A(:,l))^2;
    DLrateperAP;
    ULrateperAP;
    
    
%     mu_estimate
    
%     mu = max([double(DLEff); double(ULEff)],[],1)  %double(mu_next)

%Print
    lambda_d = double(lambda_d_next);
    lambda_u = double(lambda_u_next);
    pen_d = double(pen_d);
    pen_u = double(pen_u);
    psi_d = double(psi_d_next);
    psi_u = double(psi_u_next);
%EndPrint

    MaxMu = max([Coeff'*diag(omega), CoeffUL'*diag(p)],[],2)';
    Coeff;
    CoeffUL;
    
%     sum_u = zeros(1,L);
%     for l = 1:1:L
%         sum_u(l) = 1/(norm(ULMask(l,:)*H_u(:,l))^2)* ...
%                   ( sum( abs( ULMask(l,:)*H_u(:,(l+1):L) ).^2.*double(p_next((l+1):L)) ) + ...
%                     sum( abs(ULMask(l,:)*G_SI*DLMask).^2.*double(omega_next)) + sigma^2*norm(ULMask(l,:))^2 );
%     end
%     sum_u
%     pen_mu = double(pen_mu)
    if (~Init)
        
%     [mu, id] = max(double([DLEff; ULEff]),[],1) 
%     mu = double(xi_next)./(double(xi_next) + epspow)%mu_current + Step_mu.*Delta_mu
%     ActiveAP = max(double([DLEff; ULEff]),[],1)
    xi= double(xi_next);
%     sum_xi = sum(diag(double(omega_next)*Coeff))
%     Coeff
    phi = double(phi_next);
    phi_bar = double(phi_bar);
    chi = double(chi_next);
    obj1 = log(det( eye(K)+diag(double(lambda_d_next)) )) + log(det( eye(L)+diag(double(lambda_u_next)) )) - t_current*double(phi_next);
        if (obj1<0 && eta<1)
            disp('====> computing selection variables | obj1<0');
            PowScale.keepSolution = 1;
            PowScale.epsScale = prePowScale.epsScale + 1;
%             mu_next = double(xi_next)>epspow
            DLSel_next = [DLPowPercent>=pThres*100/M] + 10^(-5)*[DLPowPercent<pThres*100/M];
            ULSel_next = [ULPowPercent>=pThres*100/M] + 10^(-5)*[ULPowPercent<pThres*100/M];
            mu_next = floor(max(DLSel_next)+max(ULSel_next)-10^(-10));
            if ~(sum(abs(mu_next-mu_current))>0.5)
                Status.convergence = 1;
            end
            PowScale.reset = 1;
            Solution_current{5} = mu_next;
        else
            PowScale.keepSolution = 0;
            PowScale.epsScale = prePowScale.epsScale;
            PowScale.reset = 0;
            mu_next = mu_current;
            DLSel_next = DLSel_current;
            ULSel_next = ULSel_current;
        end
    else
        mu_next = mu_current;
        DLSel_next = DLSel_current;
        ULSel_next = ULSel_current;
        PowScale.keepSolution = 0;
        PowScale.epsScale = 0;
        PowScale.reset = 0;
    end
    
    if (PowScale.keepSolution)
        Solution_next = Solution_current;
    else
    Solution_next = {DLMask, ULMask, double(omega_next), ...
                                     double(p_next), ...
                                     mu_next, ... %mu,... %double(mu_next), ... %
                                     double(lambda_d_next), ...
                                     double(lambda_u_next), ...
                                     double(psi_d_next), ...
                                     double(psi_u_next), ...
                                     double(xi_next), ...
                                     {double(phi_next), double(phi_bar) }, ...
                                     double(chi_next), ...
                                     {DLSel_next, ULSel_next}, ...
                                     abs(sum(pen_d)+sum(pen_u))};
    
    end
    Status.info = diagnotics.info;


end

