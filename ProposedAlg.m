% % % Authors: Hieu V. Nguyen, Van-Dinh Nguyen, Octavia A. Dobre, Shree Krishna Sharma, Symeon Chatzinotas, Björn Ottersten, and Oh-Soon Shin, 
% % % Title: "On the Spectral and Energy Efficiencies of Full-Duplex Cell-Free Massive MIMO,"  
% % % IEEE Journal on Selected Areas in Communications (SI in Multiple Antenna Technologies for Beyond 5G), 2020


function [ OptValue, OptSolution, OptValueChain, DLRate_PerUser, ULRate_PerUser, isZeroPen ] = ProposedAlg( AllParameters, StartingPoint, strdisp, Method )
%PROPOSEDALG Summary of this function goes here
%   Detailed explanation goes here

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
Pbs = AllParameters{11};
P = AllParameters{12};
Pd = AllParameters{13};
Pu = AllParameters{14};
eta = AllParameters{15};
C_SE = AllParameters{16};
C_EE = AllParameters{17};
Rate_Threshold = AllParameters{18};

if (strfind(strdisp,'Convergence'))
    if (eta)
        MaxIteration = 30;
    else
        MaxIteration = 40;
    end
else
    MaxIteration = 30;
end



disp(['---------------------------------- ' strdisp ' - Run iterative algorithm --------------------------------']);

n = 0;
OptimalValue = -100;
OptimalValue_current = OptimalValue;
OptValueChain = [];
DLRate_PerUser = 0;
ULRate_PerUser = 0;
Solution_current = StartingPoint;

% while (n<100)
%    [ OptimalValue, DLRate_PerUser_next, ULRate_PerUser_next, Solution_next, Status] = ... 
%                     Get_optSolutionPerIteration(AllParameters, Solution_current, n); 
%         Solution_current = Solution_next;
% 
%         DLRate_PerUser = DLRate_PerUser_next;
%         ULRate_PerUser = ULRate_PerUser_next;
%         OptimalValue_current = OptimalValue
%         n=n+1;
% end

isZeroPen = 0;
Init = 1;
prePowScale.keepSolution = 0;
prePowScale.epsScale = 0;
prePowScale.reset = 0;
Backup_Point = {OptimalValue_current, Solution_current, OptValueChain, DLRate_PerUser, ULRate_PerUser, isZeroPen};
Backup_mu_current = Solution_current{5};

while (n<MaxIteration)

%     disp(['******************* ' strdisp ' --- Iteration: ' num2str(n+1) ' *********************']);
    
    % Solve 
    
    if (Method==8)
        [ OptimalValue, DLRate_PerUser_next, ULRate_PerUser_next, Solution_next, Status, PowScale] = ... 
                    Get_optSolutionPerIteration_blkmask(AllParameters, Solution_current, n, Init, prePowScale);
    else
    [ OptimalValue, DLRate_PerUser_next, ULRate_PerUser_next, Solution_next, Status, PowScale] = ... 
                    Get_optSolutionPerIteration(AllParameters, Solution_current, n, Init, prePowScale);
    end

    % Update


    if (isFeasible(OptimalValue_current, OptimalValue, Solution_current, Solution_next, Status.info, n))
        
        if (Init==0)
            Backup_Point = {OptimalValue_current, Solution_current, OptValueChain, DLRate_PerUser, ULRate_PerUser, isZeroPen};
        end
        Solution_current = Solution_next;

        DLRate_PerUser = DLRate_PerUser_next;
        ULRate_PerUser = ULRate_PerUser_next;

        prePowScale = PowScale;
        
        if (~Status.convergence && PowScale.reset)
            Init = 1;
            n = 0;
            mu_tmp = Solution_current{5};
            DLSel_tmp = Solution_current{13}{1};
            ULSel_tmp = Solution_current{13}{2};
            Solution_current = StartingPoint;
            Solution_current{5} = mu_tmp;
            Solution_current{13}{1} = DLSel_tmp;
            Solution_current{13}{2} = ULSel_tmp;
            Backup_mu_current = mu_tmp;
            continue
        end


        % Check convergence

%             if (abs(OptimalValue-OptimalValue_current)<10^-2)
        if (checkConvergence(((1-eta)*100 + eta)*OptValueChain, ((1-eta)*100 + eta)*OptimalValue, strfind(strdisp,'Convergence'), n) || Status.convergence)
            
            [mu_next, DLSel_tmp, ULSel_tmp] = varSelection(Solution_current, M, K, L, Nm, H_d, H_u);
            mu_next;
            mu_current = Solution_current{5};
            mu_next-mu_current;
            if (eta==0 && sum(abs(mu_next-Solution_current{5}))>0.5)
                Init = 1;
                n = 0;
                Solution_current = StartingPoint;
                Solution_current{5} = mu_next;
                Solution_current{13}{1} = DLSel_tmp;
                Solution_current{13}{2} = ULSel_tmp;
                Backup_mu_current = mu_next;
                continue
            else
                Solution_current{13}{1} = DLSel_tmp;
                Solution_current{13}{2} = ULSel_tmp;
                OptValueChain = [OptValueChain repmat(OptimalValue_current, 1, MaxIteration-length(OptValueChain) )];
                break;
            end
        else
            OptimalValue_current = OptimalValue;
            OptValueChain = [OptValueChain OptimalValue_current];

        end

    else

        disp('Infeasible point --> keeping the latest feasible point');
        OptValueChain = [OptValueChain repmat(OptimalValue_current, 1, MaxIteration-length(OptValueChain) )];
%         Solution_current = Backup_Point;
        OptimalValue_current = Backup_Point{1}*[Backup_Point{1}>0];
        Solution_current = Backup_Point{2};
        OptValueChain = Backup_Point{3};
        DLRate_PerUser = Backup_Point{4};
        ULRate_PerUser = Backup_Point{5};
        isZeroPen = Backup_Point{6};
        break;

    end

    n = n + 1;
    
    if (Solution_current{end}<10^(-5))
        Init = 0;
    end
    
    if (n==10)
        isZeroPen = [Solution_current{end}<10^(-5)];
        if (~isZeroPen)
            OptimalValue_current = 0;
            break
        end
    end

end

% disp('-------------------- Converging point ----------------------------');


OptValue = OptimalValue_current/log(2)*(OptimalValue_current>0);
OptSolution = Solution_current;
DLRate_PerUser = DLRate_PerUser/log(2);
ULRate_PerUser = ULRate_PerUser/log(2);



OptValueChain = OptValueChain/log(2);

% UplinkPower_PerUser = sum(p_current,2);
% Sum_Downlink_power = real(trace(W_current*W_current'));




end


function [mu_next, DLSel_next, ULSel_next] = varSelection(Solution_current, M, K, L, Nm, H_d, H_u)

% disp('====> computing selection variables');

pThres = 10^(-3);

DLMask = Solution_current{1};
ULMask = Solution_current{2};
omega = Solution_current{3};
p = Solution_current{4};

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
        Aml = ULMask(l,[(m-1)*Nm+1:m*Nm])*A((m-1)*Nm+1:m*Nm,l);
        ULPowPercent(l,m) = 100*(norm(Aml)^2)/(norm(ULMask(l,:)*A(:,l))^2);%trace(Am'*Am);
    end

end

DLSel_next = [DLPowPercent>=pThres*100/M] + 10^(-5)*[DLPowPercent<pThres*100/M];
ULSel_next = [ULPowPercent>=pThres*100/M] + 10^(-5)*[ULPowPercent<pThres*100/M];
mu_next = floor(max(DLSel_next)+max(ULSel_next)-10^(-10));

end

function [check] = isFeasible(OptimalValue_current, OptimalValue, Solution_current, Solution_next, Status, Step)

% check if the current point is feasible and the next point is infeasible

checkStatus = 1;
if (~isempty(findstr('Infeasible', Status)) && (Step>5))
    checkStatus = 0;
    disp('Out of feasible region --> Not pass range check');
end

checkAssignVar = 1;
% if (sum(sum(addVar_next{1}<0))>0 && sum(sum(addVar_next{2}<0))>0)
%     checkAssignVar = 0;
%     disp('Assignment Variables --> Not pass range check');
% end

checkObjVal = 1;

if (Step>1)
    if (Solution_current{end}<10^(-5) && Solution_next{end}>10^(-5))
        if (isempty(findstr('Successfully', Status)))
            checkObjVal = 0;
        end
    end
else
    if (abs(OptimalValue)>1000)
        checkObjVal = 0;
    end
end

% if ((OptimalValue_current>0 && OptimalValue<OptimalValue_current && Step>0) || (round(OptimalValue_current,20)==0) || (OptimalValue_current>1000))
%     checkObjVal = 0;
%     disp('OptimalValue --> Not pass range check');
%     OptimalValue_current
%     OptimalValue
% end

if ((OptimalValue_current>0 && OptimalValue<OptimalValue_current && Step>10) || (round(OptimalValue_current,20)==0) || (OptimalValue_current>1000))
    checkObjVal = 0;
    disp('OptimalValue --> Not pass range check');
    OptimalValue_current
    OptimalValue
end

checkDLPower = 1;
% if ((real(trace(W_current*W_current'))-real(trace(W_next*W_next')))/real(trace(W_current*W_current'))>0.1)
%     checkDLPower = 0;
% end

% check = ((Step<20)||(checkStatus && checkAssignVar && checkObjVal));
% check = ((Step<20) || (checkStatus && checkAssignVar && checkObjVal));
check = (checkStatus && checkAssignVar && checkObjVal && checkDLPower);



end


function [check] = checkConvergence(OptValueChain, OptimalValue, Converg, n)

check = 0;
% if (OptimalValue>80)
%     check=1;
% end

% if (n>2)
% if (abs(OptimalValue-OptValueChain(length(OptValueChain)))/abs(OptimalValue) < 0.01)
%     check = 1;
% end
% end

if (n>10)
    if (abs(OptimalValue-OptValueChain(length(OptValueChain))) < 10^-3)
        check = 1;
    end
    if (abs(OptimalValue-OptValueChain(length(OptValueChain)-5)) < 0.05)
        check = 1;
    end
    if (~Converg)
        if (abs(OptimalValue-OptValueChain(length(OptValueChain)))/abs(OptimalValue) < 0.05)
            check = 1;
        end
    else
        if (abs(OptimalValue-OptValueChain(length(OptValueChain)))/abs(OptimalValue) < 0.001)
            check = 1;
        end
    end
end

% if (check)
%     disp('<<<<< OptimalValue converges >>>>>');
% end

end
