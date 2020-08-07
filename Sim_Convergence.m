close all
clear
clc

pathCell = regexp(path, pathsep, 'split');
if (~any(strcmpi('.\yalmip', pathCell)))
    disp('_________________ Adding path ... _________________')
    addpath(genpath('./yalmip'));
    addpath(genpath('./SeDuMi_1_3'));
    addpath(genpath('./SDPT3-4.0'));
end
clear pathCell

addpath(genpath('./rician'));

pathname = fileparts('./Figures/Convergence/');
addpath(genpath('./Figures/Convergence'));

tic

%% Simulation Setting

NumOfSim = 1;

Method = 1; % 1: FD-DM-MIMO
            % 2: FD-DM-MIMO with PCA
            % 3: FD-CM-MIMO
            % 4: FD-SC (small cell)
            % 5: HD-DM-MIMO
            % 6: HD-CM-MIMO
            % 7: HD-SC (small cell)
            % 8: FD-DM-MIMO + MRT/MRC full power factors


DLULSchemes = 2; % 1: DPC/ZF-SIC
                 % 2: ZF/ZF
                 % 3: MRT/MRC

Methodname = {'FD_DM_MIMO', 'FD_DM_MIMO_PCA', 'FD_CM_MIMO', 'CRP_PF', 'Con_FD', 'FD_noNOMA', 'HD', 'OMA'};       

% global MaxIteration
% MaxIteration = 50;

filename = ['Convergence_' Methodname{Method} '_DLULSchemes' num2str(DLULSchemes) '_'  num2str(NumOfSim) '.mat'];

loadingfile = ['Layout8.mat'];%[];%
recalcD = 0;
RenewChannel = 0;

savedname = fullfile(pathname, filename);


%% Set Parameters

ParameterSetting

% if (Method==3 || Method==6)
%     Pbs = Pbs*M;
%     Nm = M*Nm;
%     M = 1;
% end

%% Build Channel Model

PathlossShadowingModel

% APs_DLs = ones(M,K);
% APs_ULs = ones(M,L);
% if (Method==4 || Method==7)
%     APs_DLs = zeros(M,K);
%     AP_NoDLs = zeros(M,1);
%     for k = 1:1:K
%         DLPos_k = repmat(positionDLUs(k,:),M,1);
%         Dist_k = sqrt(sum((DLPos_k-positionAPs).^2,2)+10^12*(AP_NoDLs>=Nm));
%         [~, AP_selected] = min(Dist_k);
%         APs_DLs(AP_selected,k) = 1;
%         AP_NoDLs(AP_selected) = AP_NoDLs(AP_selected) + 1;
%     end
% %     AP_NoDLs'
% %     sum(AP_NoDLs)
% %     sum(APs_DLs)
% %     D_A2D = D_A2D.*APs_DLs;
%     
%     APs_ULs = zeros(M,L);
%     AP_NoULs = zeros(M,1);
%     for l = 1:1:L
%         ULPos_l = repmat(positionULUs(l,:),M,1);
%         Dist_l = sqrt(sum((ULPos_l-positionAPs).^2,2)+10^12*(AP_NoULs>=Nm));
%         [~, AP_selected] = min(Dist_l);
%         APs_ULs(AP_selected,l) = 1;
%         AP_NoULs(AP_selected) = AP_NoULs(AP_selected) + 1;
%     end
% %     AP_NoULs'
% %     sum(AP_NoULs)
% %     sum(APs_ULs)
% %     D_A2U = D_A2U.*APs_ULs;
% end

% error('Layout is done !!!');

for iSim = 1:1:NumOfSim
    if (RenewChannel)
        ChannelModel
    else
        load(loadingfile,'H_DL');
        load(loadingfile,'H_UL');
        load(loadingfile,'G_AA_All');
        load(loadingfile,'G_CCI');
    end
    Channels{iSim,1} = {H_DL, H_UL, G_AA_All, G_CCI};
end

% error('Layout is done !!!');

%% Run Algorithms




% h = waitbar(0, [ str_Group ':  ' num2str(floor(0*100 / NumberOfRunning)) ' % Completed'],'Name','Percent of completion');
% h = waitbar(0, 'Please wait ...','Name','Percent of completion');

% Clus = parcluster('local');
% Clus.NumWorkers = 4;

% poolobj = parpool(Clus, Clus.NumWorkers);

figure;

Chain = cell(NumOfSim,1);

Blkrho = [];
for m = 1:1:M
    Mat_rho = sqrt(rho)*ones(Nm, Nm)-1;
    Blkrho = blkdiag(Blkrho, Mat_rho);
end
Blkrho = Blkrho + 1;

% error('Layout is done !!!');

OptValue_All = zeros(NumOfSim,1);
OptSolution_All = cell(NumOfSim,1);
OptValueChain_All = cell(NumOfSim,1);
DLRate_PerUser_All = cell(NumOfSim,1);
ULRate_PerUser_All = cell(NumOfSim,1);

iSim = 0;
while (iSim<NumOfSim)
    
    iSim = iSim + 1
    
    if (RenewChannel)
    
        H_d = Channels{iSim,1}{1,1};
        H_u = Channels{iSim,1}{1,2};
        G_SI = Channels{iSim,1}{1,3}.*Blkrho;
        G_CCI = Channels{iSim,1}{1,4};
        
%         Error_RD.error_up   = sqrt(Error_RD.error_up)*diag(D_Hu);
%         Error_RD.error_down = sqrt(Error_RD.error_down)*diag(D_Hd);
%         Error_RD.error_CCI  = sqrt(Error_RD.error_CCI).*D_glk;
    
    else
    
%     load('Parameter_convergence1.mat','Hd');
%     load('Parameter_convergence1.mat','Hu');
%     load('Parameter_convergence1.mat','G_SI');
%     load('Parameter_convergence1.mat','G_lk');
    
%     load('Convergence.mat','Hd');
%     load('Convergence.mat','Hu');
%     load('Convergence.mat','G_SI');
%     load('Convergence.mat','G_lk');

%         load(loadingfile,'Hd1');
%         load(loadingfile,'Hd2');
        load(loadingfile,'H_d');
        load(loadingfile,'H_u');
        load(loadingfile,'G_SI');
        load(loadingfile,'G_CCI');
    
    end
    
    count = 0;
    
%     percent = (iSim-1) / (NumOfSim);
%     waitbar(percent, h,[ num2str(floor(percent*100)) ' % Completed'])

% for g = Range_G
for i = length(rho_dB) % for HalfDuplex
% parfor i = 1:1:length(rho_dB)
    
    i_rho = rho_dB(i);
    

    rho = 10^(i_rho/10);
    
    iRate_Threshold = Rate_Threshold*log(2);
        
    %% Running Algorithm
    
    disp('------------------------------ Initialization ---------------------------------');
    
    
    isFeasible = 0;
    NoTry = 0;
    while (~isFeasible)
        
        NoTry = NoTry + 1;
        
        if (Method<5)
        [StartingPoint] = GetInitialization( M, K, L, Nm, sigma, H_d, H_u, G_SI, G_CCI, nu, Pd, Pu, eta, C_SE, C_EE, Pa, Ps, Pbs, Pl(1), [Method DLULSchemes], {APs_DLs, APs_ULs} );
    %     error('Layout is done !!!');
        AllParameters = {M, K, L, Nm, sigma, H_d, H_u, G_SI, G_CCI, nu, Pbs, Pl, Pd, Pu, eta, C_SE, C_EE, iRate_Threshold, [Pa, PcirAP, PcirUE], Ps, Pbh, DLULSchemes};
        end
        
        disp('------------------------------- Optimization ----------------------------------');
        if (Method<5)
        MultiInterf = H_d*StartingPoint{1};
        GSI_DLMask = G_SI*StartingPoint{1};
        end
%         ASIC_GSI_DLMask = zeros(L,K);
%         ASIC_GSI_DLMask_Norm = zeros(1,L);
%         for l = 1:1:L
%             ASIC_GSI_DLMask(l,:) = StartingPoint{2}(l,:)*GSI_DLMask;
%             ASIC_GSI_DLMask_Norm(l) = norm(ASIC_GSI_DLMask(l,:))^2;
%         end

        if (Method==1 || Method==2 || Method==3 || Method==4)
            [OptValue, OptSolution, OptValueChain, DLRate_PerUser, ULRate_PerUser, isFeasible ] = ... 
                    ProposedAlg( AllParameters, StartingPoint, ['Convergence'], Method);
        elseif (Method==5 || Method==6 || Method==7)
            [StartingPoint] = GetInitialization( M, K, L, Nm, sigma, H_d, H_u, G_SI*0, G_CCI*0, nu, Pd, Pu, eta, C_SE, C_EE, Pa, Ps, Pbs, 10^(-5), [Method DLULSchemes], {APs_DLs, APs_ULs} );
            AllParameters = {M, K, L, Nm, sigma, H_d, H_u, G_SI*0, G_CCI*0, nu, Pbs, 10^(-5)*ones(size(Pl,1),1), Pd, Pu, eta, C_SE, C_EE, [iRate_Threshold 0], [Pa, PcirAP, PcirUE], Ps, Pbh, DLULSchemes};
            [OptValueDL, OptSolutionDL, OptValueChainDL, DLRate_PerUser, ~, isFeasible ] = ... 
                    ProposedAlg( AllParameters, StartingPoint, ['Convergence'], Method);
            
            
                
            [StartingPoint] = GetInitialization( M, K, L, Nm, sigma, H_d, H_u, G_SI*0, G_CCI*0, nu, Pd, Pu, eta, C_SE, C_EE, Pa, Ps, 10^(-5), Pl(1), [Method DLULSchemes], {APs_DLs, APs_ULs} );
            AllParameters = {M, K, L, Nm, sigma, H_d, H_u, G_SI*0, G_CCI*0, nu, 10^(-5), Pl, Pd, Pu, eta, C_SE, C_EE, [0 iRate_Threshold], [Pa, PcirAP, PcirUE], Ps, Pbh, DLULSchemes};
            [OptValueUL, OptSolutionUL, OptValueChainUL, ~, ULRate_PerUser, isFeasible ] = ... 
                    ProposedAlg( AllParameters, StartingPoint, ['Convergence'], Method);
            
                
            OptValueDL
            OptValueUL
            OptValue = 0.5*(OptValueDL+OptValueUL)
            OptSolution = {OptSolutionDL, OptSolutionUL};
            OptValueChain = {OptValueChainDL, OptValueChainUL};
            DLRate_PerUser = 0.5*DLRate_PerUser;
            ULRate_PerUser = 0.5*ULRate_PerUser;
        elseif (Method==8)
            [StartingPoint] = GetInitialization_blkmask( M, K, L, Nm, sigma, H_d, H_u, G_SI, G_CCI, nu, Pd, Pu, eta, C_SE, C_EE, Pa, Ps, Pbs, Pl(1), [Method DLULSchemes], {APs_DLs, APs_ULs} );
    
            AllParameters = {M, K, L, Nm, sigma, H_d, H_u, G_SI, G_CCI, nu, Pbs, Pl, Pd, Pu, eta, C_SE, C_EE, iRate_Threshold, Pa, Ps, Pbh, DLULSchemes};
            
            [OptValue, OptSolution, OptValueChain, DLRate_PerUser, ULRate_PerUser, isFeasible ] = ... 
                    ProposedAlg( AllParameters, StartingPoint, ['Convergence'], Method);
        end
        
        
        
        if (NoTry==MaxTry)
            break;
        end
        
    end
    
    OptValue_All(iSim, i) = OptValue;
    OptSolution_All{iSim, i} = OptSolution;
    OptValueChain_All{iSim, i} = OptValueChain;
    DLRate_PerUser_All{iSim, i} = DLRate_PerUser;
    ULRate_PerUser_All{iSim, i} = ULRate_PerUser;
    
    
end
%     waitbar(i_NumOfSim / NumberOfRunning, h,[ str_Group ':  ' num2str(floor(i_NumOfSim*100 / NumberOfRunning)) ' % Completed'])
    
    

    
%     if (count==1)
%         HD_Rate(count) = mean(OptValue_HalfDuplex_Stat);
%         HD_Rate = repmat(HD_Rate(count), 1, length(rho_dB));
%     end


DLConnect = floor(sum(OptSolution{1,13}{1,1},2)')
ULConnect = floor(sum(OptSolution{1,13}{1,2},2)')
SumDLConnect = sum(DLConnect)
SumULConnect = sum(ULConnect)
NoAPsOff = length(find(OptSolution{1,5}==0))

hold on

plot([1:1:length(OptValueChain_All{iSim,1})], OptValueChain_All{iSim,1}, 'r-', 'linewidth', 2, 'markersize',9);

Chain{iSim} = OptValueChain_All{iSim,1};

figure
Layout = Plot_Layout(RadiusOfRegion, {positionAPs, positionDLUs, positionULUs, positionAPs(find(OptSolution{1,5}==0),:)}, ...
                                     {'k^', 'rs', 'bo', 'rx'}, {'AP', 'DLU', 'ULU', 'AP-Off'} )
                                 
figure
RateDis = sum(OptSolution{1,13}{1,1}) + sum(OptSolution{1,13}{1,2});
SleepID = find(RateDis<0.5);
ACtiveID = find(RateDis>=1);
Layout = Plot_Layout(RadiusOfRegion, {positionDLUs, positionULUs, positionAPs(SleepID,:) [positionAPs(ACtiveID,:), RateDis(ACtiveID)'] }, ...
                                     {'rs', 'bo', 'k^', 'rx'}, {'DLU', 'ULU', 'Sleep AP', 'Active AP'}, ...
                                     {{'none'}, {'none'}, {'none'}, ...
                                      {'ColorBar', [1 max(RateDis)], ...
                                      {'Weak Active', 'Strong Active'} } } )
                                 
figure
ModeIdx = sum(OptSolution{1,13}{1,1}) - sum(OptSolution{1,13}{1,2});
ModeDis = ModeIdx(find(OptSolution{1,5}>0));
Layout = Plot_Layout(RadiusOfRegion, {positionAPs, positionDLUs, positionULUs, [positionAPs(find(OptSolution{1,5}>0),:), ModeDis'] }, ...
                                     {'k^', 'rs', 'bo', 'rx'}, {'Sleep AP', 'DLU', 'ULU', 'HD/FD Mode'}, ...
                                     {{'none'}, {'none'}, {'none'}, ...
                                      {'ColorBar', [min(ModeDis), 0.5*(min(ModeDis) + max(ModeDis)), max(ModeDis)], ... 
                                      {'UL', 'FD', 'DL'} } }  )

end

% delete(poolobj);

% clear Clus

% close(h);

save(savedname);


%%

time = toc
hour = 0; min=0;
if (time>3600)
    hour = floor(time/3600);
    time = mod(time,3600);
end
if (time>60)
    min = floor(time/60);
    time = mod(time,60);
end
disp(['Running Time = ' num2str(hour) ' h ' num2str(min) ' m ' num2str(time) ' s.' ]);