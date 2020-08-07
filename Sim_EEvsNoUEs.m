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

pathname = fileparts('./Figures/EE_vs_NoUEs/');
addpath(genpath('./Figures/EE_vs_NoUEs'));

tic

%% Simulation Setting

Simname = 'EEvsNoUEs_ext';

Files = 1:1;
NumOfSim = 10;
NoWorkers = 6;

Method = 1; % 1: FD-DM-MIMO
            % 2: FD-DM-MIMO with PCA
            % 3: FD-CM-MIMO
            % 4: FD-SC (small cell)
            % 5: HD-DM-MIMO
            % 6: HD-SC (small cell)


DLULSchemes = 3; % 1: DPC/ZF-SIC
                 % 2: ZF/ZF
                 % 3: MRT/MRC

Methodname = {'FD_DM_MIMO', 'FD_DM_MIMO_PCA', 'FD_CM_MIMO', 'FD_SC', 'HD_DM_MIMO', 'HD_SC'};
DLULname = {'DPC_SIC','ZF','MRT_MRC'};

% global MaxIteration
% MaxIteration = 50;



loadingfile = [];%['Layout6_EE.mat'];%
recalcD = 0;
RenewChannel = 1;



%% Set Parameters

ParameterSetting

URange = 10;
lenRange = length(URange);

%% Build Channel Model

% h = waitbar(0, [ str_Group ':  ' num2str(floor(0*100 / NumberOfRunning)) ' % Completed'],'Name','Percent of completion');
h = waitbar(0, 'Please wait ...','Name','Percent of completion');

rho = 10^(rho_dB/10);
    
iRate_Threshold = Rate_Threshold*log(2);

NoEntries = length(Files)*NumOfSim;

for iFile = Files
    
    filename = ['[' Simname num2str(floor(10*Rate_Threshold)) '] Method' num2str(Method) '_' DLULname{DLULSchemes} '_' num2str(iFile) '.mat'];
    
    savedname = fullfile(pathname, filename);

    Channels = cell(NumOfSim,lenRange);


    for iU = 1:1:lenRange

        K = URange(iU);
        L = URange(iU);
        Pl_dB = 23 * ones(L,1);
        Pl = 10.^(Pl_dB/10)*NORMAL*PowNorm;
        

        PathlossShadowingModel

        for iSim = 1:1:NumOfSim
            ChannelModel
            Channels{iSim,iU} = {H_DL, H_UL, G_AA_All, G_CCI};
        end

    end

% error('Layout is done !!!');

%% Run Algorithms



% figure;

% Chain = cell(NumOfSim,1);



% error('Layout is done !!!');



OptValue_All = zeros(NumOfSim,lenRange);
OptSolution_All = cell(NumOfSim,lenRange);
OptValueChain_All = cell(NumOfSim,lenRange);
DLRate_PerUser_All = cell(NumOfSim,lenRange);
ULRate_PerUser_All = cell(NumOfSim,lenRange);


iSim = 0;
while (iSim<NumOfSim)
    
    iSim = iSim + 1
    
    percent = ((find(iFile==Files)-1)*NumOfSim + iSim-1)/NoEntries;
    waitbar(percent, h,[ num2str(floor(percent*100)) ' % Completed'])

%     Clus = parcluster('local');
%     Clus.NumWorkers = NoWorkers;
% 
%     poolobj = parpool(Clus, Clus.NumWorkers);
    
%     parfor iU = 1:1:lenRange
    for iU = 1:1:lenRange
        
        
        iPbs = Pbs;
        
        K = URange(iU);
        L = URange(iU);
        Pl_dB = 23 * ones(L,1);
        Pl = 10.^(Pl_dB/10)*NORMAL*PowNorm;
        
        Blkrho = [];
        for m = 1:1:M
            Mat_rho = sqrt(rho)*ones(Nm, Nm)-1;
            Blkrho = blkdiag(Blkrho, Mat_rho);
        end
        Blkrho = Blkrho + 1;
        
        H_d = Channels{iSim,iU}{1,1};
        H_u = Channels{iSim,iU}{1,2};
        G_SI = Channels{iSim,iU}{1,3}.*Blkrho;
        GCCI = Channels{iSim,iU}{1,4};
        
%         Error_RD.error_up   = sqrt(Error_RD.error_up)*diag(D_Hu);
%         Error_RD.error_down = sqrt(Error_RD.error_down)*diag(D_Hd);
%         Error_RD.error_CCI  = sqrt(Error_RD.error_CCI).*D_glk;
    
    
    
        
    %% Running Algorithm
    
    disp('------------------------------ Initialization ---------------------------------');
    
    
    isFeasible = 0;
    NoTry = 0;
    while (~isFeasible)
        
        NoTry = NoTry + 1;
        
        if (Method<5)
        [StartingPoint] = GetInitialization( M, K, L, Nm, sigma, H_d, H_u, G_SI, GCCI, nu, Pd, Pu, eta, C_SE, C_EE, Pa, Ps, iPbs, Pl(1), [Method DLULSchemes], {APs_DLs, APs_ULs} );
    %     error('Layout is done !!!');
        AllParameters = {M, K, L, Nm, sigma, H_d, H_u, G_SI, GCCI, nu, iPbs, Pl, Pd, Pu, eta, C_SE, C_EE, iRate_Threshold, [Pa, PcirAP, PcirUE], Ps, Pbh, [DLULSchemes 0]};
        end
        
        disp('------------------------------- Optimization ----------------------------------');
%         MultiInterf = H_d*StartingPoint{1};
%         GSI_DLMask = G_SI*StartingPoint{1};
%         ASIC_GSI_DLMask = zeros(L,K);
%         ASIC_GSI_DLMask_Norm = zeros(1,L);
%         for l = 1:1:L
%             ASIC_GSI_DLMask(l,:) = StartingPoint{2}(l,:)*GSI_DLMask;
%             ASIC_GSI_DLMask_Norm(l) = norm(ASIC_GSI_DLMask(l,:))^2;
%         end

        if (Method<5)
            [OptValue, OptSolution, OptValueChain, DLRate_PerUser, ULRate_PerUser, isFeasible ] = ... 
                    ProposedAlg( AllParameters, StartingPoint, [], Method);
        
        else
            [StartingPoint] = GetInitialization( M, K, L, Nm, sigma, H_d, H_u, G_SI*0, GCCI*0, nu, Pd, Pu, eta, C_SE, C_EE, Pa, Ps, iPbs, 10^(-5), [Method DLULSchemes], {APs_DLs, APs_ULs} );
            AllParameters = {M, K, L, Nm, sigma, H_d, H_u, G_SI*0, GCCI*0, nu, Pbs, 10^(-5)*ones(length(Pl)), Pd, Pu, eta, C_SE, C_EE, [iRate_Threshold 0], [Pa, PcirAP, PcirUE], Ps, Pbh, DLULSchemes};
            [OptValueDL, OptSolutionDL, OptValueChainDL, DLRate_PerUser, ~, isFeasible ] = ... 
                    ProposedAlg( AllParameters, StartingPoint, [], Method);
            
            
                
            [StartingPoint] = GetInitialization( M, K, L, Nm, sigma, H_d, H_u, G_SI*0, GCCI*0, nu, Pd, Pu, eta, C_SE, C_EE, Pa, Ps, 10^(-5), Pl(1), [Method DLULSchemes], {APs_DLs, APs_ULs} );
            AllParameters = {M, K, L, Nm, sigma, H_d, H_u, G_SI*0, GCCI*0, nu, 10^(-5), Pl, Pd, Pu, eta, C_SE, C_EE, [0 iRate_Threshold], [Pa, PcirAP, PcirUE], Ps, Pbh, DLULSchemes};
            [OptValueUL, OptSolutionUL, OptValueChainUL, ~, ULRate_PerUser, isFeasible ] = ... 
                    ProposedAlg( AllParameters, StartingPoint, [], Method);
            
                
            OptValueDL;
            OptValueUL;
            OptValue = 0.5*(OptValueDL+OptValueUL);
            OptSolution = {OptSolutionDL, OptSolutionUL};
            OptValueChain = {OptValueChainDL, OptValueChainUL};
            DLRate_PerUser = 0.5*DLRate_PerUser;
            ULRate_PerUser = 0.5*ULRate_PerUser;
        
        end
        
%         [StartingPoint] = GetInitialization( M, K, L, Nm, sigma, H_d, H_u, G_SI, G_CCI, nu, Pd, Pu, eta, C_SE, C_EE, Pa, Ps, Pbs, Pl(1), [Method DLULSchemes] );
%     %     error('Layout is done !!!');
%         AllParameters = {M, K, L, Nm, sigma, H_d, H_u, G_SI, G_CCI, nu, Pbs, Pl, Pd, Pu, eta, C_SE, C_EE, iRate_Threshold, Pa, Ps};
% 
%         disp('------------------------------- Optimization ----------------------------------');
% %         MultiInterf = H_d*StartingPoint{1};
% %         GSI_DLMask = G_SI*StartingPoint{1};
% %         ASIC_GSI_DLMask = zeros(L,K);
% %         ASIC_GSI_DLMask_Norm = zeros(1,L);
% %         for l = 1:1:L
% %             ASIC_GSI_DLMask(l,:) = StartingPoint{2}(l,:)*GSI_DLMask;
% %             ASIC_GSI_DLMask_Norm(l) = norm(ASIC_GSI_DLMask(l,:))^2;
% %         end
% 
%         if (Method==1 || Method==2)
%             [OptValue, OptSolution, OptValueChain, DLRate_PerUser, ULRate_PerUser, isFeasible ] = ... 
%                     ProposedAlg( AllParameters, StartingPoint, []);
%                 
%         end
        
        if (NoTry==MaxTry)
            break;
        end
        
    end
    
    OptValue_All(iSim, iU) = OptValue*BW;
    OptSolution_All{iSim, iU} = OptSolution;
    OptValueChain_All{iSim, iU} = OptValueChain*BW;
    DLRate_PerUser_All{iSim, iU} = DLRate_PerUser*BW;
    ULRate_PerUser_All{iSim, iU} = ULRate_PerUser*BW;
    
    
    end

    
% delete(poolobj);
% 
% clear Clus    

   

% hold on

% plot([1:1:length(OptValueChain_All{iSim,1})], OptValueChain_All{iSim,1}, 'r-', 'linewidth', 2, 'markersize',9);

% Chain{iSim} = OptValueChain_All{iSim,1};

end


delete(poolobj);

clear Clus 


% save(savedname);

save(savedname, '-regexp', '^(?!h$).');

end

close(h);

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