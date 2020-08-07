%% Large-scale Fading

RadiusOfRegion = 1000;
% RadiusOfNearestUser = 10;
% PLx = 145.5;
% PLy = 20.9;

PL = 140.7;
d0 = 10/1000;
d1 = 50/1000;
Shd = 8;

% FirstTier = 2*RadiusOfCell*[1 0; cos(pi/3) sin(pi/3); -cos(pi/3) sin(pi/3); -1 0; -cos(pi/3) -sin(pi/3); cos(pi/3) -sin(pi/3)];

if (isempty(loadingfile))
    Parameters = {RadiusOfRegion, PL, d0, d1, Shd, 'S2S', 0}; % type = 'S2S'
    if (Method==3 || Method==6)
        positionAPs = [0 0];
        [ D_A2A, positionAPs, ~ ] = CreateLargeScaleFading( M, M, Parameters, positionAPs )
    else
        [ D_A2A, positionAPs, ~ ] = CreateLargeScaleFading( M, M, Parameters )
    end
    
    Parameters = {RadiusOfRegion, PL, d0, d1, Shd, 'S2D'}; % type = 'S2D' : AP --> DL users
    [ D_A2D, ~, positionDLUs ] = CreateLargeScaleFading( M, K, Parameters, positionAPs )
    [ D_A2U, ~, positionULUs ] = CreateLargeScaleFading( M, L, Parameters, positionAPs );
    [ D_U2D, ~, ~ ] = CreateLargeScaleFading( L, K, Parameters, positionULUs, positionDLUs );
else

    if (Method==3 || Method==6)
        positionAPs = [0 0];
        recalcD = 1;
    else
        load(loadingfile,'positionAPs');
        recalcD = 0;
    end
    load(loadingfile,'positionDLUs');
    load(loadingfile,'positionULUs');
    if (recalcD)
        Parameters = {RadiusOfRegion, PL, d0, d1, Shd, 'S2S', 0}; % type = 'S2S'
        [ D_A2A, positionAPs, ~ ] = CreateLargeScaleFading( M, M, Parameters, positionAPs )
        Parameters = {RadiusOfRegion, PL, d0, d1, Shd, 'S2D'}; % type = 'S2D' : AP --> DL users
        [ D_A2D, ~, positionDLUs ] = CreateLargeScaleFading( M, K, Parameters, positionAPs )
        [ D_A2U, ~, positionULUs ] = CreateLargeScaleFading( M, L, Parameters, positionAPs );
        [ D_U2D, ~, ~ ] = CreateLargeScaleFading( L, K, Parameters, positionULUs, positionDLUs );
    else
        load(loadingfile,'D_A2A');
        load(loadingfile,'D_A2D');
        load(loadingfile,'D_A2U');
        load(loadingfile,'D_U2D');
    end
end

Layout = Plot_Layout(RadiusOfRegion, {positionAPs, positionDLUs, positionULUs}, {'k^', 'rs', 'bo'}, {'AP', 'DLU', 'ULU'} )


APs_DLs = ones(M,K);
APs_ULs = ones(M,L);
if (Method==4 || Method==7)
    APs_DLs = zeros(M,K);
    AP_NoDLs = zeros(M,1);
    for k = 1:1:K
        DLPos_k = repmat(positionDLUs(k,:),M,1);
        Dist_k = sqrt(sum((DLPos_k-positionAPs).^2,2)+10^12*(AP_NoDLs>=Nm));
        [~, AP_selected] = min(Dist_k);
        APs_DLs(AP_selected,k) = 1;
        AP_NoDLs(AP_selected) = AP_NoDLs(AP_selected) + 1;
    end
%     AP_NoDLs'
%     sum(AP_NoDLs)
%     sum(APs_DLs)
%     D_A2D = D_A2D.*APs_DLs;
    
    APs_ULs = zeros(M,L);
    AP_NoULs = zeros(M,1);
    for l = 1:1:L
        ULPos_l = repmat(positionULUs(l,:),M,1);
        Dist_l = sqrt(sum((ULPos_l-positionAPs).^2,2)+10^12*(AP_NoULs>=Nm));
        [~, AP_selected] = min(Dist_l);
        APs_ULs(AP_selected,l) = 1;
        AP_NoULs(AP_selected) = AP_NoULs(AP_selected) + 1;
    end
%     AP_NoULs'
%     sum(AP_NoULs)
%     sum(APs_ULs)
%     D_A2U = D_A2U.*APs_ULs;
end



