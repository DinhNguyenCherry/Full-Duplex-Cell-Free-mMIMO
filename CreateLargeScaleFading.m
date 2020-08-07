function [ D_S2D, positionSources, positionDests ] = CreateLargeScaleFading( NoSources, NoDestinations, Parameters, posSources, posDests )
%CREATED Summary of this function goes here
%   Detailed explanation goes here

Radius = Parameters{1};
% NearestS2D = Parameters{2}; % nearest distance from sources to destinations
% PLx = Parameters{2};
% PLy = Parameters{3};
L = Parameters{2};
d0 = Parameters{3};
d1 = Parameters{4};
Shd = Parameters{5};
type = Parameters{6};

if (length(Parameters)>6)
    DistLim = Parameters{7};
else
    DistLim = 0;
end

while (true)

if (nargin<nargin('CreateLargeScaleFading')-1)
    
    rvectorS = Radius*rand(1,NoSources);
    anglevectorS = 2*pi*rand(1,NoSources);
    positionSources = [(rvectorS.*cos(anglevectorS))' (rvectorS.*sin(anglevectorS))'];
    
    if (findstr(type,'S2S'))
        positionDests = positionSources;
    else
    
        rvectorD = Radius*rand(1,NoDestinations);
        anglevectorD = 2*pi*rand(1,NoDestinations);
        positionDests = [(rvectorD.*cos(anglevectorD))' (rvectorD.*sin(anglevectorD))'];
    
    end
    
    
elseif (nargin<nargin('CreateLargeScaleFading'))
    
    if (size(posSources,1)<NoSources)
        NinpSources = size(posSources,1);
        rvectorS = Radius*rand(1,NoSources-NinpSources);
        anglevectorS = 2*pi*rand(1,NoSources-NinpSources);
        positionSources = [(rvectorS.*cos(anglevectorS))' (rvectorS.*sin(anglevectorS))'];
        positionSources = [posSources; positionSources];
    else
        positionSources = posSources(1:NoSources,:);
    end
    
    if (findstr(type,'S2S'))
        positionDests = positionSources;
    else
    
        rvectorD = Radius*rand(1,NoDestinations);
        anglevectorD = 2*pi*rand(1,NoDestinations);
        positionDests = [(rvectorD.*cos(anglevectorD))' (rvectorD.*sin(anglevectorD))'];
    
    end
    
else
    
    if (size(posSources,1)<NoSources)
        NinpSources = size(posSources,1);
        rvectorS = Radius*rand(1,NoSources-NinpSources);
        anglevectorS = 2*pi*rand(1,NoSources-NinpSources);
        positionSources = [(rvectorS.*cos(anglevectorS))' (rvectorS.*sin(anglevectorS))'];
        positionSources = [posSources; positionSources];
    else
        positionSources = posSources(1:NoSources,:);
    end
    
    if (findstr(type,'S2S'))
        positionDests = positionSources;
    else
        
        if (size(posDests,1)<NoDestinations)
            NinpDests = size(posDests,1);
            rvectorD = Radius*rand(1,NoDestinations-NinpDests);
            anglevectorD = 2*pi*rand(1,NoDestinations-NinpDests);
            positionDests = [(rvectorD.*cos(anglevectorD))' (rvectorD.*sin(anglevectorD))'];
            positionDests = [posDests; positionDests];
        else
            positionDests = posDests(1:NoDestinations,:);
        end
    
    end
    
end


    position_repDests = kron(positionDests(1:NoDestinations,:),ones(NoSources,1));
    position_repSours = repmat(positionSources(1:NoSources,:), NoDestinations, 1);
    
    distance_S2D = sqrt(sum((position_repDests-position_repSours).^2 ,2));
    
    
    if (strfind(type,'S2S'))
        checkMat = reshape(distance_S2D, NoSources, NoDestinations);
        checkMat = [checkMat<DistLim]-eye(NoSources);
        if (sum(sum(checkMat))>0)
            continue
        end
    else
%         DistLim
%         Distcomp = [distance_S2D<DistLim]
        if (~isempty(find(distance_S2D<DistLim,1)))
            continue
        end
    end
    
    
    
%     PL_S2D = GetPathloss(PLx, PLy, distance_S2D/1000); % distance_S2D in km
    PL_S2D = GetPathloss(L, d0, d1, Shd, distance_S2D/1000); % distance_S2D in km
    
    if (findstr(type,'S2S'))
%         idx_self = find(PL_S2D==Inf | PL_S2D==NaN | PL_S2D==-Inf);
%         PL_S2D
        idx_self = find(isnan(PL_S2D));
        PL_S2D(idx_self) = 0;
    end
    
    D_S2D = (reshape(PL_S2D, NoSources, NoDestinations)).^(0.5);
    
    break
    
end

end


function [Pathloss] = GetPathloss(L, d0, d1, Shd, distance)

% distance (km)

% Pathloss_dB = firstPar + secondPar*log10(distance);


% distance (m)

c0 = [((d0-distance)./abs(d0-distance+10^(-20)))>=0];
c1 = [((d1-distance)./abs(d1-distance+10^(-20)))>=0];

Pathloss_dB = -L - 20 * ( c0*log10(d0) + (1-c0).*log10(distance) ) ...
                 - 15 * ( c1*log10(d1) + (1-c1).*log10(distance) );
             
ShdCheck = [distance>d1];
Pathloss_dB = Pathloss_dB - Shd*rand(size(distance)).*ShdCheck;

Pathloss = 10.^(Pathloss_dB/10);

end

