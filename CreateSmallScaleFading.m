function [Channel] = CreateSmallScaleFading( NoAntennasAtDestination, NoAntennasAtSource, dist)

% This function is to generate the channel
% NoSources : Number of source objects
% NoDestinations: Number of destination objects
% AntennaPerSource: Number of antennas per source (scalar or vector)
% AntennaPerDestination: Number of antennas per destination (scalar or vector)

if (nargin < nargin('CreateSmallScaleFading'))
    dist = {[]};
end

        
if (strcmp(dist{1},'Rice'))
    v = dist{2}*ones(NoAntennasAtDestination,NoAntennasAtSource); 
    s = dist{3}*ones(NoAntennasAtDestination,NoAntennasAtSource);
    Channel = sqrt(1/2)*(ricernd(v, s) + 1i*ricernd(v, s));
else
    Channel = sqrt(1/2)*(randn(NoAntennasAtDestination,NoAntennasAtSource) + 1i*randn(NoAntennasAtDestination,NoAntennasAtSource) );
end


	
	
end