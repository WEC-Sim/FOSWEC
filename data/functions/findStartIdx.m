% Dominic Forbush 7/13/2017
% Determine start of flap motion time from an uncut time series. Uses a
% maximum value threshold determined as thresh x the max of the "dead time" 
% test values at the start of the test.
%%%%%% 
% INPUTS
% time: the vector specifying time of samples (in seconds)
% fullTS: (size=[length(time) : k]) the data time series under  
%   consideration (equal length to time). If there are multiple time series 
%   to consider (k > 1), the operation will be carried out for each,
%   and the soonest start time will be output. 
% thresh: The multiple of max(abs(deadTime data)) at which useful data 
%   'begins'.
% leadNum: the number of wave cycles in advance of determined start of 
%   flap motion to include.
% waveT: the wave period, used to estimate the number of seconds included 
%   in leadNum.
% deadT: [OPTIONAL: default= 10 sec] the amount of dead time consistently
%   proceeding useful data
% OUTPUTS
% startIdx: the start index used for the trimmed time series
%%%% REVISION LOG
% Dominic Forbush 7/13/2017: added loop to handle k>1
% Dominic Forbush 7/17/2017: Threshold is now input argument. Modified
%   threshold to be the max(abs(deadTime data));

function startIdx=findStartIdx(time,fullTS, thresh, leadNum, waveT, deadT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for specified deadT inclusion
if nargin < 6
    deadT=10;
end

 % calculate start time for each time series under consideration
for k=1:length(fullTS(1,:));
% Find maximum of initial dead time
    [~,idxdT]=min(abs(time-deadT));
    maxdT=max(abs(fullTS(1:idxdT,k)));

% Approximate start of test; looking forward from minimum dead time
    FstartIdx= find(abs(fullTS(idxdT:end,k))> thresh*maxdT,1);
    FstartIdx=FstartIdx + idxdT;

% Step back specified number of wave cycles to find startIdx
    leadT=leadNum*waveT;
    [~,startIdxVec(k)]= min(abs(time-(time(FstartIdx)-leadT)));
end

%output the soonest start time
startIdx=min(startIdxVec);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





