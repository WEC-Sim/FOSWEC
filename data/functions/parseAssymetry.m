% Dominic Forbush 7/17/2017

% Checks for assymetry in oscillating flap motion. Removes if present.
%%%%% INPUTS
% peakLoc: vector output of peakfinder, index of local maxima.
% peakmag: vector output of peakfinder, magnitude of local maxima.
% minLoc: vector output of peakfinder, index of local minima.
% minmag: vector output of peakfinder, magnitude of local minima.
% win: number of mins/maxes to average over and compare to zero.
% thresh: the allowable mean difference from zero.
% lTS: the length of time series to parse. 
%%%%% OUTPUTS
% goodidx: outputs the indices between extrema over which the oscillation
%   is approx. zero-mean. Starts and ends on a maxima. 

function [goodidx]=parseAssymetry(peakMag,peakLoc,minMag,minLoc,thresh, win, lTS);

goodidx=1:1:lTS; % initialize output
numIt=min([length(peakLoc),length(minLoc)])-win+1; % needed iterations
badidx{numIt,1}=[]; % pre-allocate
    for k=1:numIt;
    winAvg=((sum(peakMag(k:k+win-1))./win)+(sum(minMag(k:k+win-1))./win))./2;
    if abs(winAvg)>thresh
        badidx{k,1}(:,1)=[min([peakLoc(k:k+win-1);minLoc(k:k+win-1)]):max([peakLoc(k:k+win-1);minLoc(k:k+win-1)])];
    end
    end
    % eliminate from 2nd-to-last maxima to end. Will ensure TS
    % ends on an maxima, not include 1/2 wave cycles
    badidx{numIt+1,1}(:,1)=peakLoc(end-1):lTS;
    
    % eliminate first wave cycle, ensure TS starts on maxima
    badidx{numIt+2,1}(:,1)=1:min(peakLoc);
    
    % convert cell structure to a matrix
    badidxmat=cell2mat(badidx);
    
    %delete repetitions
    badidxmat=unique(badidxmat);
    
    % eliminate bad indices from the output
    [~,~,ib]=intersect(badidxmat,goodidx);
    goodidx(ib)=[];
    
end