% Dominic Forbush 5/31/2018

% Creates a vector of coefficients corresponding to the hamming window of
% specified length

%INPUTS: 
% winLength: the length of the window to calculate
% purpose: either 'spectral' or 'filter'. If spectral, output N-1 points of
%   periodic window of N=winLength+1 points so that each adjacent window is
%   perfectly periodic. If filter, calculates a symmetric window of length
%   N = winLength, more suited for filter implementation.

%OUTPUTS:
% winVec: a vector of coefficients defined by the hamming window
%   for side-lobe attenuation. To window your time series:
%   windowedTS=winVec.*originalTS.

function winVec=hammWin(winLength,purpose);

switch purpose
    case 'spectral'
        N=winLength+1;
    case 'filter'
        N=winLength;
end

% specify for odd/even windows to ensure symmetry/periodicity
if ~rem(N,2) % even window
    halfWin=N/2;
    x=(0:halfWin-1)'/(N-1); % column vector
    winVec=0.54-0.46*cos(2*pi*x);
    winVec=[winVec; winVec(end:-1:1)];% uses back-to-back half windows
else % odd window
    halfWin=(N+1)/2;
    x=(0:halfWin-1)'/(N-1);
    winVec=0.54-0.46*cos(2*pi*x);
    winVec=[winVec; winVec(end-1:-1:1)]; % clips last point of half window
end 

if strcmp(purpose,'spectral')
    winVec(end)=[]; % delete last point of periodic window
end


end

