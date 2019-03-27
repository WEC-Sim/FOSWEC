% Dominic Forbush 2017
% This function is an open-source, simpler approximation to Matlab's
% 'filtfilt' function. This allows the processing code to be run without
% using Matlab's signal processing toolbox.
%%%%%%
% INPUTS
% b: vector of numerator coefficients of the filter
% a: vector of denominator coefficients of the filter
% x: the time series to be filtered. Filter operates along the columns of x
% pad: [OPTIONAL, 1x1] The number of zeros with which we pad the ends of the
%   time series to avoid bad endpoint outputs. Default adds time series
%   length to either side of the filter.
% method: [OPTIONAL, 1x1] The method of padding. Method 1 is standard zero
%   padding. This is the default. Method 2 is a linear ramp pad of given 
%   length from 0 to the end points of the detrended time series. This 
%   greatly improves end point behavior for moving average filters, but is 
%   not appropriate for other filter types.

%%%%%%
% OUTPUTS
% filtx: the filtered time series

function filtx=filtFB(b,a,x,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-processing
% check for specified padding
if nargin<4 || isempty(varargin{1})
    pad=zeros(floor(length(x)),1);
else
    pad=zeros(varargin{1},1);
end

% check for specified method
if nargin<5 || isempty(varargin{2})
    method=1;
end

% determine size of x, preallocate;
[row,col]=size(x);
filtx=999.*ones(row,col);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% filter procedure
for k=1:col;
    % filter one column at a time
    x2filt=x(:,k);
    xlength=length(x2filt);
    % detrend signal by the raw mean
    xtrend=mean(x2filt);
    x2filt=x2filt-xtrend;
    
    % pad ends
    switch varargin{2}
        case 1
            padfront=pad;
            padback=pad;
            x2filt=[padfront;x2filt;padback];
        case 2
            padfront=linspace(0,x2filt(1),length(pad)).';
            padback=linspace(x2filt(end),0,length(pad)).';
            x2filt=[padfront;x2filt;padback];
    end
    % apply filter forwards
    x2filt=filter(b,a,x2filt);
    
    % reverse output time series
    x2filt=flipud(x2filt);
    
    % apply filter backwards;
    x2filt=filter(b,a,x2filt);
    
    % correct time series order
    x2filt=flipud(x2filt);
    
    % remove zero-padding
    x2filt=x2filt(length(pad)+1:end-(length(pad)));
    
    %retrend and specify output column
    filtx(:,k)=x2filt+xtrend;
end
end



