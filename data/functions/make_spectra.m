

% Description: process data with a DFT to produce a windowed "spectrogram"  
% of input data via Welch's method. 
% This procedure is appropriate for data sampled at uniform intervals.

% Modified from code originally written by Dr. Brian Polagye, University of
% Wasington (2018)

% INPUTS
% u: the data vector to process
% t_input: the time-stamps of data vector u. Sample rate must be uniform!
% win_pts: the number of samples over which to window data for spectrogram
%   calculation.
% win_overlap: the amount to overlap consecutive windows, as a decimal
%   (e.g., for 50% overlap, win_overlap=0.5)
% fs: the UNIFORM sampling rate, in Hz.

% OUTPUTS
% Puu: the windowed spectral data. Each column is the spectral content for
%   an individual window. Each row is the spectral content at a given
%   frequency, corresponding to vector f.
% f: the frequency vector, corresponding to the spectral content of Puu.
%   The number of rows of Puu = length(f).
% t: the time-stamps of the individual windows (to aid visualization of
%   window overlap). Has units of t_input.

function [Puu, f, t] = make_spectra(u, t_input, win_pts, win_overlap, fs)

offset = 10;                                                            

%buffer input
t_in = buffer(t_input,win_pts,floor(win_pts*win_overlap))+offset;          % buffered input time index, offset prevents zero-wrapping
u_in = buffer(u,win_pts,floor(win_pts*win_overlap));                       % buffered velocity

%identify zero-padded buffers and prune
full_buffer = find(sum(t_in==0)==0);
t_in = t_in(:,full_buffer)-offset;                                         % subtract offset
u_in = u_in(:,full_buffer);

t = mean(t_in,1);                                                          % reset time stamps to correspond to buffered data

%calculate spectra for each sample interval
for i = 1:size(u_in,2)
    [Puu(:,i), f] = generic_fft(u_in(:,i), fs);                            % fft procedure. uses Hamming window.
    
    %initialize return matrix for subsequent intervals
    if i == 1
        Puu = horzcat(Puu(:,i),zeros(size(Puu,1),size(u_in,2)-1));
    end
end

end

