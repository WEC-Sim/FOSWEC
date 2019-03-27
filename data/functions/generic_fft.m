%% Discrete fourier processing

%Description: code to perform spectral analysis of time series data

%Inputs:
%   - y: time domain series to transform to fourier space
%   - fs: sampling frequency

%Outputs:
%   - Pyy: power spectral density (units of variance)
%   - f: frequencies associated with Pyy

%Change Log:
%   - 2/7/2016 in-line code to reduce complexity of file set for data
%   processing

% Modified from code provided by Dr. Brian Polagye, University of 
% Washington (2018)

function [Pyy, f] = generic_fft(y, fs)

%detrend time series with a linear mean
y = detrend(y,'linear');

%calculate frequency bands
nfft = size(y,1);                       %length of transform - windows are defined as power of 2
f = fs/2*linspace(0,1,nfft/2+1)';       %frequency bands - 0:1:Nyquist frequency

%apply hamming window
%note 1: a rectangular window results in contamination of the spectrum
%   (Data Analysis Methods in Physical Oceanography p. 448)
%note 2: per Harris "On the use of windows for harmonic analysis with the discrete fourier transform" 
%   a Hamming window has preferable side lobes and worst case loss of
%   signal than a Hanning window (but probably doesn't matter for signal
%   processing in acoustic applications)
wts = repmat(hammWin(size(y,1),'spectral'),1,size(y,2));
y_wts = y.*wts;

S = fft(y_wts,nfft);
%first value is sum of all values (mean) - f = 0 - should be small
%first 1:nfft/2 are positive frequencies up to Nyquist
%second nfft/2:nfft are negative frequencies (complex conjugates of positive frequencies)

%single-sided power spectra (Data Analysis Methods in Physical Oceanography (p. 421), w/ modification
% 2 * complex magnitude(S) squared/(time interval -- number of points * sampling frequency)
Pyy = 2 * abs(S(1:nfft/2+1,:)).^2/(nfft*fs);
Pyy(1,:) = Pyy(1,:)/2;      %adjust power in 0 frequency band (mean - should be close to zero)
Pyy(end,:) = Pyy(end,:)/2;  %adjust power in Nyquist band

%adjust power spectral density to satisfy Parseval's theorem E = var(s)
rms_y = (mean(y.^2)).^0.5;

E = trapz(f,Pyy);      %total energy under spectra                 
corr = E./rms_y.^2;    %empirical correcetion factor

Pyy=Pyy./repmat(corr,size(Pyy,1),1);    %corrected spectral densities with variance preserved

%checksum = sum(abs(trapz(f,Pyy)-rms_y.^2));  %diagnostic: should be close to zero if no leakage

end