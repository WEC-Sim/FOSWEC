% Dominic Forbush 2017

% approximates inverse of student's t cumulative distribution function
% via table look-up. 
% INPUTS:
% p: desired probability to calculate, 0 <= p < 1;
% nu: degrees of freedom (generally sample size - 1)
% OUTPUTS:
% x: The factor to be multiplied by standard deviation of samples
%   such that mean +/- (x*std. dev) represents a p% confidence interval
%   around mean

% This function fails for large sample sizes nu. At those sample sizes you
% are probably safe to assume the distribution is normal.
% Not exactly correct very small values of nu (e.g., 1) at high
% probabilities (e.g., P > 0.99)

function x=tinvDom(p,nu)
dt=0.001;
tvec=[-10000:dt:10000]; % vector of potential t-values, brute force inverse
coeff= (gamma((nu+1)/2)/gamma(nu/2))*(1/sqrt(pi*nu)).*dt;
intArg=(1+(tvec.^2)./nu).^(-0.5*(nu+1));
Pval=coeff*cumtrapz(intArg);
[~,idx]=min(abs(Pval-p));
x=tvec(idx);

if x<0
    error('nu too large for this table-lookup. Assume distribution normal')
end
end
