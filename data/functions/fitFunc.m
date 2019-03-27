% function intended for identification of damped harmonic
% oscillator by iteration. Assumes canonical form of damped oscillator eqn
% models system dynamics, with linear damping.

% INPUTS
% params: the parameters of interest. A 2-element vector, natural frequency
%   (rad/s) is the first element. Damping ratio (zeta) is the second
%   element. Solver will minimize fit error by varying these parameters
% t: the vector of time stamps corresponding to data, y.
% y: the position data vector against which the modeled fit will be
%   compared.
% OUTPUTS:
% err: the root mean squared error between the fitted model and the
%   position data vector. This value will be minimized when this function
%   is called within an fminsearch procedure.

% Contribution from Ben Strom, bwstrom.com. Used with permission.

function err= fitFunc(params,t,y)

w=params(1);
b=params(2);

ytest=1-exp(-b.*w.*t).*(cos(w.*sqrt(1-b^2).*t)+(b/sqrt(1-b^2))...
    .*sin(w.*sqrt(1-b^2).*t));

% least-squared error criteria
err=sqrt(sum(ytest-y).^2); 