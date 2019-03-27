% Dominic Forbush 7/21/17

% Get hydrostatic stiffness fit as based on WECSIM phase 1 testing static
% flap 1 tests. The measured force is from a pulley, always applied to the
% end of the flap in a horizontal direction. Theta angle measured at flap
% base between flap centerline and vertical. The cosine of this angle times
% the flap length is the vertical distance(i.e., the lever arm) between the
% flap rotation axis and the application of the measured force.

% The fit relates displacement in radians to developed torque in N-m.

clear; 
% Values populated manually from WECSIM phase 1, median of steady-state
% measurements at each maintained flap displacement.
Theta=[3.278, 5.336, 6.919, 8.619, 10.38, 12.36, 14.16, 15.78, 17.56, 19.85, 21.11].*pi/180; % displacement (radians)
F=[7.638, 10.61, 13.27, 15.9, 18.51, 22.82, 25.37, 29.16, 34.1, 41.08, 47.42]; % measured force (Newtons)
R=0.58;                                                                     % flap length, (m)
Torque=F.* 0.58.*cos(Theta);                                                % torque on flap (N-m)

% Quadratic fit, constrained such that F=0 when displacement = 0 
param=[Theta.^2; Theta];
fitCoeff=mrdivide(Torque,param) ;
fitCoeff=[fitCoeff, 0];
% fit

figure(1);
clf
plot(Theta,Torque,'ko','MarkerFaceColor','k')
hold on
grid on
Theta_test=linspace(0,max(Theta)*1.1,100);
plot(Theta_test,polyval(fitCoeff,Theta_test),'--k')

grid on
xlabel('\Delta \Theta (rad)');
ylabel('Torque (N m)')

%% Calculate the moment of inertia about the flap hinge axis
M=23.14;                                                                    % flap mass in Kg
IcgYY=1.19;                                                                 % flap inertia in Kg m^2 from center of mass 
dist=0.46-0.29;                                                             % perp distance from IcgYY to flap axis
IYY=IcgYY + M*dist^2;                                                       % inertia about hinge