% Dominic Forbush 7/31/2017
% CompareBEMForcedOsc code compares linear damping, added mass, and
% excitation coefficients as a function of frequency between the experiment
% and BEM code

clear; close all; clc;

% directory setup
base_folder=pwd;
ForcOsc_folder='./final/ForcedOscillation';
BEMpathname='..\..\..\..\FOSWEC_Sims\sims\hydroData';
filename='flap_WAMIT_KR.h5';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load experiment data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(ForcOsc_folder);
load('ForcedOsc_final.mat');                                                % Forced oscillation data structure
wExp=2*pi./ForcedOsc.T;                                                     % Experimental radian frequencies evaluated
[wExpU,ia,ic]=unique(wExp);
[magU,id,ie]=unique(ForcedOsc.theta_targ);
for k=1:length(wExpU);
    goodIdx=find(k==ic);
    for k2=1:length(magU)
        goodIdx2=find(k2==ie);
        idx=intersect(goodIdx,goodIdx2);                                    % find oscillations at a given frequency and amplitude
        for k3=1:length(idx)
            Bexp(k,k2,k3)=ForcedOsc.final.CvTonly(idx(k3));                 % linear damping coeffs (assuming just linear damping)
            Bexpd(k,k2,k3)=ForcedOsc.final.CdTonly(idx(k3));                % quad. damping coeffs (assuming just quadratic damping)
            BexpStd(k,k2,k3)=ForcedOsc.final.CvTonlystd(idx(k3));           % linear damping est. standard deviations (cycle-to-cycle)
            BexpStdd(k,k2,k3)=ForcedOsc.final.CdTonlystd(idx(k3));          % quad damping est. standard deviations (cycle-to-cycle)
            Aexp(k,k2,k3)=ForcedOsc.final.AT(idx(k3));                      % added mass estimate
            AexpStd(k,k2,k3)=ForcedOsc.final.Astd(idx(k3));                 % added mass est. standard deviations
        end
        Bexp(k,k2)=mean(Bexp(k,k2,Bexp(k,k2,:)>0));                         % average statistics of repeated trials at a given freq and amp.
        Bexpd(k,k2)=mean(Bexpd(k,k2,Bexpd(k,k2,:)>0));                      
        BexpStd(k,k2)=mean(BexpStd(k,k2,Bexp(k,k2,:)>0));
        BexpStdd(k,k2)=mean(BexpStdd(k,k2,Bexpd(k,k2,:)>0));
        Aexp(k,k2)=mean(Aexp(k,k2,Aexp(k,k2,:)>0));
        AexpStd(k,k2)=mean(AexpStd(k,k2,Aexp(k,k2,:)>0));
    end
end
Bexp=Bexp(:,:,1);                                                           % utilize averages for plotting
Bexpd=Bexpd(:,:,1);
BexpStd=BexpStd(:,:,1);
BexpStdd=BexpStdd(:,:,1);
Aexp=Aexp(:,:,1);
AexpStd=AexpStd(:,:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load BEM data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(base_folder)
cd(BEMpathname)

% for forced osc tests only flap one was moving: load flap 1 coeff. data
bBEM=h5load(filename,'/body1/hydro_coeffs/radiation_damping/all');          % BEM radiation damping estimates
aBEM=h5load(filename,'/body1/hydro_coeffs/added_mass/all');                 % BEM added mass estimates
wBEM=h5load(filename,'/simulation_parameters/w');                           % BEM frequency domain

% Dimensionalize
rho=1000;                                                                   % fluid density kg/m^3
g=9.807;                                                                    % gravitational acceleration m/s^2
m= 23.14;                                                                   % flap mass (kg)

aBEM=aBEM*rho;                                                              % Corrects non-dimensionalization used in BEM output
bBEMsc=bBEM;
for k=1:length(wBEM)
    bBEMsc(:,:,k)=bBEM(:,:,k)*rho*wBEM(k);
end

bBEM=bBEMsc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate BEM comparison vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aBEM11=squeeze(aBEM(1,1,:));                                                % X-translation mode                                               
aBEM55=squeeze(aBEM(5,5,:));                                                % Flap pitching mode
aBEM51=squeeze(aBEM(5,1,:));                                                % cross-term (e.g., pitch added mass developed from x-translation mode)
bBEM11=squeeze(bBEM(1,1,:));
bBEM55=squeeze(bBEM(5,5,:));
bBEM51=squeeze(bBEM(5,1,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot for comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
off=0;                                                                      % a shift in x-direction to add so to facilitate comparison of different amplitude data
for k2=1:length(magU)
    switch k2                                                               % plot format selected based on amplitude
        case 1
            cStr='bo';
            fStr='b';
        case 2
            cStr='rs';
            fStr='r';
        case 3
            cStr='g^';
            fStr='g';
    end
    subplot(1,3,1)                                                          % plot experiment linear damping estimate
    ax(k2)=plot(wExpU(~isnan(Bexp(:,k2)))+off,Bexp(~isnan(Bexp(:,k2)),k2)...
        ,cStr,'MarkerFaceColor',fStr,'LineWidth',1.2);
    if k2==1;
        hold on
        grid on;
    end
        for k=1:length(wExpU)                                               % add error bars
            line([wExpU(k)+off,wExpU(k)+off],[Bexp(k,k2)+2*BexpStd(k,k2),...% Lines bound +/- 2 sigma
            Bexp(k,k2)-BexpStd(k,k2)],'Color',fStr)
        end
        
    subplot(1,3,2)                                                          % plot experiment quadratic damping estimate
    ax(k2)=plot(wExpU(~isnan(Bexpd(:,k2)))+off,Bexpd(~isnan(Bexpd(:,k2)),k2)...
        ,cStr,'MarkerFaceColor',fStr,'LineWidth',1.2);
    if k2==1;
        hold on
        grid on;
    end
        for k=1:length(wExpU)                                               % add error bars
            line([wExpU(k)+off,wExpU(k)+off],[Bexpd(k,k2)+2*BexpStdd(k,k2),...% Lines bound +/- 2 sigma
            Bexpd(k,k2)-BexpStdd(k,k2)],'Color',fStr)
        end   
        
    subplot(1,3,3)                                                          % plot added mass estimate
    ax(k2)=plot(wExpU(~isnan(Aexp(:,k2)))+off,Aexp(~isnan(Aexp(:,k2)),k2)...
        ,cStr,'MarkerFaceColor',fStr,'LineWidth',1.2);
    if k2==1;
        hold on
        grid on
    end
    for k=1:length(wExpU)                                                   % add error bars
        line([wExpU(k)+off,wExpU(k)+off],[Aexp(k,k2)+2*AexpStd(k,k2),...    % lines bound +/- 2 sigma 
        Aexp(k,k2)-2*AexpStd(k,k2)],'Color',fStr)                           
    end
    off=off+0.075;                                                          % adds offset for visualization
end


%% Add in experiment estimates from regression and numerical estimates 
% linear damping
subplot(1,3,1)
ax(k2+1)=plot(wExpU(1:end-1),fliplr(ForcedOsc.final.clin(2:end)),'d',...    % linear damping estimate from regression method
    'MarkerFaceColor',[0.5 0.5 0.5]);
ax(k2+2)=plot(wBEM,bBEM55,'s-k');                                           % BEM estimate of linear damping
xlabel('\omega [rad/s]')
ylabel('Total Damping, C_{Tot} [N m s]')
xlim([0 8]);

% quadratic damping
subplot(1,3,2)
ax(k2+1)=plot(wExpU(1:end),fliplr(ForcedOsc.final.cquad(1:end)),'d',...     % quad. damping estimate from regression method
    'MarkerFaceColor',[0.5 0.5 0.5]);                                       % there is no BEM estimate of quadratic damping
xlim([0 8]);
xlabel('\omega [rad/s]')
ylabel('Drag, C_D [N m s^2]')

% added mass
subplot(1,3,3)
ax(k2+1)=plot(wExpU(1:end),fliplr(ForcedOsc.final.Amass(1:end)),'d',...     % added mass damping from regression method (average over all oscillation amplitudes)
    'MarkerFaceColor',[0.5 0.5 0.5]);
ax2(k2+2)=plot(wBEM,aBEM55,'s-k');                                          % BEM estimate of added mass
xlabel('\omega [rad/s]')
ylabel('Added Mass, A [kg m^2]')
xlim([0 8]);
legend([ax],{'\Delta \Theta=10^o', '\Delta \Theta=15^o',...
    '\Delta \Theta=20^o','Regression','WAMIT'},'Location', 'South')





