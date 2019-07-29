% Dominic Forbush 7/31/2017
% CompareBEMForcedOsc code compares linear damping, added mass, and
% excitation coefficients as a function of frequency between the experiment
% and BEM code
clear; close all; clc;

% directory setup
base_folder=pwd;
ForcOsc_folder='./final/ForcedOscillation';
BEMpathname='..\..\..\FOSWEC_Sims\sims\hydroData';
filename='flap_WAMIT.h5';

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
            Bexp(k,k2,k3)=ForcedOsc.final.CvTonly(idx(k3));                 % linear damping coeffs
            BexpStd(k,k2,k3)=ForcedOsc.final.CvTonlystd(idx(k3));           % linear damping est. standard deviations (cycle-to-cycle)
            Aexp(k,k2,k3)=ForcedOsc.final.AT(idx(k3));                      % added mass estimate
            AexpStd(k,k2,k3)=ForcedOsc.final.Astd(idx(k3));                 % added mass est. standard deviations
        end
        Bexp(k,k2)=mean(Bexp(k,k2,Bexp(k,k2,:)>0));                         % average statistics of repeated trials at a given freq and amp.
        BexpStd(k,k2)=mean(BexpStd(k,k2,Bexp(k,k2,:)>0));
        Aexp(k,k2)=mean(Aexp(k,k2,Aexp(k,k2,:)>0));
        AexpStd(k,k2)=mean(AexpStd(k,k2,Aexp(k,k2,:)>0));
    end
end
Bexp=Bexp(:,:,1);                                                           % utilize averages for plotting
BexpStd=BexpStd(:,:,1);
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
            cStr='o';
            fStr=[0 0.4470 0.7410];
        case 2
            cStr='s';
            fStr=[0.8500 0.3250 0.0980];
        case 3
            cStr='^';
            fStr=[0.929 0.694 0.125];
    end
    subplot(1,2,1)                                                          % plot experiment linear damping estimate
    ax(k2)=plot(wExpU(~isnan(Bexp(:,k2)))+off,Bexp(~isnan(Bexp(:,k2)),k2)...
        ,cStr,'MarkerFaceColor',fStr,'MarkerEdgeColor',fStr,'LineWidth',1.2);
    if k2==1;
        hold on
        grid on;
    end
        for k=1:length(wExpU)                                               % add error bars
            line([wExpU(k)+off,wExpU(k)+off],[Bexp(k,k2)+2*BexpStd(k,k2),...% Lines bound +/- 2 sigma
            Bexp(k,k2)-BexpStd(k,k2)],'Color',fStr)
        end
        
        
    subplot(1,2,2)                                                          % plot added mass estimate
    ax(k2)=plot(wExpU(~isnan(Aexp(:,k2)))+off,Aexp(~isnan(Aexp(:,k2)),k2)...
        ,cStr,'MarkerFaceColor',fStr,'MarkerEdgeColor',fStr,'LineWidth',1.2);
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
subplot(1,2,1)
ax(k2+1)=plot(wBEM,bBEM55,'s-k');                                           % BEM estimate of linear damping
xlabel('\omega [rad/s]')
ylabel('Total Damping, C_{Tot} [N m s]')
xlim([0 8]);

% added mass
subplot(1,2,2)
ax(k2+1)=plot(wBEM,aBEM55,'s-k');                                          % BEM estimate of added mass
xlabel('\omega [rad/s]')
ylabel('Added Mass, A [kg m^2]')
xlim([0 8]);
legend([ax],{'\Delta \Theta=10^o', '\Delta \Theta=15^o',...
    '\Delta \Theta=20^o','WAMIT'},'Location', 'South')

%% Create three fits 
lowAmpA=Aexp(:,1);
medAmpA=Aexp(:,2);
hiAmpA=Aexp(:,3);

% linear fit
lowAmpfit1=polyfit(wExpU(~isnan(lowAmpA)),lowAmpA(~isnan(lowAmpA)),1);
medAmpfit1=polyfit(wExpU(~isnan(medAmpA)),medAmpA(~isnan(medAmpA)),1);
hiAmpfit1=polyfit(wExpU(~isnan(hiAmpA)),hiAmpA(~isnan(hiAmpA)),1);
% fit errors
lAfErr(1)=norm((lowAmpA(~isnan(lowAmpA))-polyval(lowAmpfit1,wExpU(~isnan(lowAmpA)))));
mAfErr(1)=norm((medAmpA(~isnan(medAmpA))-polyval(medAmpfit1,wExpU(~isnan(medAmpA)))));
hAfErr(1)=norm((hiAmpA(~isnan(hiAmpA))-polyval(hiAmpfit1,wExpU(~isnan(hiAmpA)))));

% quad fit
lowAmpfit2=polyfit(wExpU(~isnan(lowAmpA)),lowAmpA(~isnan(lowAmpA)),2)
medAmpfit2=polyfit(wExpU(~isnan(medAmpA)),medAmpA(~isnan(medAmpA)),2)
hiAmpfit2=polyfit(wExpU(~isnan(hiAmpA)),hiAmpA(~isnan(hiAmpA)),2)
% fit errors
lAfErr(2)=norm((lowAmpA(~isnan(lowAmpA))-polyval(lowAmpfit2,wExpU(~isnan(lowAmpA)))));
mAfErr(2)=norm((medAmpA(~isnan(medAmpA))-polyval(medAmpfit2,wExpU(~isnan(medAmpA)))));
hAfErr(2)=norm((hiAmpA(~isnan(hiAmpA))-polyval(hiAmpfit2,wExpU(~isnan(hiAmpA)))));

% cubic fit
lowAmpfit3=polyfit(wExpU(~isnan(lowAmpA)),lowAmpA(~isnan(lowAmpA)),3)
medAmpfit3=polyfit(wExpU(~isnan(medAmpA)),medAmpA(~isnan(medAmpA)),3)
hiAmpfit3=polyfit(wExpU(~isnan(hiAmpA)),hiAmpA(~isnan(hiAmpA)),3)
% fit errors
lAfErr(3)=norm((lowAmpA(~isnan(lowAmpA))-polyval(lowAmpfit3,wExpU(~isnan(lowAmpA)))));
mAfErr(3)=norm((medAmpA(~isnan(medAmpA))-polyval(medAmpfit3,wExpU(~isnan(medAmpA)))));
hAfErr(3)=norm((hiAmpA(~isnan(hiAmpA))-polyval(hiAmpfit3,wExpU(~isnan(hiAmpA)))));

% quart fit
lowAmpfit4=polyfit(wExpU(~isnan(lowAmpA)),lowAmpA(~isnan(lowAmpA)),4)
medAmpfit4=polyfit(wExpU(~isnan(medAmpA)),medAmpA(~isnan(medAmpA)),4)
hiAmpfit4=polyfit(wExpU(~isnan(hiAmpA)),hiAmpA(~isnan(hiAmpA)),4)
% fit errors
lAfErr(4)=norm((lowAmpA(~isnan(lowAmpA))-polyval(lowAmpfit4,wExpU(~isnan(lowAmpA)))));
mAfErr(4)=norm((medAmpA(~isnan(medAmpA))-polyval(medAmpfit4,wExpU(~isnan(medAmpA)))));
hAfErr(4)=norm((hiAmpA(~isnan(hiAmpA))-polyval(hiAmpfit4,wExpU(~isnan(hiAmpA)))));

% quint fit (highest order possible)
lowAmpfit5=polyfit(wExpU(~isnan(lowAmpA)),lowAmpA(~isnan(lowAmpA)),5)
medAmpfit5=polyfit(wExpU(~isnan(medAmpA)),medAmpA(~isnan(medAmpA)),5)
hiAmpfit5=polyfit(wExpU(~isnan(hiAmpA)),hiAmpA(~isnan(hiAmpA)),5)
% fit errors
lAfErr(5)=norm((lowAmpA(~isnan(lowAmpA))-polyval(lowAmpfit5,wExpU(~isnan(lowAmpA)))));
mAfErr(5)=norm((medAmpA(~isnan(medAmpA))-polyval(medAmpfit5,wExpU(~isnan(medAmpA)))));
hAfErr(5)=norm((hiAmpA(~isnan(hiAmpA))-polyval(hiAmpfit5,wExpU(~isnan(hiAmpA)))));

figure; clf;
scatter([1:1:5],lAfErr)
hold on; grid on;
scatter([1:1:5],mAfErr)
scatter([1:1:5],hAfErr)
legend('Low Amp','Med Amp','High Amp')
xlabel('Polynomial Fit Order')
ylabel('Norm. Fit error')

% apply fits at BEM frequencys

AlowFit=polyval(lowAmpfit2,wBEM);
AmedFit=polyval(medAmpfit2,wBEM);
AhiFit=polyval(hiAmpfit2,wBEM);

%figure(1); 
%subplot(1,2,2)
%plot(wBEM,AlowFit,'--b')
%plot(wBEM,AmedFit,'--r')
%plot(wBEM,AhiFit,'--g')

save('AddedMassFits', 'AlowFit', 'AmedFit', 'AhiFit','wBEM')

