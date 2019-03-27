% Generate power-spectral-density via Welch's method of specified loads, 
% position and power for the desired RAO directory.

% This processing is appropriate for irregular wave cases. Load data were
% not trimmed or filtered in ./WECSIM2/Config#Irr; that processing is
% performed here for load data.

clear; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Directory information/setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

homeDir=cd;
RAOdir={'./final/Config4Irr'};                                              % Irregular wave Configuration directory
filName={'Config4Irr_final.mat'};                                           % File name of final data structure saved in above directory                                                                 % window for PSD determination, in seconds
loadQuery={'lcTyF1','lcFxF1','lcTyF2','lcFxF2','power'};                    % load fields to investigate. Also accepts 'power', to take PSD of power generated at the flaps.
savedir='./Config4Irr';                                                     % specify a save directory

% User inputs
ratel=100;                                                                  % load sampling rate (100 Hz)
ratew=50;                                                                   % wave gauge 6 sampling rate (50 Hz) from WaveTuning.mat cases
LfiltNum=0.125*ones(8,1) ;                                                  % FIR filter coefficients to be applied to load and position data (12.5 Hz cutoff)
WfiltNum=0.25*ones(4,1);                                                    % FIR filter coeff. to be app'd to wave data (12.5 Hz cutoff)
Tn=50;                                                                      % the number of peak periods to be considered in analysis
winN=15;                                                                    % length of spectral windows to use, in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% iterate through source directories
figIdx=1; % initialize figure count
for k=1:length(RAOdir)
    f1name=strrep(filName{k},'_final.mat',''); % the first field name of the loaded structure
    
    % enter directory, load data
    cd(RAOdir{k});
    data=load(filName{k});
    temp=fieldnames(data.(f1name).inter);
    i_min=[3000, 3000, 3000, 3000];                                         % visual estimate of start times for load TS
    w_min=round(i_min.*(ratew/ratel));                                      % corresponding start of wave TS
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% perform PSD calculation (Hamming window, 50% overlap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % iterate through trial numbers
    for k2=1:length(temp)
        f2name=strcat('Trial0', num2str(k2));                               % there are never more than 9 trials, so this syntax is robust
        Tq=data.(f1name).T(k2);                                             % peak period of this trial
        titStr=strcat('T_p=', num2str(Tq),'  H_m=', num2str(data.(f1name).H(k2)));
       
        % load wave data for this trial from calibration runs
        Wts=data.(f1name).final.WaveTS{k2,1}(w_min(k,k2):w_min(k,k2)+Tq*Tn*ratew);
        tw=[1/ratew:1/ratew:length(Wts)/ratew];
       
        % calculate PSD for wave time series
        [Pwwtemp,fw(:,k2),tw_win]=make_spectra(Wts, tw, winN*ratew, 0.5, ratew);
        Pww(:,k2)=median(Pwwtemp,2);
        
        % perform spectral analysis of the specified loads/powers
        for k3=1:length(loadQuery)
            if strcmp(loadQuery{k3},'power') || strcmp(loadQuery{k3},'Power') % perform power calculation and processing if specified in loadQuery
                
                % apply filter to flap position, calculate velocity
                filtPos1=filtFB(LfiltNum,1,data.(f1name).inter.(f2name).flapPosF1(i_min(k,k2):i_min(k,k2)+Tq*Tn*ratel)...
                    ,[],2)*pi/180;                                          % convert to radians
                filtPos2=filtFB(LfiltNum,1,data.(f1name).inter.(f2name).flapPosF2(i_min(k,k2):i_min(k,k2)+Tq*Tn*ratel)...
                    ,[],2)*pi/180;                                           % convert to radians
                vel1=diff(filtPos1)./(1/ratel);                            % Flap 1 velocity
                vel2=diff(filtPos2)./(1/ratel);                            % Flap 2 velocity
                T1=data.(f1name).inter.(f2name).lcTyF1(i_min(k,k2):i_min(k,k2)+(Tq*Tn*ratel)-1);
                T2=data.(f1name).inter.(f2name).lcTyF2(i_min(k,k2):i_min(k,k2)+(Tq*Tn*ratel)-1);
                            
                Power1=abs(vel1.*T1);                                       % Flap 1 power; device power production assumed agnostic to flap direction
                Power2=abs(vel2.*T2);                                       % Flap 2 power
                PowerTot=Power1+Power2;                                     % Combined flap powers
                tp=[1/ratel:1/ratel:length(Lts)/ratel];
                
                % Power PSD calculation
                [Pp1temp,fp(:,k3),tp_win]=make_spectra(Power1,tp,winN*ratel, 0.5, ratel); % Welch's method using hamming windows. 
                [Pp2temp,~,~]=make_spectra(Power2,tp,winN*ratel, 0.5, ratel);
                [Pp12temp,~,~]=make_spectra(PowerTot,tp,winN*ratel, 0.5, ratel);
                [Vel1temp,~,~]=make_spectra(vel1,tp,winN*ratel, 0.5, ratel);
                
                % median of spectral power taken over all windows
                Pp1=median(Pp1temp,2);                                      % flap 1 PSD power
                Pp2=median(Pp2temp,2);                                      % flap 2 PSD power
                PpTot=median(Pp12temp,2);                                   % combined flap PSD power
                VelS=median(Vel1temp,2);                                    % Flap 1 velocity PSD, for comparison
                
                % Make power PSD plots
                figure(figIdx); clf;
                semilogx(fp*2*pi,20*log10(Pp1), 'LineWidth',1.2,'Color','b');
                hold on; grid on;
                semilogx(fp*2*pi,20*log10(Pp2), 'LineWidth',1.2,'Color','r');
                semilogx(fp*2*pi,20*log10(PpTot), 'LineWidth',1.2,'Color','g');
                %legend('Flap 1 Power','Flap 2','Combined')                  % color code may be distorted, depending on Matlab version
                xlim([0 10])
                xlabel('Freq. (rad/s)')
                ylabel('\Phi_{PP}, dB')
                title(titStr)
                figIdx=figIdx+1;
                
                % legend color error bug work around
                figure(figIdx-1);
                p1=semilogy(NaN,NaN, '-b', 'LineWidth',1.2);
                p2=semilogy(NaN,NaN, '-r', 'LineWidth',1.2);
                p3=semilogy(NaN,NaN, '-g', 'LineWidth',1.2);
                legend([p1 p2 p3],{'Flap 1 Power', 'Flap 2', 'Combined'})
                % legend color error bug work around
                
            else                                                            % handles load query fields that are not power
                Lts=data.(f1name).inter.(f2name).(loadQuery{k3})(i_min(k,k2):i_min(k,k2)+Tq*Tn*ratel);
                tl=[1/ratel:1/ratel:length(Lts)/ratel];     
                [Plltemp,fl(:,k3),tl_win]=make_spectra(Lts, tl, winN*ratel, 0.5, ratel); % welchs method using hamming windows
                Pll(:,k3)=median(Plltemp,2);                                % saves each load as column, the median over all spectral windows
            end
        end
        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(figIdx); clf;
        subplot(1,2,1)                                                      % PSD of wave height
        semilogx(2*pi*fw(:,k2),20*log10(Pww(:,k2)),'LineWidth',1.2);
        grid on;
        xlabel('Freq. (rad/s)')
        ylabel('\Phi{\eta \eta}, dB');
        xlim([0 10])
        title(titStr)
        subplot(1,2,2)
        semilogx(fl*2*pi,20*log10(Pll),'LineWidth',1.2);                    % PSD of load
        grid on;
        xlabel('Freq. (rad/s)')
        ylabel('\Phi_{LL}, dB'); 
        xlim([0 10])
        legend(loadQuery);
        figIdx=figIdx+1;
        
        
    end
    clear temp data Pll Plltemp Pww Pwwtemp fw fl Lts Wts Pp1 Pp1temp Pp2 ...
        Pp2temp PpTot Pp12temp Power1 Power2 PowerTot tp tl tw vel1 vel2 T1 T2 filtPos1...
        filtPos2
end


