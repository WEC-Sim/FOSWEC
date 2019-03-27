% CompareBEMexc plots excitation coefficients from regular wave testing,
% irregular wave testing, and computational boundary element estimation.

% User inputs specify regular wave, irregular wave, or both, plotting, the
% field names fo the data to plot, and the desired plot titles.

clear; close all; clc;

% Directory information
RegWdir='./final/WaveExcitationReg';
IrregWdir='./final/WaveExcitationIrr';
expPathreg=pwd;

% User Inputs
queryName={'cFxF1','cTyF1','cFxF2','cTyF2'};                               % a field name specifying the phase/mag plot to present (accepts a vector)
titStr={'f_{Fx} F1','f_{Ty} F1','f_{Fx} F2','f_{Ty} F2'};                  % string vector of plot titles correspoding to the above
plotWhat=[1 2];                                                            % 1 to plot regular wave, 2 to plot irregular waves, [1,2] or [1:2] for both
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load BEM data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BEMpathname='..\..\..\..\FOSWEC_Sims\sims\hydroData';
filename='flap_WAMIT_KR.h5';                                                % This is a BEM run of the Flap with the origin at the hinge
cd(BEMpathname)

% extract data of interest
flap1reL=h5load(filename,'/body1/hydro_coeffs/excitation/re');
flap1imL=h5load(filename,'/body1/hydro_coeffs/excitation/im');
flap1magL=h5load(filename,'/body1/hydro_coeffs/excitation/mag');
wBEM=h5load(filename,'/simulation_parameters/w');                          % frequency domain of BEM calculation

% Dimensionalize
rho=1000;                                                                   % fluid density kg/m^3
g=9.807;                                                                    % gravitational acceleration m/s^2

% specify outputs to plot
headingNum=1;                                                               % the wave direction heading desired (inspect H5 file)
for k=1:length(queryName)                                                   % the mode numbers under consideration (a vector: [Fx Fy Fz Tx Ty Tz])
    if isempty(strfind(queryName{1,k},'Fx'))==0
        modeNum(k)=1;
        axisStr{k}='Magnitude (N/m)';                                       % the axis label and unit for given query
    elseif isempty(strfind(queryName{1,k},'Fy'))==0
        modeNum(k)=2;
        axisStr{k}='Magnitude (N/m)';
    elseif isempty(strfind(queryName{1,k},'Fz'))==0
        modeNum(k)=3;
        axisStr{k}='Magnitude (N/m)';
    elseif isempty(strfind(queryName{1,k},'Tx'))==0
        modeNum(k)=4;
        axisStr{k}='Magnitude (N-m/m)';
    elseif isempty(strfind(queryName{1,k},'Ty'))==0
        modeNum(k)=5;
        axisStr{k}='Magnitude (N-m/m)';
    elseif isempty(strfind(queryName{1,k},'Tz'))==0
        modeNum(k)=6;
        axisStr{k}='Magnitude (N-m/m)';
    end
    
    % create data matrices
    flap1re(:,k)=squeeze(flap1reL(modeNum(k),headingNum,:));
    flap1im(:,k)=squeeze(flap1imL(modeNum(k),headingNum,:));
    flap1mag(:,k)=squeeze(flap1magL(modeNum(k),headingNum,:));
    
    % calculate phase
    flap1ph(:,k)= atan2(flap1im(:,k),flap1re(:,k));
    
end

% Dimensionalize
flap1mag=flap1mag.*rho.*g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load experiment data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% regular waves
cd(expPathreg);
cd(RegWdir);
load('WaveExcReg_final.mat')
wReg=(2*pi./WaveExcReg.T).';                                                % regular wave frequency domain

% Irregular waves
cd(expPathreg);
cd(IrregWdir);
load('WaveExcIrr_final.mat')
freqw=WaveExcIrr.final.FreqDomain.freqw;                                   % irregular wave frequency domain

colorvec={'b','r','g','m','c'};                                             % color coding by wave height

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot irregular waves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if max(plotWhat)==2
    % extract averaged spectra, magnitude, and phase
    [Hunq,ia,ic]=unique(WaveExcIrr.Hm0(WaveExcIrr.Flag < 1));                   % finds unique wave heights of unflagged irregular cases
    for k3=1:length(queryName);
        figure(2*(k3-1)+1); clf;                                                % a new figure generated for each element of queryName vector
        for k=1:length(Hunq)
            goodidx=find(ic==k);                                                % creates temporary variables to plot
            avfreqw=WaveExcIrr.final.FreqDomain.freqAvgw{k};                    % band-averaged spectra
            avmeanMagSpec=WaveExcIrr.final.(queryName{k3}).MagAvg{1,k};         % weighted averaged, band averaged magnitude
            avmeanPhSpec=WaveExcIrr.final.(queryName{k3}).phAvg{1,k};           % weighted averaged, band averaged phase
            plot(avfreqw,avmeanMagSpec,'LineWidth',1.2,'Color',colorvec{k},...
                'DisplayName',strcat('IrrH=',num2str(Hunq(k)))); % plot irregular wave excitation mag., height-by-height
            ax(k)=gca;
            if k==1;                                                            % apply axis labels and grid on the first iteration
                hold on; grid on;
                xlabel('Frequency (rad/s)')
                ylabel(axisStr{k3})
                figure(2*(k3-1)+2); clf;
            end
            figure(2*(k3-1)+2)
            polarplot(avmeanPhSpec,avfreqw,'o','Color',colorvec{k},...          % plot irregular wave excitation phase, height-by-height
                'MarkerFaceColor',colorvec{k},'DisplayName',strcat('IrrH=',num2str(Hunq(k))));
            ax2(k)=gca;
            if k==1;                                                            % limit r-axis to meaningful frequency range
                hold on; grid on;
                rlim([0 10]);
            end
            
            figure(2*(k3-1)+1);
            clear avfreqw avmeanMagSpec avmeanPhSpec                            % clear temporary variables
        end
        figure(2*(k3-1)+1)
        title([titStr{k3}])
        xlim([0 10])
        figure(2*(k3-1)+2)
        title(['\Phi' titStr{k3}])
        rlim([0 10])
        legidx=k;
    end
else
    legidx=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot regular wave cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if min(plotWhat)==1
    [Treg,ia,ic]=unique(WaveExcReg.T);
    [Hreg,id,ie]=unique(WaveExcReg.H);
    for k3=1:length(queryName);                                                 % adds regular wave data to the plots created above
        for k=1:length(Hreg)                                                    % goes height-by-height.
            goodidx=find(ie==k);
            
            figure(2*(k3-1)+1);                                                 % regular wave magnitude plotted
            ax3(legidx+k)=plot(wReg(goodidx),WaveExcReg.final.(strcat(queryName{k3},'mag'))(goodidx),...
                '^','Color',colorvec{ie(goodidx(1))},'MarkerFaceColor',colorvec{ie(goodidx(1))},...
                'DisplayName',strcat('RegH=',num2str(Hreg(k))));
            if k==1;                                                            % apply axis labels and grid on the first iteration
                hold on; grid on;
                xlabel('Frequency (rad/s)')
                ylabel(axisStr{k3})
            end
            
            figure(2*(k3-1)+2)                                                  % regular wave phase plotted
            ax2(legidx+k)=polarplot(WaveExcReg.final.(strcat(queryName{k3},'mag'))(goodidx),...
                wReg(goodidx),'^','Color',colorvec{ie(goodidx(1))},'MarkerFaceColor',colorvec{ie(goodidx(1))}...
                ,'DisplayName',strcat('RegH=',num2str(Hreg(k))));
            if k==1;                                                            % limit r-axis to meaningful frequency range
                hold on; grid on;
                rlim([0 10]);
            end
            
        end
    end
    
    
    % % plots data out-of-range to make legend label generic
    %figure(2*(k3-1)+1)
    %xl=xlim;
    %ax3=plot(-1,-1,'k^','MarkerFaceColor','k');
    %legvec{length(Hunq)+1}='Regular Waves';                                     % make legend labels for future plots
    %xlim(xl);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot BEM cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find body number in BEM output based upon queryName
for k=1:length(queryName)
    if queryName{k}(end)=='1'
        tempmag=flap1mag;
        tempph=flap1ph;
    elseif queryName{k}(end)=='2'                                           % if querying flap 2, loads flap 1 info (BEM run was single-body)
        tempmag=flap1mag;
        tempph=flap1ph;
        disp('Discount Flap 2 phase information. BEM run was single-body')
    else
        tempmag=platmag;
        tempph=platph;
    end
    figure(2*(k-1)+1);                                                      % plot BEM magnitudes
    ax4=plot(wBEM,tempmag(:,k),'k','LineWidth',1.2,'DisplayName','BEM');
    %xlim([2 8])
    
    figure(2*(k-1)+2);                                                      % plot BEM phases
    polarplot(tempph(:,k),wBEM,'sk','LineWidth',1.2,'DisplayName','BEM');
    ax5=gca;
    ax5.RTick=[2 4 6 8 10];
    ax5.RTickLabel={'\omega=2','\omega=4','\omega=6','\omega=8','\omega=10'};
    
    clear tempmag tempph
    
end

figure(2*(k3-3)+1);                                                         % show legend for last mag figure
legend('show')

figure(2*(k3-3)+2);                                                         % show legend for last phase figure
legend('show')






