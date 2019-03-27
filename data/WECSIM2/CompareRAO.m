% Compare Response Amplitude Operators for regular + irregular,
% and device configuration number, facilitating comparison to simulation.

% User inputs specify regular wave, irregular wave, or both plotting, the
% device configuration number of interest, the field names of the data
% that is of interest, and the corresponding plot titles.


clear; close all; clc;

% Directory Information
homedir=pwd;
basedir='./final';

% User Inputs
plotWhat=[1 2]; % 1 to plot regular wave cases, 2 to plot irregular wave cases, [1:2] or [1,2] to plot both
configNum=4; % a single selection of [1:4], specifies configuration on which to plot simulation data. Just use one
% 1: only flap 1 unlocked. 2: flaps 1 and 2 unlocked. 3: flaps and platform
% heave unlocked. 4: flaps and platform heave, surge and pitch are
% unlocked.
queryName={'flapPosF1','flapPosF2','platPosz','platPosx','platPosRy'};      % a field name specifying the phase/mag plot to show, if present
titStr={'RAO Flap 1','RAO Flap 2','RAO Plat Z','RAO Plat X','RAO Plat Ry'}; % magnitude plot title strings that correspond to the above
ptitStr={'\angle RAO Flap 1','\angle RAO Flap 2','\angle RAO Plat Z',...    % phase plot title strings
    '\angle RAO Plat X','\angle RAO Plat Ry'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load experiment data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(basedir)

% load each configuration
RegStr=strcat('./Config',num2str(configNum),'Reg');
RegFile=strcat('Config',num2str(configNum),'Reg_final.mat');
IrrStr=strcat('./Config',num2str(configNum),'Irr');
IrrFile=strcat('Config',num2str(configNum),'Irr_final.mat');
cd(RegStr)
load(RegFile)                                                               % load regular wave runs
cd ..
cd(IrrStr)
load(IrrFile)                                                               % load irregular wave runs
cd ..
cd(homedir)
cd(basedir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step through query fields by Configuration number.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch configNum                                                            % selects top-level field name based on device configurations
    case 1
        regFile=Config1Reg;
        irrFile=Config1Irr;
        cVal={[0 0.447 0.7410],[0.85 0.325 0.0980],[0.929 0.694 0.125],'k'}; % color encodes wave height. 4 wave heights for Config. 1, 3 for others.
        cValIrr={[0.85 0.325 0.0980],[0.929 0.694 0.125],'k'};
    case 2
        regFile=Config2Reg;
        irrFile=Config2Irr;
        cVal={[0.85 0.325 0.0980],[0.929 0.694 0.125],'k'};
        cValIrr={[0.85 0.325 0.0980],[0.929 0.694 0.125],'k'};
    case 3
        regFile=Config3Reg;
        irrFile=Config3Irr;
        cVal={[0.85 0.325 0.0980],[0.929 0.694 0.125],'k'};
        cValIrr={[0.85 0.325 0.0980],[0.929 0.694 0.125],'k'};
    case 4
        regFile=Config4Reg;
        irrFile=Config4Irr;
        cVal={[0.85 0.325 0.0980],[0.929 0.694 0.125],'k'};
        cValIrr={[0.85 0.325 0.0980],[0.929 0.694 0.125],'k'};
end

for k2=1:length(queryName)                                                  % query name loop
    if isfield(regFile.final,queryName{k2})==1                              % if the specified query name exists for a device config., it will be plotted
        makePlot=1;
        figure(3*(k2-1)+1);
    else;
        makePlot=0;
    end
    if makePlot==1;                                                         % plots specified field (if exist)
        [Hreg,ia,ic]=unique(regFile.H);
        [Hirr,id,ie]=unique(irrFile.H);
        switch k2                                                           % selects axis label with units based on mode of motion specified by query vec.
            case 1
                axisLabel='RAO deg/m';
            case 2
                axisLabel='RAO deg/m';
            case 3
                axisLabel='RAO m/m';
            case 4
                axisLabel='RAO m/m';
            case 5
                axisLabel='RAO deg/m';
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Plot Regular wave cases
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if min(plotWhat)==1;
            for k3=1:length(Hreg);                                              % regular wave height loop
                goodidx=find(ic==k3);
                formatString={'MarkerFaceColor',cVal{k3},'MarkerEdgeColor',cVal{k3},'marker','s'};  % color code based on wave height
                formatStringSpec={'MarkerFaceColor',cVal{k3},'MarkerEdgeColor',cVal{k3},'marker','^'};
                formatStringPolar={'MarkerFaceColor',cVal{k3},'MarkerEdgeColor',cVal{k3},'marker','^'};
                figure(3*(k2-1)+1);
                ax(k3)=scatter(2*pi./regFile.T(goodidx),regFile.final.(queryName{k2}).RAOmean(goodidx)... % Plot RAO mag mean from time series calculation
                    ,36,formatString{:});
                if k3==1;                                                       % add gridlines, label axis
                    hold on
                    grid on
                    title(titStr{k2});
                    xlim([0 12])
                    xlabel('Frequency (rad/s)')
                    ylabel(axisLabel)
                end
                % make error bars
                for k5=1:length(2*pi./regFile.T(goodidx))
                    figure(3*(k2-1)+1);
                    line([2*pi./regFile.T(goodidx(k5)), 2*pi./regFile.T(goodidx(k5))],[(regFile.final.(queryName{k2}).RAOmean(goodidx(k5))... % error bars show +/- 2 sigma
                        -2.*regFile.final.(queryName{k2}).RAOstd(goodidx(k5))), (regFile.final.(queryName{k2}).RAOmean(goodidx(k5))...
                        +2.*regFile.final.(queryName{k2}).RAOstd(goodidx(k5)))],'LineStyle','-','LineWidth',1.2,'Color',cVal{k3});
                end
                
                % spectral estimate
                ax2(k3)=scatter(2*pi./regFile.T(goodidx),regFile.final.(queryName{k2}).RAOmag(goodidx)... % plot RAO mag mean from spectral method
                    ,36,formatStringSpec{:});
                
                % phase plot
                figure(3*(k2-1)+2);
                ax4(k3)=polarscatter(regFile.final.(queryName{k2}).RAOph(goodidx),2*pi./regFile.T(goodidx)...
                    ,36,formatStringPolar{:});
                if k3==1
                    hold on;
                    grid on;
                    temp=gca;
                    temp.RLim=[0 12];
                    temp.RTick=[2:2:12];
                    %temp.RTickLabel=[2:2:12);
                    title(ptitStr{k2})
                end
                
            end
            figure(3*(k2-1)+1)
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Plot Irregular wave cases
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if max(plotWhat)==2
            for k4=1:length(Hirr);
                goodidx=find(ie==k4);
                % format string based upon the significant wave height.
                formatStringIrr={'Color',cValIrr{k4},'LineWidth',1.4};
                formatStringPolar={'marker','o','MarkerFaceColor',cValIrr{k4},'MarkerEdgeColor',cValIrr{k4}};
                if k4==1;                                                       % add gridlines, label axis
                    hold on
                    grid on
                    title(titStr{k2});
                    xlim([0 12])
                    xlabel('Frequency (rad/s)')
                    ylabel(axisLabel)
                end
                % magnitude plot
                figure(3*(k2-1)+1);
                ax3(k4)=plot(irrFile.final.freqAvgw{1,k4},irrFile.final.(queryName{k2}).RAOMagAvg{k4},formatStringIrr{:}); % plot RAO mag from spectral method (no time-domain method for irregular wave)
                
                % phase plot
                figure(3*(k2-1)+2);
                polarscatter(irrFile.final.(queryName{k2}).RAOphAvg{k4},irrFile.final.freqAvgw{1,k4},...
                    36,formatStringPolar{:});
                if k4==1
                    hold on;
                    grid on;
                    temp=gca;
                    temp.RLim=[0 12];
                    temp.RTick=[2:2:12];
                    %temp.RTickLabel=[2:2:12);
                    title(ptitStr{k2})
                end
                
            end
        else
            legStr3=[];                                                         % makes empty vector label in case irregular wave cases are not to be plotted.
        end
    else
        continue
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make plot legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the order and marker style here corresponds to color vector and the
% plotting order (regular time-domain, regular frequency domain, then irregular).

if min(plotWhat)==1;
    if length(cVal)==4                                                         % when 4 unique wave heights for regular waves
        figure(1);                                                          % magnitude plot
        xl=xlim;
        h1(1)=scatter(NaN,NaN,36,'MarkerFaceColor',cVal{1},'MarkerEdgeColor',cVal{1},'marker','s',... % time-domain estimates, reg waves
            'DisplayName',strcat('Hreg=',num2str(Hreg(1)))); hold on;
        h1(2)=scatter(NaN,NaN,36,'MarkerFaceColor',cVal{2},'MarkerEdgeColor',cVal{2},'marker','s',...
            'DisplayName',strcat('Hreg=',num2str(Hreg(2))));
        h1(3)=scatter(NaN,NaN,36,'MarkerFaceColor',cVal{3},'MarkerEdgeColor',cVal{3},'marker','s',...
            'DisplayName',strcat('Hreg=',num2str(Hreg(3))));
        h1(4)=scatter(NaN,NaN,36,'MarkerFaceColor',cVal{1},'MarkerEdgeColor',cVal{1},'marker','^',... % spectral estimates, reg waves
            'DisplayName',strcat('Spec Hreg=',num2str(Hreg(1))));
        h1(5)=scatter(NaN,NaN,36,'MarkerFaceColor',cVal{2},'MarkerEdgeColor',cVal{2},'marker','^',...
            'DisplayName',strcat('Spec Hreg=',num2str(Hreg(2))));
        h1(6)=scatter(NaN,NaN,36,'MarkerFaceColor',cVal{3},'MarkerEdgeColor',cVal{3},'marker','^',...
            'DisplayName',strcat('Spec, Hreg=',num2str(Hreg(3))));
        xlim(xl)
        
        figure(2); xl2=rlim;                                                % phase plot
        h2(1)=polarscatter(NaN,NaN,36,'MarkerFaceColor',cVal{1},'MarkerEdgeColor',cVal{1},'marker','^',...
            'DisplayName',strcat('Spec Hreg=',num2str(Hreg(1))));
        h2(2)=polarscatter(NaN,NaN,36,'MarkerFaceColor',cVal{2},'MarkerEdgeColor',cVal{2},'marker','^',...
            'DisplayName',strcat('Spec Hreg=',num2str(Hreg(2))));
        h2(3)=polarscatter(NaN,NaN,36,'MarkerFaceColor',cVal{3},'MarkerEdgeColor',cVal{3},'marker','^',...
            'DisplayName',strcat('Spec Hreg=',num2str(Hreg(3))));
        rlim(xl2)
        axIdx=6;                                                            % mag. plot axis definition index
        axIdx2=3;                                                           % phase plot axis definition index
        
    elseif length(cVal)==3
        figure(1); xl=xlim;                                                 % magnitude plot
        h1(1)=scatter(NaN,NaN,36,'MarkerFaceColor',cVal{1},'MarkerEdgeColor',cVal{1},'marker','s',... % time-domain estimates, reg waves
            'DisplayName',strcat('Hreg=',num2str(Hreg(1)))); hold on;
        h1(2)=scatter(NaN,NaN,36,'MarkerFaceColor',cVal{2},'MarkerEdgeColor',cVal{2},'marker','s',...
            'DisplayName',strcat('Hreg=',num2str(Hreg(2))));
        h1(3)=scatter(NaN,NaN,36,'MarkerFaceColor',cVal{1},'MarkerEdgeColor',cVal{1},'marker','^',... % spectral estimates, reg waves
            'DisplayName',strcat('Spec Hreg=',num2str(Hreg(1))));
        h1(4)=scatter(NaN,NaN,36,'MarkerFaceColor',cVal{2},'MarkerEdgeColor',cVal{2},'marker','^',...
            'DisplayName',strcat('Spec Hreg=',num2str(Hreg(2))));
        xlim(xl)
        
        figure(2); xl2=rlim;                                                % phase plot
        h2(1)=polarscatter(NaN,NaN,36,'MarkerFaceColor',cVal{1},'MarkerEdgeColor',cVal{1},'marker','^',...
            'DisplayName',strcat('Spec Hreg=',num2str(Hreg(1))));
        h2(2)=polarscatter(NaN,NaN,36,'MarkerFaceColor',cVal{2},'MarkerEdgeColor',cVal{2},'marker','^',...
            'DisplayName',strcat('Spec Hreg=',num2str(Hreg(2))));
        rlim(xl2);
        
        axIdx=4;                                                            % mag. plot axis definition index
        axIdx2=2;                                                           % phase plot axis definition index
    end
else
    axIdx=0;
    axIdx2=0;
    
end
if max(plotWhat)==2;                                                        % spectral estimates, irregular waves
    figure(1); xl=xlim;                                                     % mag. plot
    h1(axIdx+1)=plot(NaN,NaN,'-','LineWidth',1.4,'Color',cValIrr{1},...
        'DisplayName',strcat('Hirr=',num2str(Hirr(1))));
    h1(axIdx+2)=plot(NaN,NaN,'-','LineWidth',1.4,'Color',cValIrr{2},...
        'DisplayName',strcat('Hirr=',num2str(Hirr(2))));
    xlim(xl)
    
    figure(2); xl2=rlim;                                                   % phase plot
    h2(axIdx2+1)=polarscatter(NaN,NaN,36,'MarkerFaceColor',cValIrr{1},'MarkerEdgeColor',cValIrr{1},'marker','o',...
        'DisplayName',strcat('Hirr=',num2str(Hirr(1))));
    h2(axIdx2+2)=polarscatter(NaN,NaN,36,'MarkerFaceColor',cValIrr{2},'MarkerEdgeColor',cValIrr{2},'marker','o',...
        'DisplayName',strcat('Hirr=',num2str(Hirr(2))));
    rlim(xl2);
    
end
figure(1);
legend(h1)
xlim(xl);

figure(2);
legend(h2)
rlim(xl2)
