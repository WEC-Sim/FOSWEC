%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear;
%% User Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
process_inter=0; %1:process *.txt to *.mat, 0: loads *.mat
process_final=1; % 1:process final data structure, 0: loads final *.mat
plot_data=1; % 1 to generate plots

KsPoly=[34.399565513094930,53.132452294804445,0];                          % Nonlinear hydrostatic force resisting flap motion as func of theta (rad)
Iyy= 1.8587;                                                               % flap dry MOI about hinge (for added mass calc (kg m^2)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Directory info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_folder = pwd;
final_folder='./ForcedOscillation';
inter_folder='../inter/ForcedOscillation';
inter_file='ForcedOsc_inter.mat';
final_file='ForcedOsc_final.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process intermediate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if process_inter==1
    % import tst log *.xlsx
    cd(inter_folder)
    [num,txt,raw]=xlsread('WECSIM2_exp2_ForcedOscillation.xlsx','Log');
    data.Exp='ForcedOsc';
    data.Header=txt(18,2:end);
    data.Trial=num(~isnan(num(:,1)),1);                                     % not all number labels, see data.Notes
    data.T=num(:,2);
    data.TauCom=num(:,3);
    data.T_fs=num(:,4);
    data.theta_targ=num(:,8);
    data.theta_obs=num(:,9);
    data.Notes=txt(19:end,8);
    data.Flag=num(:,11);
    
    % corrections (for the NA labels; trial label changes)
    data.theta_targ(10:11)=num([22:23],8);
    data.T(isnan(num(:,1)))=[];
    data.TauCom(isnan(num(:,1)))=[];
    data.T_fs(isnan(num(:,1)))=[];
    data.theta_targ(isnan(num(:,1)))=[];
    data.theta_obs(isnan(num(:,1)))=[];
    data.Notes(isnan(num(:,1)))=[];
    data.Flag(isnan(num(:,1)))=[];
    clear txt num raw
    
    % imports excel data, adds to ForcedOsc file the commanded torque
    % period (T) and amplitude (TauCom), scaled period, target and obs.
    % amplitude of angular positions, notes, and error flag.
    
    % define summary fields
    numTrials=length(data.Trial);
    for k=1:numTrials
        if data.Flag(k)==1;                                                 % if error flag, cell=999 from data struct
            data.T(k)=999;
            data.TauCom(k)=999;
            data.T_fs(k)=999;
            data.theta_targ(k)=999;
            data.theta_obs(k)=999;
        else
            continue
        end
    end
    
    ForcedOsc.T=data.T;
    ForcedOsc.TauCom=data.TauCom;
    ForcedOsc.theta_targ=data.theta_targ;
    ForcedOsc.theta_obs=data.theta_obs;
    ForcedOsc.Flag=data.Flag;
    
    % clear errors
    ForcedOsc.T(ForcedOsc.T==999)=[];
    ForcedOsc.TauCom(ForcedOsc.TauCom==999)=[];
    ForcedOsc.theta_targ(ForcedOsc.theta_targ==999)=[];
    ForcedOsc.theta_obs(ForcedOsc.theta_obs==999)=[];
    
    %% Load data from text files
    for k=1:numTrials
        if ForcedOsc.Flag(k)==0;                                            % does not process flagged cases
            trial_str=sprintf('%02d', data.Trial(k));
            Trial=['Trial' trial_str];
            cd(['./' Trial]);                                               % move to specific trial directory
            contents=dir('*.txt');                                          % find text files
            for j=1:length(contents)
                sensor_file=contents(j).name;
                sensor_name=strrep(sensor_file,'.txt','');
                ForcedOsc.inter.(Trial).(sensor_name)=load(sensor_file);    % load data here
            end
            date_utc=datetime(ForcedOsc.inter.(Trial).time,'convertfrom','datenum');
            time_utc=datetime(date_utc,'TimeZone','UTC');
            datetime_local = datetime(time_utc,'TimeZone','America/Los_Angeles');
            datetime_local.Format = 'dd-MMM-yyyy';
            ForcedOsc.inter.(Trial).Date=datetime_local(1);
            datetime_local.Format = 'HH:mm:ss';
            ForcedOsc.inter.(Trial).TimeStart=datetime_local(1);
            ForcedOsc.inter.(Trial).TimeEnd=datetime_local(end);
            timeStart=ForcedOsc.inter.(Trial).time(1);
            ForcedOsc.inter.(Trial).time=(ForcedOsc.inter.(Trial).time-timeStart)*60*60*24; % convert time to seconds
            cd ..
            clear timeStart
        else
            continue
        end
    end
    save(inter_file,'ForcedOsc');
    clearvars -except ForcedOsc final_folder final_file inter_file inter_folder process_final plot_data base_folder
    cd(base_folder)
else
    load([inter_folder,'\',inter_file]);                                    % load 'inter' *.mat
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process Final Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if process_final==1;
    if exist(final_folder) ==  0
        mkdir(final_folder)
    end
    cd(final_folder)
    addpath('../../../functions')
    assymThresh=0.1;                                                        % allowable assymetry in oscillation; given as fraction of target amplitude
    assymWin=5;                                                             % number of adjacent peaks to average to consider assymetry
    %% check consistency of sampling frequency
    Fnames=fieldnames(ForcedOsc.inter);
    numFields=length(Fnames);
    dtFlag=zeros(numFields,1);
    for k=1:numFields                                                       % flags >1% sampling time error
        dterr=max(abs(diff(ForcedOsc.inter.(Fnames{k}).time)-0.01));
        if dterr > 10E-4
            dtFlag(k)=1;
        end
        dt(k)=mean(diff(ForcedOsc.inter.(Fnames{k}).time));
    end
    
    %% Find flap velocity, acceleration, torque signals. Trim data.
    filtNum=0.125*ones(8,1);                                                % FOSWEC sensor filter coefficients
    
    % pre-allocation
    startIdx=zeros(numFields,1);
    peakLoc{numFields,1}=zeros(1,1);
    peakmag{numFields,1}=zeros(1,1);
    minLoc{numFields,1}=zeros(1,1);
    minmag{numFields,1}=zeros(1,1);
    goodidx{numFields,1}=zeros(1,1);
    contFlag=zeros(numFields,1);
    
    for k=1:numFields;
        if dtFlag(k)==0;                                                    % check constant sampling time
            FiltPos=filtFB(filtNum,1,ForcedOsc.inter.(Fnames{k}).flapPosF1,[],2).*pi/180; % convert to radians
            ForcedOsc.final.TimeDomain.(Fnames{k}).flapVelF1=diff(FiltPos)./dt(k);
            ForcedOsc.final.TimeDomain.(Fnames{k}).flapAccF1=filtFB(filtNum,1, ...
                diff(ForcedOsc.final.TimeDomain.(Fnames{k}).flapVelF1)./dt(k),[],2);
            ForcedOsc.final.TimeDomain.(Fnames{k}).flapT1=filtFB(filtNum,1, ...
                ForcedOsc.inter.(Fnames{k}).lcTyF1,[],2);
        end
        
        % find the start of oscillations
        startIdx(k)=findStartIdx(ForcedOsc.inter.(Fnames{k}).time,...           % Trim by start of flap forcing.
            ForcedOsc.inter.(Fnames{k}).flapPosF1, 2, -4, ForcedOsc.T(k),4);    % Negative leadNum argument ensure full amp osc.
        ForcedOsc.final.TimeDomain.(Fnames{k}).flapVelF1=ForcedOsc.final.TimeDomain.(Fnames{k}).flapVelF1(startIdx(k):end);
        ForcedOsc.final.TimeDomain.(Fnames{k}).flapAccF1=ForcedOsc.final.TimeDomain.(Fnames{k}).flapAccF1(startIdx(k):end);
        ForcedOsc.final.TimeDomain.(Fnames{k}).flapT1=ForcedOsc.final.TimeDomain.(Fnames{k}).flapT1(startIdx(k):end);
        ForcedOsc.final.TimeDomain.(Fnames{k}).flapPosF1=FiltPos(startIdx(k):end);
        
        %% Detect and remove drift from neutral mean position
        % theta from the resulting calculations. Uses matlab file exchange
        % contribution 'peakfinder'.
        
        [peakLoc{k,1},peakmag{k,1}]=peakfinder(ForcedOsc.final.TimeDomain.(Fnames{k}).flapPosF1,...
            [],-10,1,false,false);
        [minLoc{k,1},minmag{k,1}]=peakfinder(ForcedOsc.final.TimeDomain.(Fnames{k}).flapPosF1,...
            [],10,-1,false,false);
        
        % detect points within the assymetry threshold
        assymThreshActual=assymThresh*ForcedOsc.theta_targ(k)*pi/180;
        goodidx{k,1}=parseAssymetry(peakmag{k,1},peakLoc{k,1},minmag{k,1},minLoc{k,1},...
            assymThreshActual, assymWin, length(ForcedOsc.final.TimeDomain.(Fnames{k}).flapPosF1));
        
        if max(diff(goodidx{k,1})) > 1;                                    % flag for discontinuous retained time series
            contFlag(k)=1;
        end
        
        % restrict consideration to good points
        ForcedOsc.final.TimeDomain.(Fnames{k}).flapVelF1=ForcedOsc.final.TimeDomain.(Fnames{k}).flapVelF1(goodidx{k,1});
        ForcedOsc.final.TimeDomain.(Fnames{k}).flapAccF1=ForcedOsc.final.TimeDomain.(Fnames{k}).flapAccF1(goodidx{k,1});
        ForcedOsc.final.TimeDomain.(Fnames{k}).flapT1=ForcedOsc.final.TimeDomain.(Fnames{k}).flapT1(goodidx{k,1}).*-1; % correct flipped load cell
        ForcedOsc.final.TimeDomain.(Fnames{k}).flapPosF1=ForcedOsc.final.TimeDomain.(Fnames{k}).flapPosF1(goodidx{k,1});
        
        % find peaks of retained points for windowing analysis
        [peakLoc{k,1},peakmag{k,1}]=peakfinder(ForcedOsc.final.TimeDomain.(Fnames{k}).flapPosF1,...
            [],-10,1,true,false);
        
        minL=min(floor(diff(peakLoc{k,1})/2));
        
        
        %% Analyze parameters in time domain
        
        if ~isempty(peakLoc{k,1})
            for k2=1:length(peakLoc{k,1})-1
                
                % window for velocity; find min acceleration each half-period
                window=peakLoc{k,1}(k2):peakLoc{k,1}(k2+1);
                winL=length(window);
                
                %for each window, define matrix x,b to solve for A in Ax=b.
                Theta=ForcedOsc.final.TimeDomain.(Fnames{k}).flapPosF1(window);
                HydroStat=sign(Theta).*polyval(KsPoly,abs(Theta));
                StateMat(:,2)=ForcedOsc.final.TimeDomain.(Fnames{k}).flapVelF1(window);
                StateMat(:,1)=ForcedOsc.final.TimeDomain.(Fnames{k}).flapAccF1(window);
                TauTotMat(:,1)=(ForcedOsc.final.TimeDomain.(Fnames{k}).flapT1(window)...
                    -HydroStat);
                StateMat(:,3)= ForcedOsc.final.TimeDomain.(Fnames{k}).flapVelF1(window)...
                    .*abs(ForcedOsc.final.TimeDomain.(Fnames{k}).flapVelF1(window)); % for quadratic damping
                
                %% Calculation assuming Cv or Cd damping sources, solve Ax=b
                CvOnly{k,1}(:,k2)=StateMat(:,[1,2])\TauTotMat;
                CdOnly{k,1}(:,k2)=StateMat(:,[1,3])\TauTotMat;
                
                % calculate residuals
                normErrCv{k,1}(k2,1)=norm(StateMat(:,[1,2])*CvOnly{k,1}(:,k2)-TauTotMat)./length(TauTotMat);
                normErrCd{k,1}(k2,1)=norm(StateMat(:,[1,3])*CdOnly{k,1}(:,k2)-TauTotMat)./length(TauTotMat);
                
                
                clear TauTotMat StateMat
            end
            
            % Asign to final data structure
            ForcedOsc.final.AT(k,1)=mean(CvOnly{k}(1,:));
            ForcedOsc.final.Astd(k,1)=std(CvOnly{k}(1,:));
            ForcedOsc.final.CvTonly(k,1)=mean(CvOnly{k}(2,:));
            ForcedOsc.final.CvTonlystd(k,1)=std(CvOnly{k}(2,:));
            ForcedOsc.final.CdTonly(k,1)=mean(CdOnly{k}(2,:));
            ForcedOsc.final.CdTonlystd(k,1)=std(CdOnly{k}(2,:));
            
            
        end
    end
    %% Calculate Cd based off of energy decrement linearization from Cv estimate.
    % for each oscillation frequency
    
    % Method of Inman (2008), linearization of quadratic damping system.
    
    [Tunique,ia,ic]=unique(ForcedOsc.T);
    for k=1:length(Tunique)
        goodidx=find(ic==k);
        
        fitVal=8.*(2*pi./ForcedOsc.T(goodidx)).*(ForcedOsc.theta_obs(goodidx).*(pi/180))./(3*pi);
        figure;clf
        plot(fitVal,ForcedOsc.final.CvTonly(goodidx),'ob');
        hold on; grid on;
        dampFit=polyfit(fitVal,ForcedOsc.final.CvTonly(goodidx),1);
        plot(fitVal,polyval(dampFit,fitVal),'--k','LineWidth',1.2);
        xlabel('8\omega_d \theta/3 \pi')
        ylabel('c_{tot}')
        legend('Experimental Data','Fit','Location','SouthEast')
        
        % cquad is the quadratic damping component; clin is the total linear
        % damping component (wave radiation damping and linear mechanical damping)
        % for this frequency
        ForcedOsc.final.cquad(k)=mean(dampFit(1)./fitVal);
        ForcedOsc.final.clin(k)=dampFit(2);
        ForcedOsc.final.Amass(k)=median(ForcedOsc.final.AT(goodidx,1));
        ForcedOsc.final.AmassStd(k)=std(ForcedOsc.final.AT(goodidx,1));
    end
    save(final_file,'ForcedOsc')
else
    cd(final_folder)
    load(final_file);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_data==1;
    %% Plot of frequency vs Cv, Cd, Added mass for independent determination.
    
    close all;
    figure(1);clf;
    figure(2);clf;
    figure(3);clf;
    [Tunique,ia,ic]=unique(ForcedOsc.T);
    for k=1:length(Tunique)
        goodidx=find(ic==k);
        % offset plotted value along x axis to aid visualization for
        % distinct amplitudes
        off10=0;                                                            % offset for 10 dg oscillations
        off15=0;                                                            % offset for 15 dg oscillations
        off20=0;                                                            % offset for 20 dg oscillations
        for k2=1:length(goodidx)
            if ForcedOsc.final.CvTonly(goodidx(k2)) ~= 0;                   % uncalculated cases were left to zero
                switch ForcedOsc.theta_targ(goodidx(k2))
                    case 10
                        off=off10;
                        cString='b';
                        mType='o';
                        
                    case 15
                        off=off15;
                        cString='r';
                        mType='^';
                        
                    case 20
                        off=off20;
                        cString='g';
                        mType='s';
                        
                end
                %% Using Cv Calculation
                figure(1)
                subplot(2,1,2)
                scatter(2*pi/(ForcedOsc.T(goodidx(k2)))+off...              % plot data
                    ,ForcedOsc.final.AT(goodidx(k2))-Iyy,mType,'MarkerFaceColor',cString)
                if k2==1;
                    hold on;
                    grid on;
                    xlabel('\omega (rad/s)')
                    ylabel('A (kg m^2)')
                end
                line([2*pi/(ForcedOsc.T(goodidx(k2)))+off, 2*pi/(ForcedOsc.T(goodidx(k2)))+off],... % error bars (+/- 2 std)
                    [2*ForcedOsc.final.Astd(goodidx(k2))+ForcedOsc.final.AT(goodidx(k2)), ...
                    ForcedOsc.final.AT(goodidx(k2))-2*ForcedOsc.final.Astd(goodidx(k2))], 'LineStyle','-',...
                    'Color',cString,'LineWidth',1.2);
                
                subplot(2,1,1)
                switch ForcedOsc.theta_targ(goodidx(k2))                    % format data by amplitude, plot
                    case 10
                        ax10=scatter(2*pi/(ForcedOsc.T(goodidx(k2)))+off...
                            ,ForcedOsc.final.CvTonly(goodidx(k2)),mType,'MarkerFaceColor',cString);
                    case 15
                        ax15=scatter(2*pi/(ForcedOsc.T(goodidx(k2)))+off...
                            ,ForcedOsc.final.CvTonly(goodidx(k2)),mType,'MarkerFaceColor',cString);
                    case 20
                        ax20=scatter(2*pi/(ForcedOsc.T(goodidx(k2)))+off...
                            ,ForcedOsc.final.CvTonly(goodidx(k2)),mType,'MarkerFaceColor',cString);
                end
                
                if k2==1;
                    hold on;
                    grid on;
                    xlabel('\omega (rad/s)')
                    ylabel('C_v')
                end
                line([2*pi/(ForcedOsc.T(goodidx(k2)))+off, 2*pi/(ForcedOsc.T(goodidx(k2)))+off],... % error bar (+/- 2 std)
                    [2*ForcedOsc.final.CvTonlystd(goodidx(k2))+ForcedOsc.final.CvTonly(goodidx(k2)), ...
                    -2*ForcedOsc.final.CvTonlystd(goodidx(k2))+ForcedOsc.final.CvTonly(goodidx(k2))], 'LineStyle','-',...
                    'Color',cString,'LineWidth',1.2);
                
                %% Using Cd Calculation
                figure(2)
                subplot(2,1,2)
                scatter(2*pi/(ForcedOsc.T(goodidx(k2)))+off...              % plot data
                    ,ForcedOsc.final.AT(goodidx(k2))-Iyy,mType,'MarkerFaceColor',cString)
                if k2==1;
                    hold on;
                    grid on;
                    xlabel('\omega (rad/s)')
                    ylabel('A (kg m^2)')
                end
                line([2*pi/(ForcedOsc.T(goodidx(k2)))+off, 2*pi/(ForcedOsc.T(goodidx(k2)))+off],... % error bars (+/- 2 std)
                    [2*ForcedOsc.final.Astd(goodidx(k2))+ForcedOsc.final.AT(goodidx(k2)), ...
                    ForcedOsc.final.AT(goodidx(k2))-2*ForcedOsc.final.Astd(goodidx(k2))], 'LineStyle','-',...
                    'Color',cString,'LineWidth',1.2);
                
                subplot(2,1,1)
                switch ForcedOsc.theta_targ(goodidx(k2))                    % format data by amplitude, plot
                    case 10
                        ax101=scatter(2*pi/(ForcedOsc.T(goodidx(k2)))+off...
                            ,ForcedOsc.final.CdTonly(goodidx(k2)),mType,'MarkerFaceColor',cString);
                    case 15
                        ax151=scatter(2*pi/(ForcedOsc.T(goodidx(k2)))+off...
                            ,ForcedOsc.final.CdTonly(goodidx(k2)),mType,'MarkerFaceColor',cString);
                    case 20
                        ax201=scatter(2*pi/(ForcedOsc.T(goodidx(k2)))+off...
                            ,ForcedOsc.final.CdTonly(goodidx(k2)),mType,'MarkerFaceColor',cString);
                end
                
                if k2==1;
                    hold on;
                    grid on;
                    xlabel('\omega (rad/s)')
                    ylabel('C_d')
                end
                line([2*pi/(ForcedOsc.T(goodidx(k2)))+off, 2*pi/(ForcedOsc.T(goodidx(k2)))+off],... % error bar (+/- 2 std)
                    [2*ForcedOsc.final.CdTonlystd(goodidx(k2))+ForcedOsc.final.CdTonly(goodidx(k2)), ...
                    -2*ForcedOsc.final.CdTonlystd(goodidx(k2))+ForcedOsc.final.CdTonly(goodidx(k2))], 'LineStyle','-',...
                    'Color',cString,'LineWidth',1.2);
                
            else
                continue
            end
        end
        
        
    end
    figure(1)
    subplot(2,1,1)
    legend([ax10 ax15 ax20], {'\Delta \Theta=10^o','\Delta \Theta=15^o',...
        '\Delta \Theta=20^o'},'Location','SouthEast');
    
    figure(2)
    subplot(2,1,1)
    legend([ax101 ax151 ax201], {'\Delta \Theta=10^o','\Delta \Theta=15^o',...
        '\Delta \Theta=20^o'},'Location','NorthEast');
    
    %% Linear approximation plot of linear and quadratic damping ,combined calculation
    
    figure(3)
    plot((2*pi)./Tunique,ForcedOsc.final.clin,'ob');
    hold on; grid on;
    plot((2*pi)./Tunique,ForcedOsc.final.cquad,'^r');
    xlabel('Freq. (rad/s)')
    ylabel('Damping')
    legend('c_{lin}','c_{quad}','Location','NorthWest')
end
