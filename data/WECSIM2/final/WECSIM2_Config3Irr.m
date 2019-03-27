% Calculates Response Amplitude operators for Configuration 3 irregular
% waves. Flaps are unlocked, platform is free in heave.

%Uses phase space data when available. When not, throws warning and uses
% FOSWEC mounted sensors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
process_inter = 0;   % 1:processes inter data, 0:loads inter *.mat
process_final = 0;   % 1:processes final data, 0:loads final *.mat
plot_data = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Directory Info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_folder = pwd;
final_folder = './Config3Irr';
inter_folder = '../inter/Config3Irr';
inter_file = 'Config3Irr_inter.mat';
final_file = 'Config3Irr_final.mat';
addpath(genpath(strrep(pwd,'\WECSIM2\final\Config3Irr','')))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import Test Log and 'inter' Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if process_inter == 1
    % import test log *.xlsx
    cd(inter_folder);
    [num,txt,raw]   = xlsread('WECSIM2_Config3Irr.xlsx','Log');
    data.Exp        = 'Config3Irr';
    data.Header     = txt(5,2:end);
    data.Trial      = num(:,2);
    data.T          = num(:,3);
    data.H          = num(:,4);
    data.T_fs       = num(:,5);
    data.H_fs       = num(:,6);
    data.Notes      = txt(6:end,9);
    data.Flag       = num(:,10);
    clear txt num raw
    
    Config3Irr.T          = data.T;
    Config3Irr.H          = data.H;
    Config3Irr.Flag       = data.Flag;
    
    % import 'inter' data for all trials
    numTrials = length(data.Trial);
    for i = 1:numTrials
        trial_str = sprintf('%02d',data.Trial(i));
        Trial = ['Trial' trial_str];
        cd(['./' Trial])
        Config3Irr.inter.(Trial).T          = data.T(i);
        Config3Irr.inter.(Trial).H          = data.H(i);
        Config3Irr.inter.(Trial).Notes      = data.Notes(i);
        Config3Irr.inter.(Trial).Flag       = data.Flag(i);
        Config3Irr.inter.(Trial).T_fs       = data.T_fs(i);
        Config3Irr.inter.(Trial).H_fs       = data.H_fs(i);
        contents = dir( '*.txt' );
        for j = 1:length(contents)           % import all sensor data
            sensor_file = contents(j).name;
            sensor_name = strrep(sensor_file,'.txt','');
            Config3Irr.inter.(Trial).(sensor_name) = load(sensor_file);     % load data here
        end
        date_utc = datetime(Config3Irr.inter.(Trial).time,'convertfrom','datenum');
        time_utc = datetime(date_utc,'TimeZone','UTC');
        datetime_local = datetime(time_utc,'TimeZone','America/Los_Angeles');
        datetime_local.Format = 'dd-MMM-yyyy';
        Config3Irr.inter.(Trial).Date = datetime_local(1);
        datetime_local.Format = 'HH:mm:ss';
        Config3Irr.inter.(Trial).TimeStart = datetime_local(1);
        Config3Irr.inter.(Trial).TimeEnd = datetime_local(end);
        timeStart = Config3Irr.inter.(Trial).time(1);
        Config3Irr.inter.(Trial).time = (Config3Irr.inter.(Trial).time - timeStart)*60*60*24;
        cd ..
    end
    save(inter_file,'Config3Irr')                                          % save 'inter' data to *.mat
    
    clearvars -except Config3Irr final_folder final_file  inter_file inter_folder process_final plot_data base_folder
    cd(base_folder)
else
    load([inter_folder,'/', inter_file])                                   % load 'inter' *.mat
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if process_final == 1
    if exist(final_folder) ==  0
        mkdir(final_folder)
    end
    cd(final_folder)
    addpath('../../../functions')
    dur=20;                                                                 % number of peak wave periods to consider
    numTrials=length(fieldnames(Config3Irr.inter));
    dtFlag=zeros(numTrials,1);
    filtNum=0.125*ones(8,1);                                                % time domain filter coefficients for 100 Hz data
    filtNumW=0.25*ones(4,1);                                                % time domain filter coefficients for 50 Hz data
    thresh=0.95;                                                            % Above this percentile (as a decimal) will be retained
    [Hunq,ia,ic]=unique(Config3Irr.H);
    n=3;                                                                    % number of adjacent frequency bands to merge 
    smoothNum=ones(n,1).*(1/n);   
    % filter entire inter dataset
    i_min = [3000; 3000; 3000; 3000];
    
    % load undisturbed wave time series
    cd('../../../WECSIM/final');
    load('WaveTuning');
    cd('../../WECSIM2/final/Config3Irr');
    
    for i = 1:numTrials
        trial_str = sprintf('%02d',i);
        Trial = ['Trial' trial_str];
        
        Fs = 1./mean(diff(Config3Irr.inter.(Trial).time));                 % sampling frequency
        dt(i) = 1./Fs;
        dterr=max(abs(diff(Config3Irr.inter.(Trial).time)-dt(i)));         % sampling period
        
        % check consistency of sampling interval, flag if inconsistent
        if dterr>10E-4
            dtFlag(i)=1;
        end
        
        if dtFlag<1;
            
            % filter entire inter dataset; uses phase-space data for
            % platform heave. 
            time = Config3Irr.inter.(Trial).time;
            flapPosF1 = filtFB(filtNum,1,Config3Irr.inter.(Trial).flapPosF1,[],2);
            flapPosF2 = filtFB(filtNum,1,Config3Irr.inter.(Trial).flapPosF2,[],2); 
            platPosz = filtFB(filtNumW,1,Config3Irr.inter.(Trial).heave,[],2); % interpolate to faster time. this may shift by up to 0.01 seconds!
            tidx=0:0.02:0.02*(length(platPosz)-1);
            platPosz=interp1(tidx,platPosz,Config3Irr.inter.(Trial).time);
            
            %% import wave gauge data from calibration runs
            Tq=Config3Irr.T(i);                                             % query wave period
            Hq=Config3Irr.H(i);                                             % query wave height
            goodidx=find(abs(Tq-WaveTuning.RandWave.T)<0.05);               % find calibrated wave runs matching T,H.
            goodidx2=find(abs(Hq-WaveTuning.RandWave.H)<0.01);
            goodidx=intersect(goodidx,goodidx2);
            if length(goodidx)>1                                            % if multiple matching runs, take the first
                goodidx=goodidx(1);
            end
            
            % wave elevations
            fname=strcat('Test',num2str(goodidx));
            TS=WaveTuning.RandWave.Tests.(fname).time;
            dtW(i)=mean(diff(TS));                                          % sampling interval for waves
            Config3Irr.final.WaveTS{i,1}=filtFB(filtNumW,1,detrend(WaveTuning.RandWave.Tests.(fname).wg6,'constant'),[],2);
            wg6= Config3Irr.final.WaveTS{i,1};
            Tn=Tq./dt(i);                                                  % number of encoder samples per wave period
            Tw=Tq./dtW(i);                                                 % number of wave gauge samples per wave period
            
            % trim start time, cut off end based upon specified wave duration
            flapPosF1 = flapPosF1(i_min(i):i_min(i)+dur*Tn+1);              % known a priori by observation
            flapPosF2 = flapPosF2(i_min(i):i_min(i)+dur*Tn+1);
            platPosz= platPosz(i_min(i):i_min(i)+dur*Tn+1);
            [~,Widx(i)]=min(abs(TS-time(i_min(i))));                        % idx of corresponding start time in wave TS
            wg6 = wg6(Widx(i):Widx(i)+dur*Tw);                              % calculated by differing sample times
            
            % interpolate to wave times
            flapPosF1=interp1(time(1:dur.*Tn+2),flapPosF1,TS(1:dur.*Tw+1),'linear','extrap');
            flapPosF2=interp1(time(1:dur.*Tn+2),flapPosF2,TS(1:dur.*Tw+1),'linear','extrap');
            platPosz=interp1(time(1:dur.*Tn+2),platPosz,TS(1:dur.*Tw+1),'linear','extrap');
            
            %% Fourier transform procedure
            L = length(flapPosF1);                                          % length of signal
            L2=2^nextpow2(L);
            Y=fft(detrend(flapPosF1,'linear'),L2)./L2;                      % fourier transform
            f = (1/dtW(i))*(0:(L2/2))/L2;                                   % frequency domain
            PP = 2*pi.*f;                                                   % Period (rad/s)
            % mean_y = mean(flapPosF1);
            temp=Y(1:L2/2+1);
            temp(2:end-1)=2*temp(2:end-1);
            Y=temp;
            Ymag = abs(Y);                                                  % magnitude
            Yph = atan2(imag(Y),real(Y));                                   % phase
            clear temp
            
            % for flap 2
            Y2=fft(detrend(flapPosF2,'linear'),L2)./L2;                     % fourier transform
            temp=Y2(1:L2/2+1);
            temp(2:end-1)=2*temp(2:end-1);
            Y2=temp;
            Y2mag = abs(Y2);                                                % magnitude
            Y2ph = atan2(imag(Y2),real(Y2));                                % phase
            clear temp
            
            % for platform
            Y3=fft(detrend(platPosz,'linear'),L2)./L2;                      % fourier transform
            temp=Y3(1:L2/2+1);
            temp(2:end-1)=2*temp(2:end-1);
            Y3=temp;
            Y3mag = abs(Y3);                                                % magnitude
            Y3ph = atan2(imag(Y3),real(Y3));                                % phase
            clear temp
            
            % fft procedure for wave probe
            Y_p = fft(wg6,L2)./L2;                                          % fourier transform
            temp=Y_p(1:L2/2+1);
            temp(2:end-1)=2*temp(2:end-1);
            Y_p=temp;
            Y_pmag=abs(Y_p);                                                % magnitude
            Y_pph=atan2(imag(Y_p),real(Y_p));                               % phase
            clear temp
            
            % save spectra to data structure
            Config3Irr.final.flapPosF1.PosSpec{i,1}=Y;
            Config3Irr.final.flapPosF2.PosSpec{i,1}=Y2;
            Config3Irr.final.platPosz.PosSpec{i,1}=Y3;
            Config3Irr.final.flapPosF1.offset(i,1)=median(flapPosF1);
            Config3Irr.final.flapPosF2.offset(i,1)=median(flapPosF2);
            Config3Irr.final.platPosz.offset(i,1)=median(platPosz);
            Config3Irr.final.WaveSpec{i,1}=Y_p;
            Config3Irr.final.freqw{i,1}=PP;
            
            %RAO Calculation
            Config3Irr.final.flapPosF1.RAOmag{i,1} = Ymag./Y_pmag;          % Units: deg/m
            Config3Irr.final.flapPosF1.RAOph{i,1} = Yph-Y_pph;
            Config3Irr.final.flapPosF2.RAOmag{i,1} = Y2mag./Y_pmag;
            Config3Irr.final.flapPosF2.RAOph{i,1} = Y2ph-Y_pph;
            Config3Irr.final.platPosz.RAOmag{i,1} = Y3mag./Y_pmag;          % units: m/m
            Config3Irr.final.platPosz.RAOph{i,1} = Y3ph-Y_pph;
        end
    end
    queryName={'flapPosF1','flapPosF2','platPosz'};                         % desired signals to analyze
      for k3=1:length(queryName)

        for k=1:length(Hunq)
                goodidx=find(ic==k);
            for k2=1:length(goodidx)
                if Config3Irr.T(goodidx(k2))>2  % interpolate same significant wave height to common sample freq
                    combMagSpec{k,k3}(:,k2)=interp1(Config3Irr.final.freqw{2},...
                        Config3Irr.final.(queryName{k3}).RAOmag{goodidx(k2)},Config3Irr.final.freqw{1});
                    combPhSpec{k,k3}(:,k2)=interp1(Config3Irr.final.freqw{2}...
                        ,Config3Irr.final.(queryName{k3}).RAOph{goodidx(k2)},Config3Irr.final.freqw{1});
                    DSWave(:,k2)=interp1(Config3Irr.final.freqw{2}...
                        ,abs(Config3Irr.final.WaveSpec{goodidx(k2)}),Config3Irr.final.freqw{1});
                else
                    combMagSpec{k,k3}(:,k2)=Config3Irr.final.(queryName{k3}).RAOmag{goodidx(k2)};
                    combPhSpec{k,k3}(:,k2)=Config3Irr.final.(queryName{k3}).RAOph{goodidx(k2)};
                    DSWave(:,k2)=abs(Config3Irr.final.WaveSpec{goodidx(k2)});
                end
                                
                % threshold by wave forcing, send all below thresh %ile to zero
                ordn=ceil(thresh*length(DSWave(:,k2)));
                sortW=sort(DSWave(:,k2));
                Wthresh=sortW(ordn);
                lpidx=find(DSWave(:,k2)<Wthresh);
                badidx=find(isnan(DSWave(:,k2)));
                DSWave(badidx,k2)=0;                                        % occasionally the last value (dropped anyway) goes to NaN
                
                % weight by wave power at given frequency
                combMagSpec{k,k3}(:,k2)= combMagSpec{k,k3}(:,k2).*DSWave(:,k2);
                combMagSpec{k,k3}(badidx,k2)=0;
                combPhSpec{k,k3}(badidx,k2)=0;
                
                % take sin, cos of phase angles for angular average
                sincombPhSpec{k,k3}(:,k2)=sin(combPhSpec{k,k3}(:,k2)).*DSWave(:,k2);
                coscombPhSpec{k,k3}(:,k2)=cos(combPhSpec{k,k3}(:,k2)).*DSWave(:,k2);
                    
            end
            % assign temporary variable to wave frequency vector
            avfreqw=Config3Irr.final.freqw{1};
            
            % calculate weighted average
            meanMagSpec{k,k3}=sum(combMagSpec{k,k3},2)./sum(DSWave,2);
         
            % calculate average radian angle (prevents averaging error due
            % to pi-wrapping
            meanPhSpec{k,k3}=atan2(sum(sincombPhSpec{k,k3},2)./sum(DSWave,2),...
                sum(coscombPhSpec{k,k3},2)./sum(DSWave,2));
            
             % bin average every n points; moving average filter (basically
             % 'smooth', but does not require DSP package)
            avmeanPhSpec{k,k3}=filtFB(smoothNum,1,meanPhSpec{k,k3},[],2);
            avmeanMagSpec{k,k3}=filtFB(smoothNum,1,meanMagSpec{k,k3},[],2); 
            
            % cut by wave power threshold;
            avmeanMagSpec{k,k3}(lpidx)=[];
            avmeanPhSpec{k,k3}(lpidx)=[];
            avfreqw(lpidx)=[];
            
            % cut off wave power frequency to excited range
            hiidx=find(avfreqw > 50);
            avmeanMagSpec{k,k3}(hiidx)=[];
            avmeanPhSpec{k,k3}(hiidx)=[];
            avfreqw(hiidx)=[]; 
            
            % Save averaged, thresholded spectra
            Config3Irr.final.(queryName{k3}).RAOMagAvg{1,k}=avmeanMagSpec{k,k3};
            Config3Irr.final.(queryName{k3}).RAOphAvg{1,k}=avmeanPhSpec{k,k3};%           
            Config3Irr.final.freqAvgw{1,k}=avfreqw;
        end  
      end
    save(final_file,'Config3Irr');
else
    cd(final_folder);
    load(final_file);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_data == 1
    queryName={'flapPosF1','flapPosF2','platPosz'};
    numTrials=length(Config3Irr.T);
    [Hunq,ia,ic]=unique(Config3Irr.H);
    colorvec={'b','r','g','m','c'};                                         % color code by wave height
    for k3=1:length(queryName)
        % initialize figure
        figure(2*(k3-1)+1); clf;
        for k=1:length(Hunq)
            goodidx=find(ic==k);
            
            plot(Config3Irr.final.freqAvgw{1,k},...
                Config3Irr.final.(queryName{k3}).RAOMagAvg{1,k},'o','Color',colorvec{k})       % magnitude plot
            legvec{k}=strcat('H=',num2str(Hunq(k)));
            if k==1;
                hold on; grid on;
                xlabel('Frequency (rad/s)')
                ylabel('RAO (deg/m)')
                figure(2*(k3-1)+2); clf;
                yl=ylim;
                ylim([0 yl(2)])
            end
            
            figure(2*(k3-1)+2)
            polarplot(Config3Irr.final.(queryName{k3}).RAOphAvg{1,k},...
                Config3Irr.final.freqAvgw{1,k},'o','MarkerFaceColor',colorvec{k}) % phase plot
            if k==1;
                hold on; grid on;
                rlim([0 10]);

            end
            figure(2*(k3-1)+1)
            legend(legvec,'Location','NorthEast')
            figure(2*(k3-1)+2)
            legend(legvec,'Location','NorthEast')
            figure(2*(k3-1)+1);
         end
        figure(2*(k3-1)+1)
        title(queryName{k3})
        savefig(['Config2Irr_RAO_',queryName{k3},'.fig'])
        figure(2*(k3-1)+2)
        title(queryName{k3})
    end
    cd ..
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


