% Calculate wave excitation force coefficients for the regular wave cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear;
%% User Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
process_inter=1; %1:process *.txt to *.mat, 0: loads *.mat
process_final=1; % 1:process final data structure, 0: loads final *.mat
plot_data=1; % 1 to generate plots
dur=20; % number of wave cycles to consider
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Directory info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_folder = pwd;
final_folder='./WaveExcitationReg';
inter_folder='../inter/WaveExcitationReg';
log_folder = '../logs';
inter_file='WaveExcReg.mat';
final_file='WaveExcReg_final.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process intermediate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if process_inter==1
    % import tst log *.xlsx
    cd(log_folder)
    [num,txt,~]=xlsread('WECSIM2_WaveExcitationReg.xlsx','Log');
    data.Exp='WaveExcReg';
    data.Header=txt(5,2:end);
    data.Trial=num(:,2);
    data.T=num(:,3);
    data.H=num(:,4);
    data.T_fs=num(:,5);
    data.H_fs=num(:,6);
    data.Notes=txt(20:end,9);
    data.Flag=num(:,10);
    clear txt num raw
    % imports excel data
    
    % define summary fields    
    WaveExcReg.Flag=data.Flag;
    WaveExcReg.T=data.T(data.Flag<1);
    WaveExcReg.H=data.H(data.Flag<1);
    WaveExcReg.H_fs=data.H_fs(data.Flag<1);
    WaveExcReg.T_fs=data.T_fs(data.Flag<1);
    numTrials=length(WaveExcReg.T);
    
    %   Note that all flagged cases are towards at the end of the time series.
    cd(inter_folder)        
    %% Load data from text files
    for k=1:numTrials
        if WaveExcReg.Flag(k)==0;                                           % does not process flagged cases
            trial_str=sprintf('%02d', data.Trial(k));
            Trial=['Trial' trial_str];
            cd(['./' Trial]);
            contents=dir('*.txt');                                          % find text files
            for j=1:length(contents)
                sensor_file=contents(j).name;
                sensor_name=strrep(sensor_file,'.txt','');
                WaveExcReg.inter.(Trial).(sensor_name)=load(sensor_file);    % load data here
            end
            date_utc=datetime(WaveExcReg.inter.(Trial).time,'convertfrom','datenum');
            time_utc=datetime(date_utc,'TimeZone','UTC');
            datetime_local = datetime(time_utc,'TimeZone','America/Los_Angeles');
            datetime_local.Format = 'dd-MMM-yyyy';
            WaveExcReg.inter.(Trial).Date=datetime_local(1);
            datetime_local.Format = 'HH:mm:ss';
            WaveExcReg.inter.(Trial).TimeStart=datetime_local(1);
            WaveExcReg.inter.(Trial).TimeEnd=datetime_local(end);
            timeStart=WaveExcReg.inter.(Trial).time(1);
            WaveExcReg.inter.(Trial).time=(WaveExcReg.inter.(Trial).time-timeStart)*60*60*24; % convert time to seconds
            cd ..
            clear timeStart
        else
            continue
        end
    end
    save(inter_file,'WaveExcReg');
    clearvars -except WaveExcReg final_folder final_file inter_file inter_folder process_final plot_data base_folder dur
    cd(base_folder)
else
    load([inter_folder,'\',inter_file]);                                    % load 'inter' *.mat
end
clear data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process Final Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check consistency of sampling frequency
if process_final==1;
    if exist(final_folder) ==  0
        mkdir(final_folder)
    end
    cd(final_folder)
    addpath('../../../functions')
    % select field names containing load data
    FnameString={'lcFxF1','lcFxF2','lcFxPlat','lcFyF1','lcFyF2','lcFyPlat','lcFzF1','lcFzF2'...
        ,'lcFzPlat','lcTxF1','lcTxF2','lcTxPlat','lcTyF1','lcTyF2','lcTyPlat','lcTzF1','lcTzF2'...
        ,'lcTzPlat'};
    Fnames=fieldnames(WaveExcReg.inter);
    numFields=length(Fnames);
    dtFlag=zeros(numFields,1);
    for k=1:numFields                                                       % flags non-uniform sampling frequency
        dterr=max(abs(diff(WaveExcReg.inter.(Fnames{k}).time)-0.01));
        if dterr > 10E-4
            dtFlag(k)=1;
        end
        dt(k)=mean(diff(WaveExcReg.inter.(Fnames{k}).time));
    end
    
    %% Import system loads
    
    % Filter configuration
    filtNum=0.125*ones(8,1);                                                 % filter coeffs. for FOSWEC sensors (100 Hz)
    filtNumW=0.25*ones(4,1);                                                 % filter coeffs. for wave gauges (50 Hz)

    % pre-allocation
    startIdx=zeros(numFields,1);
    peakLoc{numFields,1}=zeros(1,1);
    peakmag{numFields,1}=zeros(1,1);
    goodidx{numFields,1}=zeros(1,1);
    % load undisturbed wave time series
    cd('../../../WECSIM/final');
    load('WaveTuning');
    cd('../../WECSIM2/final/WaveExcitationReg');
    for k=1:numFields
        if dtFlag(k)==0;
            Tq=WaveExcReg.T(k);                                             % query wave period
            Hq=WaveExcReg.H(k);                                             % query wave height
            goodidx=find(abs(Tq-WaveTuning.RegWave.T)<0.05);               % find calibrated wave runs matching T,H.
            goodidx2=find(abs(Hq-WaveTuning.RegWave.H)<0.01);
            goodidx=intersect(goodidx,goodidx2);
            goodidx=goodidx(1);                                             % if multiple matching runs, take the first
            
            %% Wave elevations
            fname=strcat('Test',num2str(goodidx));
            TS=WaveTuning.RegWave.Tests.(fname).time;
            dtW(k)=mean(diff(TS));                                          % sampling interval for waves
            WaveExcReg.final.TimeDomain.(Fnames{k}).WaveTS=filtFB...
                (filtNum,1,detrend(WaveTuning.RegWave.Tests.(fname).wg6,'constant'),[],2);
            Tn=Tq./dt(k);                                                   % number of FOSWEC samples (100 Hz) per wave period
            Tw=Tq./dtW(k);                                                  % number of wave gauge samples (50 Hz) per wave period
            
            %% Trim start times
            thresh=0.75*max(WaveExcReg.final.TimeDomain.(Fnames{k}).WaveTS);    % specify minimum peak height
            [peakLoc{k,1},peakmag{k,1}]=peakfinder(WaveExcReg.final.TimeDomain.(Fnames{k}).WaveTS...
                ,[],thresh,1,false,false);
            startidx=peakLoc{k,1}(2);                                       % start calculations from second peak, ensuring fully-developed wave field 
            
            % find start time, use to find startidx for forces
            tW=startidx*dtW;
            
            % calculate end times: dur wave cycles
            WaveExcReg.final.TimeDomain.(Fnames{k}).WaveTS=...
                WaveExcReg.final.TimeDomain.(Fnames{k}).WaveTS(startidx:round(startidx+Tw*dur));
            
            % min/max of wave heights for amplitude estimate
            [maxWaveIdx,maxWave]=peakfinder(WaveExcReg.final.TimeDomain.(Fnames{k}).WaveTS,[],0.001,1,0);
            [minWaveIdx,minWave]=peakfinder(WaveExcReg.final.TimeDomain.(Fnames{k}).WaveTS,[],-0.001,-1,0);
            maxWave=maxWave(1:min(length(maxWave),length(minWave)));
            maxWaveIdx=maxWaveIdx(1:min(length(maxWave),length(minWave)));
            minWave=minWave(1:min(length(maxWave),length(minWave)));
            minWaveIdx=minWaveIdx(1:min(length(maxWave),length(minWave)));
            WaveAmp=maxWave-minWave;
             
            %% Fourier transform wave data
            
             Lw=length(WaveExcReg.final.TimeDomain.(Fnames{k}).WaveTS);     % length of retained wave time series
             Lw=2^nextpow2(Lw);
             P2=fft(WaveExcReg.final.TimeDomain.(Fnames{k}).WaveTS,Lw)./Lw; % fourier transform
             temp=P2(1:Lw/2+1);
             temp(2:end-1)=2*temp(2:end-1);
             WaveExcReg.final.FreqDomain.(Fnames{k}).WaveF=temp;            
             clear temp
             freq=(1/dtW(k))*(0:(Lw/2))./Lw;                                % wave frequency domain   
             freqw=freq*2*pi;                                               % convert to radians
                        
            % find dominant frequency peak (single for regular waves)
            [~,idxDom(k)]=max(abs(WaveExcReg.final.FreqDomain.(Fnames{k}).WaveF));
            wDom(k)=freqw(idxDom(k));
            
            
           
            %% Loads
            
            for m=1:18                                                      % 6 modes on 3 bodies
                endname=strrep(FnameString{m},'lc','');                     % deletes lc from the processed field name.
                coeffname=strcat('c',endname);                              % field name for excitation coefficients
                magname=strcat(coeffname,'mag');
                magnameTS=strcat(magname, 'TS');
                phname=strcat(coeffname,'phase');
                
                % invert load cells as needed based upon prescribed
                % coordinate system
                if m==4|| m==7 || m==3 || m==13 || m==16;                   % flip sign of y, z directions for flap 1 load cell
                WaveExcReg.final.TimeDomain.(Fnames{k}).(endname)=filtFB...
                    (filtNum,1,-1*detrend(WaveExcReg.inter.(Fnames{k}).(FnameString{m}),'constant'),[],2);
                else
                    WaveExcReg.final.TimeDomain.(Fnames{k}).(endname)=filtFB...
                    (filtNum,1,detrend(WaveExcReg.inter.(Fnames{k}).(FnameString{m}),'constant'),[],2); 
                end
                L=length(WaveExcReg.final.TimeDomain.(Fnames{k}).(endname));
                
                % find startidx
                TSf=[0:dt(k):dt(k).*L];
                [~,startidxf]=min(abs(tW(k)-TSf));
                
                % calculate end times: by wave cycles
                WaveExcReg.final.TimeDomain.(Fnames{k}).(endname)=...
                    WaveExcReg.final.TimeDomain.(Fnames{k}).(endname)([startidxf:round(startidxf+Tn*dur)]);
                L=length(WaveExcReg.final.TimeDomain.(Fnames{k}).(endname)); % length of retained force time series
                
                % find minima/maxima of loads
               [maxForceIdx,maxForce]=peakfinder(WaveExcReg.final.TimeDomain.(Fnames{k}).(endname),[],1,1,0);
               [minForceIdx,minForce]=peakfinder(WaveExcReg.final.TimeDomain.(Fnames{k}).(endname),[],-1,-1,0);
               maxForce=maxForce(1:min(length(maxForce),length(minForce)));
               minForce=minForce(1:min(length(maxForce),length(minForce)));
               maxForceIdx=maxForceIdx(1:min(length(maxForce),length(minForce)));
               minForceIdx=minForceIdx(1:min(length(maxForce),length(minForce)));
            
               %% Excitation coefficient calculation based on time-domain
               % maxima (at single excited frequency)
               ForceAmp=maxForce-minForce;
               ForceAmp=ForceAmp(1:min(length(WaveAmp),length(ForceAmp)));             
               WaveAmpN=WaveAmp(1:min(length(WaveAmp),length(ForceAmp)));
               excMedTS=median(ForceAmp./WaveAmpN);
               excStdTS=std(ForceAmp./WaveAmpN);
               
               % save the time series results
               WaveExcReg.final.(magnameTS)(k,1)=excMedTS;
               WaveExcReg.final.(magnameTS)(k,2)=excStdTS;
               
                %% Fourier Transform force data 
                
                L=2^nextpow2(L);
                P2=fft(WaveExcReg.final.TimeDomain.(Fnames{k}).(endname),L)./L; % fourier transform
                temp=P2(1:L/2+1);
                temp(2:end-1)=2*temp(2:end-1);
                WaveExcReg.final.FreqDomain.(Fnames{k}).(endname)=temp;
                clear temp
                freqf=(1/dt(k))*(0:(L/2))./L;                               % force measurement frequency domain
                freqfw=freqf*2*pi;                                          % convert to radians
                
                % find closest force frequency that agrees with dominant wave
                % freq.
                [minfreqDiff(m,k),idxDomf(m,k)]=min(abs(wDom(k)-freqfw));
                               
                %% Calculate excitation coefficient by dividing spectra
                
                % magnitude (from frequency domain estimate)
                WaveExcReg.final.(magname)(k,1)=abs(WaveExcReg.final.FreqDomain.(Fnames{k}).(endname)(idxDomf(m,k)))/...
                    abs(WaveExcReg.final.FreqDomain.(Fnames{k}).WaveF(idxDom(k)));
                
                % force phase
                tempf=atan2(imag(WaveExcReg.final.FreqDomain.(Fnames{k}).(endname)(idxDomf(m,k)))...
                    ,real(WaveExcReg.final.FreqDomain.(Fnames{k}).(endname)(idxDomf(m,k))));
                
                % wave phase
                tempw=atan2(imag(WaveExcReg.final.FreqDomain.(Fnames{k}).WaveF(idxDom(k)))...
                    ,real(WaveExcReg.final.FreqDomain.(Fnames{k}).WaveF(idxDom(k))));
                
                % relative phase to wave forcing
                 WaveExcReg.final.(phname)(k,1)=tempw-tempf;
                 clear tempf; clear tempw;
                
            end
        end
    end
    save(final_file, 'WaveExcReg');
else
    cd(final_folder);
    load(final_file)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_data==1;
    queryName={'cFxF1','cTyF1','cFxF2','cTyF2'}; % a field name specifying the phase/mag plot present (accepts a vector)
    [Tuniq,ia,ic]=unique(WaveExcReg.T);
    [Huniq,id,ie]=unique(WaveExcReg.H);
    for m=1:length(queryName)
        figure(m); clf;
        magname=strcat(queryName{m},'mag');
        phname=strcat(queryName{m},'phase');
        legpre=strrep(queryName{m},'c','');
        for k=1:length(Tuniq)
            goodidx=find(ic==k);
            numGood=length(goodidx);
            
            % offset to aid visualization of each wave height at common
            % frequencies
            off1=0;
            off2=0;
            off3=0;
            for k2=1:numGood;
                switch ie(goodidx(k2));
                    case 1
                        off=off1;
                        cString='b';
                        mType='-o';
                        subplot(2,1,1) % magnitude plot
                        ax1=plot(2*pi/(WaveExcReg.T(goodidx(k2)))+off,...
                            WaveExcReg.final.(magname)(goodidx(k2)),mType,'MarkerFaceColor',cString);
                        
                        subplot(2,1,2) % phase plot
                        ax12=polarplot(WaveExcReg.final.(phname)(goodidx(k2)),...
                            2*pi/(WaveExcReg.T(goodidx(k2)))+off,mType,'MarkerFaceColor',cString);
                        
                    case 2
                        off=off2;
                        cString='r';
                        mType='--^';
                        subplot(2,1,1)
                        ax2=plot(2*pi/(WaveExcReg.T(goodidx(k2)))+off,...
                            WaveExcReg.final.(magname)(goodidx(k2)),mType,'MarkerFaceColor',cString);
                        
                        subplot(2,1,2)
                        ax22=polarplot(WaveExcReg.final.(phname)(goodidx(k2)),...
                            2*pi/(WaveExcReg.T(goodidx(k2)))+off,mType,'MarkerFaceColor',cString);
                        
                    case 3
                        off=off3;
                        cString='g';
                        mType='-.s';
                        subplot(2,1,1)
                        ax3=plot(2*pi/(WaveExcReg.T(goodidx(k2)))+off,...
                            WaveExcReg.final.(magname)(goodidx(k2)),mType,'MarkerFaceColor',cString);
                        
                        subplot(2,1,2)
                        ax32=polarplot(WaveExcReg.final.(phname)(goodidx(k2)),...
                            2*pi/(WaveExcReg.T(goodidx(k2)))+off,mType,'MarkerFaceColor',cString);
                end
                
                if k2==1;
                    subplot(2,1,1)
                    xlabel('\omega(rad/s)')
                    ylabel(strcat('mag(',legpre,')'))
                    hold on;
                    grid on;
                    subplot(2,1,2)
                    hold on;
                    grid on;
                end
                switch ie(goodidx(k2));
                    case 1
                        off1=off1+0.04;
                    case 2
                        off2=off2+0.04;
                    case 3
                        off3=off3+0.04;
                end
            end
        end
        subplot(2,1,1)
        legend([ax1 ax2 ax3],{'H=0.015 m','H=0.045 m','H=0.136 m'},'Location','SouthEast')
        title(queryName{m});
    end    
end


