%Calculates force excitation coefficient for irregular wave cases. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear;
%% User Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
process_inter=1;    % 1:process *.txt to *.mat, 0: loads *.mat
process_final=1;    % 1:process final data structure, 0: loads final *.mat
plot_data=1;        % 1 to generate plots
dur=20;             % the number of wave cycles to consider
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Directory info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_folder= pwd;
final_folder='.\WaveExcitationIrr';
inter_folder='..\inter\WaveExcitationIrr';
log_folder = '../logs';
inter_file='WaveExcIrr.mat';
final_file='WaveExcIrr_final.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process intermediate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if process_inter==1
    % import tst log *.xlsx
    cd(log_folder)
    [num,txt,~]=xlsread('WECSIM2_WaveExcitationIrr.xlsx','Log');
    data.Exp='WaveExcIrr';
    data.Header=txt(6,2:end);
    data.Trial=num(:,2);
    data.Tp=num(:,3);
    data.Hm0=num(:,4);
    data.Tp_fs=num(:,5);
    data.Hm0_fs=num(:,6);
    data.Notes=txt(18:end,9);
    data.Flag=num(:,10);
    clear txt num raw
    % imports excel data
    
    % define summary fields
    WaveExcIrr.Flag=data.Flag;
    WaveExcIrr.Tp=data.Tp;
    WaveExcIrr.Hm0=data.Hm0;
    WaveExcIrr.Hm0_fs=data.Hm0_fs;
    WaveExcIrr.Tp_fs=data.Tp_fs;
    numTrials=length(WaveExcIrr.Tp);
    
    %   Note that all flagged cases are towards at the end of the time series.
    cd(inter_folder)            
    %% Load data from text files
    for k=1:numTrials
        if WaveExcIrr.Flag(k)==0;                                           % does not process flagged cases
            trial_str=sprintf('%02d', data.Trial(k));
            Trial=['Trial' trial_str];
            cd(['./' Trial]);
            contents=dir('*.txt');                                          % find text files
            for j=1:length(contents)
                sensor_file=contents(j).name;
                sensor_name=strrep(sensor_file,'.txt','');
                WaveExcIrr.inter.(Trial).(sensor_name)=load(sensor_file);    % load data here
            end
            date_utc=datetime(WaveExcIrr.inter.(Trial).time,'convertfrom','datenum');
            time_utc=datetime(date_utc,'TimeZone','UTC');
            datetime_local = datetime(time_utc,'TimeZone','America/Los_Angeles');
            datetime_local.Format = 'dd-MMM-yyyy';
            WaveExcIrr.inter.(Trial).Date=datetime_local(1);
            datetime_local.Format = 'HH:mm:ss';
            WaveExcIrr.inter.(Trial).TimeStart=datetime_local(1);
            WaveExcIrr.inter.(Trial).TimeEnd=datetime_local(end);
            timeStart=WaveExcIrr.inter.(Trial).time(1);
            WaveExcIrr.inter.(Trial).time=(WaveExcIrr.inter.(Trial).time-timeStart)*60*60*24; % convert time to seconds
            cd ..
            clear timeStart
        else
            continue
        end
    end
    save(inter_file,'WaveExcIrr');
    clearvars -except WaveExcIrr final_folder final_file inter_file inter_folder process_final plot_data base_folder dur
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
    Fnames=fieldnames(WaveExcIrr.inter);
    numFields=length(Fnames);
    dtFlag=zeros(numFields,1);
    for k=1:numFields                                                       % flags non-uniform sampling frequency
        dterr=max(abs(diff(WaveExcIrr.inter.(Fnames{k}).time)-0.01));
        if dterr > 10E-4
            dtFlag(k)=1;
        end
        dt(k)=mean(diff(WaveExcIrr.inter.(Fnames{k}).time));
    end
    
    %% Import system loads
    
    % Filter configuration
    filtNum=0.125*ones(8,1);                                                % filter coeffs. for FOSWEC sensors (100 Hz)
    filtNumW=0.25*ones(4,1);                                                % filter coeffs. for wave gauges (50 Hz)
    n=3;                                                                    % number of adjacent frequency bands to merge
    smoothNum=ones(n,1).*(1/n);
    thresh2=0.95;                                                           % wave spectral amplitude percentile (as a decimal) to retain
    
    % pre-allocation
    startIdx=zeros(numFields,1);
    peakLoc{numFields,1}=zeros(1,1);
    peakmag{numFields,1}=zeros(1,1);
    % load undisturbed wave time series
    cd('../../../WECSIM/final');
    load('WaveTuning');
    cd('../../WECSIM2/final/WaveExcitationIrr');
    it=1;                                                                   % iteration counter                                                          
    for k=1:length(WaveExcIrr.Tp)
        if dtFlag(it)< 1 && WaveExcIrr.Flag(k) < 1
            Tq=WaveExcIrr.Tp(k);                                            % query wave period
            Hq=WaveExcIrr.Hm0(k);                                           % query wave height
            goodidx=find(abs(Tq-WaveTuning.RandWave.T)<0.05);               % find calibrated wave runs matching T,H.
            goodidx2=find(abs(Hq-WaveTuning.RandWave.H)<0.01);
            goodidx=intersect(goodidx,goodidx2);
            if length(goodidx)>1                                            % if multiple matching runs, take the first
                goodidx=goodidx(1);
            end
            
            
            %% Wave elevations
            fname=strcat('Test',num2str(goodidx));
            TS=WaveTuning.RandWave.Tests.(fname).time;
            dtW(it)=mean(diff(TS));                                          % sampling interval for waves
            WaveExcIrr.final.TimeDomain.(Fnames{it}).WaveTS=filtFB...
                (filtNumW,1,detrend(WaveTuning.RandWave.Tests.(fname).wg6,'linear'),[],2); % filter wave gauge data
            Tn=Tq./dt(it);                                                   % number of FOSWEC samples (100 Hz) per wave period
            Tw=Tq./dtW(it);                                                  % number of wave gauge samples (50 Hz) per wave period
            
            %% Trim start times
            thresh=0.5*max(WaveExcIrr.final.TimeDomain.(Fnames{it}).WaveTS);    % specify minimum peak height
            [peakLoc{it,1},peakmag{it,1}]=peakfinder(WaveExcIrr.final.TimeDomain.(Fnames{it}).WaveTS... % find peaks above specified height
                ,[],thresh,1,false,false);
            startidx=peakLoc{it,1}(2);                                      % start calculations from second peak, ensuring wave field fully developed
            
            % find start time, use to find startidx for forces
            tW=startidx*dtW;
            
            % calculate end times: dur wave cycles
            WaveExcIrr.final.TimeDomain.(Fnames{it}).WaveTS=...
                WaveExcIrr.final.TimeDomain.(Fnames{it}).WaveTS(startidx:round(startidx+Tw*dur));

            
            %% Fourier transform wave data
            Lw=length(WaveExcIrr.final.TimeDomain.(Fnames{it}).WaveTS);     % length of retained wave time series
            Lw=2^nextpow2(Lw);
            P2=fft(WaveExcIrr.final.TimeDomain.(Fnames{it}).WaveTS,Lw)./Lw; % fourier transform
            temp=P2(1:Lw/2+1);
            temp(2:end-1)=2*temp(2:end-1);
            WaveExcIrr.final.FreqDomain.(Fnames{it}).WaveF=temp;
            clear temp;
            freq=(1/dtW(it))*(0:(Lw/2))./Lw;                                % wave frequency domain
            freqw{it}(:,1)=freq*2*pi;                                       % convert to radian frequency
            Lw=length(WaveExcIrr.final.TimeDomain.(Fnames{it}).WaveTS);
            %% Loads
            
            for m=1:18                                                      % 6 modes on 3 bodies -> 18 total to consider
                endname=strrep(FnameString{m},'lc','');                     % deletes lc from the processed field name.
                coeffname=strcat('c',endname);                              % field name for excitation coefficients
                magname=strcat(coeffname,'mag');
                phname=strcat(coeffname,'phase');
                
                % invert load cells as needed based upon prescribed
                % coordinate system
                if m==4|| m==7 || m==3 || m==14 || m==16;                   % flip sign of y, z directions for flap 1 load cell
                    WaveExcIrr.final.TimeDomain.(Fnames{it}).(endname)=filtFB...
                        (filtNum,1,-1*WaveExcIrr.inter.(Fnames{it}).(FnameString{m}),[],2);
                else
                    WaveExcIrr.final.TimeDomain.(Fnames{it}).(endname)=filtFB...
                        (filtNum,1,WaveExcIrr.inter.(Fnames{it}).(FnameString{m}),[],2);
                end
                L=length(WaveExcIrr.final.TimeDomain.(Fnames{it}).(endname));
                
                % find startidx
                TSf=[0:dt(it):dt(it).*L];
                [~,startidxf]=min(abs(tW(it)-TSf));
                
                % calculate end times: by wave cycles
                WaveExcIrr.final.TimeDomain.(Fnames{it}).(endname)=...
                    WaveExcIrr.final.TimeDomain.(Fnames{it}).(endname)([startidxf:round(startidxf+Tn*dur)]);
                L=length(WaveExcIrr.final.TimeDomain.(Fnames{it}).(endname));
                
                % interpolate force time series to wave sampling intervals
                WaveExcIrr.final.TimeDomain.(Fnames{it}).(endname)=interp1([0:dt(it):dt(it).*(L-1)]...
                    ,WaveExcIrr.final.TimeDomain.(Fnames{it}).(endname),[0:dtW(it):dtW(it).*(Lw-1)],'linear','extrap').';

                
                %% Fourier Transform force data
                L=length(WaveExcIrr.final.TimeDomain.(Fnames{it}).(endname)); % interpolated force data time series length
                L2=2^nextpow2(L);
                P2=fft(detrend(WaveExcIrr.final.TimeDomain.(Fnames{it}).(endname),'linear'),L2)./L2; % fourier transform
                temp=P2(1:L2/2+1); 
                temp(2:end-1)=2*temp(2:end-1);
                WaveExcIrr.final.FreqDomain.(Fnames{it}).(endname)=temp;
                clear temp
                freqf=(1/dtW(it))*(0:(L2/2))./L2;                           % force measurement frequency domain
                freqfw{it}(:,1)=freqf*2*pi;                                 % convert to radians
                
                
                %% Calculate excitation coefficient by dividing spectra
                
                % magnitude (from frequency domain estimate)
                WaveExcIrr.final.(magname){it,1}=abs(WaveExcIrr.final.FreqDomain.(Fnames{it}).(endname))./...
                    abs(WaveExcIrr.final.FreqDomain.(Fnames{it}).WaveF);
                
                % force phase
                tempf=atan2(imag(WaveExcIrr.final.FreqDomain.(Fnames{it}).(endname)),real(WaveExcIrr.final.FreqDomain.(Fnames{it}).(endname)));
                
                % wave phase
                tempw=atan2(imag(WaveExcIrr.final.FreqDomain.(Fnames{it}).WaveF),real(WaveExcIrr.final.FreqDomain.(Fnames{it}).WaveF));
                
                % relative phase to wave forcing
                WaveExcIrr.final.(phname){it,1}=tempw-tempf;
                clear tempf; clear tempw;
                
            end
            it=it+1;
        end
    end
    queryName={'cFxF1','cTyF1','cFxF2','cTyF2','cFzPlat'};                  % a field name specifying the phase/mag plot present (accepts a vector)
    [Hunq,ia,ic]=unique(WaveExcIrr.Hm0(WaveExcIrr.Flag < 1));
    temp=WaveExcIrr.Tp(WaveExcIrr.Flag < 1);
    for k3=1:length(queryName);
        figure(2*(k3-1)+1); clf;
        for k=1:length(Hunq)
            goodidx=find(ic==k);
            for k2=1:length(goodidx) % interpolate distinct wave period fft data to common frequency resolution
                if temp(goodidx(k2)) > 2
                    combMagSpec{k,k3}(:,k2)=interp1(freqw{2},WaveExcIrr.final.(strcat(queryName{k3},'mag')){goodidx(k2),1},freqw{1});
                    combPhSpec{k,k3}(:,k2)=interp1(freqw{2},WaveExcIrr.final.(strcat(queryName{k3},'phase')){goodidx(k2),1},freqw{1});
                    DSWave(:,k2)=interp1(freqw{2},abs(WaveExcIrr.final.FreqDomain.(Fnames{goodidx(k2)}).WaveF),freqw{1});
                    
                else
                    combMagSpec{k,k3}(:,k2)=WaveExcIrr.final.(strcat(queryName{k3},'mag')){goodidx(k2),1};
                    combPhSpec{k,k3}(:,k2)=WaveExcIrr.final.(strcat(queryName{k3},'phase')){goodidx(k2),1};
                    DSWave(:,k2)=abs(WaveExcIrr.final.FreqDomain.(Fnames{goodidx(k2)}).WaveF);
                end
                
                % threshold by wave forcing spectral amplitude 
                ordn=ceil(thresh2*length(DSWave(:,k2)));
                sortW=sort(DSWave(:,k2));
                Wthresh=sortW(ordn);
                lpidx=find(DSWave(:,k2)<Wthresh);
                badidx=find(isnan(DSWave(:,k2)));                           % also delete any NaN indices
                lpidx=[lpidx;badidx];                                       % combines indices to send to zero
                DSWave(badidx,k2)=0;                                        % occasionally the last value (dropped anyway) goes to NaN
                
                % weight by wave power at given frequency
                combMagSpec{k,k3}(:,k2)=combMagSpec{k,k3}(:,k2).*DSWave(:,k2);
                combMagSpec{k,k3}(badidx,k2)=0;
                combPhSpec{k,k3}(badidx,k2)=0;
                
                % take sin, cos of angle for angular average calculations
                sincombPhSpec{k,k3}(:,k2)=sin(combPhSpec{k,k3}(:,k2)).*DSWave(:,k2);
                coscombPhSpec{k,k3}(:,k2)=cos(combPhSpec{k,k3}(:,k2)).*DSWave(:,k2);
            end
            % assign temporary variable for wave frequency
            avfreqw=freqw{1};
            
            % calculate weighted average
            meanMagSpec{k,k3}=sum(combMagSpec{k,k3},2)./sum(DSWave,2);
            
            % calculate average radian angle (prevents averaging error due
            % to pi-wrapping)
            meanPhSpec{k,k3}=atan2(sum(sincombPhSpec{k,k3},2)./sum(DSWave,2),...
                sum(coscombPhSpec{k,k3},2)./sum(DSWave,2));
            
            % bin average every n points
            avmeanPhSpec{k,k3}=filtFB(smoothNum,1,meanPhSpec{k,k3},[],2);
            avmeanMagSpec{k,k3}=filtFB(smoothNum,1,meanMagSpec{k,k3},[],2);
            
            % cut by wave power threshold
            avmeanMagSpec{k,k3}(lpidx)=[];
            avmeanPhSpec{k,k3}(lpidx)=[];
            avfreqw(lpidx)=[];
            
            % cut off wave power frequency to excited range
            hiidx=find(avfreqw > 50);
            avmeanMagSpec{k,k3}(hiidx)=[];
            avmeanPhSpec{k,k3}(hiidx)=[];
            avfreqw(hiidx)=[]; 
            
            % save the interpolated wave frequency (lower resolution) to the data structure
            WaveExcIrr.final.FreqDomain.freqw=freqw;
            
            % save the averaged, thresholded spectra
            WaveExcIrr.final.(queryName{k3}).MagAvg{1,k}=avmeanMagSpec{k,k3};
            WaveExcIrr.final.(queryName{k3}).phAvg{1,k}=avmeanPhSpec{k,k3};
            WaveExcIrr.final.FreqDomain.freqAvgw{1,k}=avfreqw;
        end
    end
    %     cd(final_folder)
    save(final_file, 'WaveExcIrr');
else
    cd(final_folder);
    load(final_file)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Merge spectra, calculate weighted averages

if plot_data==1;
    queryName={'cFxF1','cTyF1','cFxF2','cTyF2','cFzPlat'};                  % a field name specifying the phase/mag plot present (accepts a vector)
    [Hunq,ia,ic]=unique(WaveExcIrr.Hm0(WaveExcIrr.Flag < 1));
    colorvec={'b','r','g','m','c'};                                         % color coding is by wave height
    for k3=1:length(queryName)
        figure(2*(k3-1)+1); clf;
        for k=1:length(Hunq)
            goodidx=find(ic==k);                                            
            plot(WaveExcIrr.final.FreqDomain.freqAvgw{1,k}, WaveExcIrr.final... % magnitude plot
                .(queryName{k3}).MagAvg{1,k},'LineWidth',1.2,'Color',colorvec{k})
            legvec{k}=strcat('H=',num2str(Hunq(k)));
            if k==1;
                hold on; grid on;
                xlabel('Frequency (rad/s)')
                ylabel('Magnitude')
                figure(2*(k3-1)+2); clf;
            end
            figure(2*(k3-1)+2)                                              % phase plot
            polarplot(WaveExcIrr.final.(queryName{k3}).phAvg{1,k},WaveExcIrr.final...
                .FreqDomain.freqAvgw{1,k},'o','MarkerFaceColor',colorvec{k})
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
        title(['Mag' queryName{k3}])
        figure(2*(k3-1)+2)
        title(['\Phi' queryName{k3}])
    end
end

