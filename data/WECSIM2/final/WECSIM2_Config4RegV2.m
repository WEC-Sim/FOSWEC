
% Update 8/16/2018: Uses the "Heave", "Surge", and "Pitch" sensor names
% from the 'inter' structure to perform the RAO calculation. These seem to
% have a more realistic trend.
% V2: Uses FOSWEC mounted sensors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all;
%close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
process_inter = 1;   % 1:processes inter data, 0:loads inter *.mat
process_final = 1;   % 1:processes final data, 0:loads final *.mat
plot_data =1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Directory Info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_folder = pwd;
final_folder = './Config4Reg';
inter_folder = '../inter/Config4Reg';
log_folder = '../logs';
inter_file = 'Config4Reg_inter.mat';
final_file = 'Config4Reg_final.mat';
addpath(genpath(strrep(pwd,'\WECSIM2\final\Config4Reg','')))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import Test Log and 'inter' Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if process_inter == 1
    % import test log *.xlsx
    cd(log_folder)
    [num,txt,raw] = xlsread('WECSIM2_Config4Reg.xlsx','Log');
    data.Exp = 'Config4Reg';
    data.Header = txt(13,2:end);
    data.Trial =num(:,2);
    data.T =    num(:,3);
    data.H =    num(:,4);
    data.T_fs = num(:,5);
    data.H_fs = num(:,6);
    data.Notes = txt(:,9);
    data.Flag = num(:,10);
    clear txt num raw
    % brings in excel data, adds to Config4Reg the period, height, scaled
    % period and height, notes and then the flag
    
    Config4Reg.Flag       = data.Flag;
    
    cd(inter_folder)
    % import 'inter' data for all trials
    numTrials = length(data.Trial);
    %flag count
    count = 0;
    for i = 1:numTrials
        if Config4Reg.Flag(i) == 1
            count=count+1;
        elseif Config4Reg.Flag(i) == 0
            continue
        end
    end
    
    % added to remove rows if flagged from Config4Reg.T,H
    for i = 1:numTrials-count
        if Config4Reg.Flag(i) == 1
            data.H(i) = [];
            data.T(i) = [];
            data.Notes(i) = [];
            data.Flag(i)=[];
            data.T_fs(i)=[];
            data.H_fs(i)=[];
            data.Trial(i)=[];
            Config4Reg.Flag = data.Flag;
        elseif Config4Reg.Flag(i) == 0
            continue
        end
    end
    
    Config4Reg.T          = data.T;
    Config4Reg.H          = data.H;
    numTrials = length(data.Trial);
    
    for i = 1:numTrials
        
        if data.Flag(i)==0
            trial_str = sprintf('%02d',data.Trial(i));
            Trial = ['Trial' trial_str];
            cd(['./' Trial])
            Config4Reg.inter.(Trial).T          = data.T(i);
            Config4Reg.inter.(Trial).H          = data.H(i);
            Config4Reg.inter.(Trial).Notes      = data.Notes(i);
            Config4Reg.inter.(Trial).Flag       = data.Flag(i);
            Config4Reg.inter.(Trial).T_fs       = data.T_fs(i);
            Config4Reg.inter.(Trial).H_fs       = data.H_fs(i);
            contents = dir( '*.txt' );
        end
        if data.Flag(i)==1
            continue
        end
        for j = 1:length(contents)                                          % import all sensor data
            sensor_file = contents(j).name;
            sensor_name = strrep(sensor_file,'.txt','');
            Config4Reg.inter.(Trial).(sensor_name) = load(sensor_file);     % load data here
        end
        date_utc = datetime(Config4Reg.inter.(Trial).time,'convertfrom','datenum');
        time_utc = datetime(date_utc,'TimeZone','UTC');
        datetime_local = datetime(time_utc,'TimeZone','America/Los_Angeles');
        datetime_local.Format = 'dd-MMM-yyyy';
        Config4Reg.inter.(Trial).Date = datetime_local(1);
        datetime_local.Format = 'HH:mm:ss';
        Config4Reg.inter.(Trial).TimeStart = datetime_local(1);
        Config4Reg.inter.(Trial).TimeEnd = datetime_local(end);
        timeStart = Config4Reg.inter.(Trial).time(1);
        Config4Reg.inter.(Trial).time = (Config4Reg.inter.(Trial).time - timeStart)*60*60*24;
        cd ..
    end
    save(inter_file,'Config4Reg')                                          % save 'inter' data to *.mat
    clearvars -except Config4Reg final_folder final_file inter_file inter_folder process_final plot_data base_folder
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
    Fnames=fieldnames(Config4Reg.inter);
    numTrials=length(Fnames);
    dtFlag=zeros(numTrials,1);
    filtNum=0.125*ones(8,1);
    filtNumW=0.25*ones(4,1);
    n=3; smoothNum=(1/n).*ones(n,1);% adjacent frequency bands to merge in spectra
    % by inspection: when body motion begins
    i_min = 5000.*ones(numTrials,1);
    queryVec={'flapPosF1','flapPosF2','platPosz','platPosx','platPosRy'};    % desired signals to analyze
    
    % load wave calibration time series for this operating condition
    cd('../../../WECSIM/final');
    load('WaveTuning');
    cd('../../WECSIM2/final/Config4Reg');
    
    for i = 1:numTrials
        Tq=Config4Reg.T(i);                                                 % wave period
        Hq=Config4Reg.H(i);                                                 % wave height
        Trial = Fnames{i};
        
        Fs = 1./mean(diff(Config4Reg.inter.(Trial).time));                 % sampling frequency
        dt(i) = 1./Fs;
        dterr=max(abs(diff(Config4Reg.inter.(Trial).time)-dt(i)));         % sampling period
        
        if dterr>10E-4
            dtFlag(i)=1;
        end
        if dtFlag(i)<1 %&& Config4Reg.inter.(Trial).flag<1;
            time = Config4Reg.inter.(Trial).time;
            Tn=round(Tq./dt(i));                                               % number of samples per wave period
            
            % import wave gauge data from calibration run
            goodidx=find(abs(Tq-WaveTuning.RegWave.T)<0.05);               % find calibrated wave runs matching T,H.
            goodidx2=find(abs(Hq-WaveTuning.RegWave.H)<0.01);
            goodidx=intersect(goodidx,goodidx2);
            if length(goodidx)>1                                            % if multiple matching runs, take the first
                goodidx=goodidx(1);
            end
            
            % wave elevations
            fname=strcat('Test',num2str(goodidx));
            TS=WaveTuning.RegWave.Tests.(fname).time;
            dtW(i)=mean(diff(TS));                                          % sampling interval for waves
            Config4Reg.final.WaveTS{i,1}=filtFB(filtNumW,1,detrend(WaveTuning.RegWave.Tests.(fname).wg6,'constant'),[],2);
            wg6= Config4Reg.final.WaveTS{i,1};
            % restrict wave time series to inquiry window
            [~,Widx(i)]=min(abs(TS-time(i_min(i))));
            Tw=round(Tq./dtW(i));
            wg6=wg6(Widx(i):Widx(i)+dur.*Tw);
            
            %% wave spectral calculations (for spectral RAO calc
            L=length(wg6);
            L2=2^nextpow2(L);
            Y_p = fft(wg6,L2)./L2;
            f = (1/dtW(i))*(0:(L2/2))/L2;                                   % frequency domain
            PP = 2*pi.*f;
            temp=Y_p(1:L2/2+1);
            temp(2:end-1)=2*temp(2:end-1);
            Y_p=temp;
            Y_pmag=abs(Y_p);
            Y_pmag=filtFB(smoothNum,1,Y_pmag,[],2);
            Y_pph=atan2(imag(Y_p),real(Y_p));
            [maxMag(i),idxDom(i)]=max(Y_pmag);                             % find the (single) excited wave frequency
            wDom(i)=PP(idxDom(i));                                         % excited radian wave frequency
            Config4Reg.final.WaveSpec{i,1}=Y_p;
            Config4Reg.final.Wavefreq{i,1}=PP;
            
            clear temp
            
            for k=1:length(queryVec);
                switch queryVec{k}
                    case 'flapPosF1'
                        varname = filtFB(filtNum,1,Config4Reg.inter.(Trial).flapPosF1,[],2); % retain trend
                        thresh=[];
                    case 'flapPosF2'
                        varname = filtFB(filtNum,1,Config4Reg.inter.(Trial).flapPosF2,[],2);
                        thresh=[];
                    case 'platPosz'
                        varname = filtFB(filtNumW,1,Config4Reg.inter.(Trial).platPosz,[],2);
                        %tidx=0:0.02:0.02*(length(varname)-1);
                        %varname=interp1(tidx,varname,Config4Reg.inter.(Trial).time); % interpolate to faster time. this may shift by up to 0.01 seconds!
                        thresh=(max(varname)-min(varname))/2;
                    case 'platPosx'
                        varname = filtFB(filtNum,1,Config4Reg.inter.(Trial).platPosx,[],2);
                        %tidx=0:0.02:0.02*(length(varname)-1);
                        %varname=interp1(tidx,varname,Config4Reg.inter.(Trial).time); % interpolate to faster time. this may shift by up to 0.01 seconds!
                        thresh=(max(varname)-min(varname))/2;
                    case 'platPosRy'
                        varname = filtFB(filtNum,1,Config4Reg.inter.(Trial).platPosRy,[],2);
                        %tidx=0:0.02:0.02*(length(varname)-1);
                        %varname=interp1(tidx,varname,Config4Reg.inter.(Trial).time); % interpolate to faster time. this may shift by up to 0.01 seconds!
                        thresh=(max(varname)-min(varname))/2;
                end
                
                %% clip to 'steady state' data
                varname=varname(i_min:i_min+dur.*Tn+1);
                time_c=time(i_min:i_min+dur.*Tn+1);
                %% Time-domain RAO Calculation: find peaks
                [peakLoc_max,peakMag_max]=peakfinder(varname,thresh,-15,1,false,false);
                [peakLoc_min,peakMag_min]=peakfinder(varname,thresh,15,-1,false,false);
                offset=(median(peakMag_max)-median(peakMag_min))/2 + median(peakMag_min); % check for assymetry
                
                % calculate statistics
                magmean=median(peakMag_max)-median(peakMag_min);
                magstd=std(peakMag_max)+std(peakMag_min);
                Config4Reg.final.(queryVec{k}).RAOmean(i,1)=magmean./Hq;
                Config4Reg.final.(queryVec{k}).RAOstd(i,1)=magstd./Hq;
                Config4Reg.final.(queryVec{k}).offset(i,1)=offset;
                
                % save timeseries
                Config4Reg.final.(queryVec{k}).TS{i,1}=varname;
                Config4Reg.final.time{i,1}=time_c;
                
                %% Spectral analysis- for comparison (should agree roughly with the predicted mean)
                
                % fft procedure
                L=length(varname);
                L2=2^nextpow2(L);
                Y=fft(detrend(varname,'linear'),L2)./L2;
                f = (1/dt(i))*(0:(L2/2))/L2;
                freqw=2*pi*f;
                temp=Y(1:L2/2+1);
                temp(2:end-1)=2*temp(2:end-1);
                Y=temp;
                Ymag = abs(Y);
                Ymag=filtFB(smoothNum,1,Ymag,[],2);
                Yph= atan2(imag(Y),real(Y));
                
                % save spectra
                Config4Reg.final.(queryVec{k}).PosSpec{i,1}=Y;
                
                % RAO calculation at dominant frequency
                [~,fidx]=min(abs(wDom(i)-freqw));
                Config4Reg.final.(queryVec{k}).RAOmag(i,1)=Ymag(fidx)./maxMag(i);
                Config4Reg.final.(queryVec{k}).RAOph(i,1)=Yph(fidx)- Y_pph(idxDom(i));
                Config4Reg.final.freqw{:,i}=freqw;
                
            end
        end
    end
    save(final_file,'Config4Reg')
    
else
    cd(final_folder)
    load(final_file)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_data == 1
    [Hunq,ia,ic]=unique(Config4Reg.H);
    queryVec={'flapPosF1','flapPosF2','platPosz','platPosx','platPosRy'};    % desired signals to analyze
    colorvec={'b','r','g','k'};
    labelVec={'RAO deg/m','RAO deg/m','RAO m/m','RAO m/m','RAO deg/m'};
    for k=1:length(queryVec)
        % color code based upon Hsig
        for k2=1:length(Hunq)
            goodidx=find(ic==k2);
            figure(2*k-1);
            errorbar(2*pi./Config4Reg.T(goodidx),Config4Reg.final.(queryVec{k}).RAOmean(goodidx),Config4Reg.final.(queryVec{k}).RAOstd(goodidx)...
                ,'s','Color',colorvec{ic(goodidx(1))},'MarkerEdgeColor',colorvec{ic(goodidx(1))});
            hold on
            %scatter(2*pi./Config4Reg.T(goodidx),Config4Reg.final.(queryVec{k}).RAOmag(goodidx),'^',...
             %   colorvec{ic(goodidx(1))});
            if k2==1;
            xlabel('Freq (rad/s)')
            ylabel(labelVec{k})
            grid on
            end
            savefig(['Config4Reg_RAO_',queryVec{k},'.fig'])
            
            figure(2*k);
            polarplot(Config4Reg.final.(queryVec{k}).RAOph(goodidx),2*pi./Config4Reg.T(goodidx),'o',...
                'MarkerEdgeColor',colorvec{ic(goodidx(1))},'MarkerFaceColor',colorvec{ic(goodidx(1))})
            if k2==1;
                hold on
                grid on
            end           
        end
    end
    cd ..    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   