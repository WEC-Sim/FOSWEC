% Dominic Forbush 2017
% Process the wave tuning cases into a single data structure. Check error
% between repeated cases on provided plots.

% Run process_inter for each sub-directory (e.g. RegularWaveTuning#) prior
% to running process_final.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear;
%% User Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
process_inter=1; %1:process *.txt to *.mat for inter_folder directory, 0: loads *.mat
process_final=0; % combine all (REG/IRR) wave tuning runs into a single structure. Run only after all wave case directories have been run with process_inter.
plot_final=0; % plot runs from process_final with common settings to check repeatability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Directory info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
homeDir=pwd;                                                                % home directory (facilitates relative paths)
final_folder='./final/';                                                    % location where final_file will be saved (only used if process_final==1)
inter_folder='./inter/RegularWaveTuning1';                                   % update inter_folder for wave case directory. This is where inter_file is saved
inter_file='WaveDat.mat';                                                   % intermediate '.mat' files containing all data for single wave case
final_file='WaveTuning';                                                    % name of final structure containing all wave data
%   Generated iff WaveDat.mat exists for each wave case directory and process_final==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process intermediate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if process_inter==1
    
    % Create output folder if it doesn't exist
    if exist(final_folder) == 0
        mkdir(final_folder)
    end
    
    %% Load data from text files
    cd(inter_folder);
    numTrials=length(dir(pwd))-3;
    for k=1:numTrials
        % does not process flagged cases
        trial_str=sprintf('%02d', k);
        Trial=['Trial' trial_str];
        cd(['./' Trial]);
        contents=dir('*.txt');                                              % find text files
        for j=1:length(contents)
            sensor_file=contents(j).name;
            if j==1;
                tmp=fileread(contents(j).name);
                Hstr=regexp(tmp,'H=\d.\d\d\d','match');
                if isempty(Hstr)
                    Hstr=regexp(tmp,'H=\d.\d\d','match');
                end
                Tstr=regexp(tmp,'T=\d.\d\d','match');
                if isempty(Tstr);
                    Tstr=regexp(tmp,'T=\d.\d','match');
                end
                Hstr=strrep(Hstr,'H=','');
                Tstr=strrep(Tstr,'T=','');
                WaveDat.H(k,1)=str2double(Hstr);
                WaveDat.T(k,1)=str2double(Tstr);
            end
            sensor_name=strrep(sensor_file,'.txt','');
            WaveDat.inter.(Trial).(sensor_name)=load(sensor_file);          % load data here
        end
        date_utc=datetime(WaveDat.inter.(Trial).time,'convertfrom','datenum');
        time_utc=datetime(date_utc,'TimeZone','UTC');
        datetime_local = datetime(time_utc,'TimeZone','America/Los_Angeles');
        datetime_local.Format = 'dd-MMM-yyyy';
        WaveDat.inter.(Trial).Date=datetime_local(1);
        datetime_local.Format = 'HH:mm:ss';
        WaveDat.inter.(Trial).TimeStart=datetime_local(1);
        WaveDat.inter.(Trial).TimeEnd=datetime_local(end);
        timeStart=WaveDat.inter.(Trial).time(1);
        WaveDat.inter.(Trial).time=(WaveDat.inter.(Trial).time-timeStart)*60*60*24; % convert time to seconds
        cd ..
        clear timeStart
        
    end
    save(inter_file,'WaveDat');
    cd(homeDir)
    clearvars WaveDat
else
    cd(homeDir)
    load(fullfile(final_folder,final_file))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Combine Trials Into Single, Large Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if process_final==1;
    cd(inter_folder)
    cd('..')
    Runs=dir(pwd);
    randidx=zeros(length(Runs),1);
    regidx=zeros(length(Runs),1);
    for k=1:length(Runs)                                                    % find directories of Random and Regular wave cases
        if ~isempty(regexp(Runs(k).name,'Random')) && Runs(k).isdir > 0
            randidx(k)=regexp(Runs(k).name,'Random');
        end
        if ~isempty(regexp(Runs(k).name,'Regular')) && Runs(k).isdir > 0
            regidx(k)=regexp(Runs(k).name,'Regular');
        end
    end
    cntrnd=0;
    cntreg=0;
    for k=1:length(randidx);
        if randidx(k)>0;
            cntrnd=cntrnd+1;
            cd(strcat('./',Runs(k).name));
            load('WaveDat');
            m=length(WaveDat.H);                                            % Extracts peak height, period for each trial and adds to summary structure.
            if cntrnd==1
                idx=m;
                WaveTuning.RandWave.H(1:m,1)=WaveDat.H(:);
                WaveTuning.RandWave.T(1:m,1)=WaveDat.T(:);
                for k2=1:m
                    Fname=strcat('Test',num2str(k2));
                    trial_str=sprintf('%02d', k2);
                    Qname=['Trial' trial_str];
                    WaveTuning.RandWave.Tests.(Fname)=WaveDat.inter.(Qname);
                end
                idx=idx+1;
            else
                
                WaveTuning.RandWave.H(idx:idx+m-1)=WaveDat.H(:);
                WaveTuning.RandWave.T(idx:idx+m-1)=WaveDat.T(:);
                for k2=1:m
                    Fname=strcat('Test',num2str(idx));
                    trial_str=sprintf('%02d', k2);
                    Qname=['Trial' trial_str];
                    WaveTuning.RandWave.Tests.(Fname)=WaveDat.inter.(Qname);
                    idx=idx+1;
                end
            end
        elseif regidx(k)>0
            cntreg=cntreg+1;
            cd(strcat('./',Runs(k).name));
            load('WaveDat');
            m=length(WaveDat.H);
            if cntreg==1
                idx=m;
                WaveTuning.RegWave.H(1:m,1)=WaveDat.H(:);
                WaveTuning.RegWave.T(1:m,1)=WaveDat.T(:);
                for k2=1:m
                    Fname=strcat('Test',num2str(k2));
                    trial_str=sprintf('%02d', k2);
                    Qname=['Trial' trial_str];
                    WaveTuning.RegWave.Tests.(Fname)=WaveDat.inter.(Qname);
                end
                idx=idx+1;
            else
                WaveTuning.RegWave.H(idx:idx+m-1)=WaveDat.H(:);
                WaveTuning.RegWave.T(idx:idx+m-1)=WaveDat.T(:);
                for k2=1:m
                    Fname=strcat('Test',num2str(idx));
                    trial_str=sprintf('%02d', k2);
                    Qname=['Trial' trial_str];
                    WaveTuning.RegWave.Tests.(Fname)=WaveDat.inter.(Qname);
                    idx=idx+1;
                end
                
            end
            
        else
            continue
        end
        cd('../');
    end
    cd(homeDir)
    cd(final_folder)
    save('WaveTuning',final_file);
else
    cd(homeDir)
    cd(final_folder)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Identical Runs, calculate error.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_final==1;
    
    load(final_file);
    % find unique periods, wave heights for regular and random waves
    [Tun,ia,ic]=unique(WaveTuning.RegWave.T);
    [Hun,ib,ie]=unique(WaveTuning.RegWave.H);
    [Tpun,ka,kc]=unique(WaveTuning.RandWave.T);
    [Hsun,kb,ke]=unique(WaveTuning.RandWave.H);
    %% Regular waves
    for k=1:length(Tun);
        goodidx=find(ic==k);
        goodidx2=find(ie==k);
        goodidx=intersect(goodidx,goodidx2);
        if ~isempty(goodidx)
            figure(k); clf;
            for k2=1:length(goodidx);
                fname=strcat('Test',num2str(goodidx(k2)));
                plot(WaveTuning.RegWave.Tests.(fname).time,...
                    WaveTuning.RegWave.Tests.(fname).wg6);
                if k2==1;
                    xlabel('Time (s)')
                    ylabel('WG6 Displacement (m)')                          % Wave Gauge 6 was at the prospective device location
                    grid on;
                    hold on;
                    tstr= ['Reg. Wave', ' ', 'H=',num2str(Hun(k)), ' ',...
                        'T=',num2str(Tun(k)), ' ', ...
                        'Tests=',num2str(length(goodidx))];
                    title(tstr);
                end
            end
        end
    end
    maxFigNum=k;
    clear goodidx goodidx2
    
    %% Irregular waves
    for k=1:length(Tpun)
        goodidx=find(kc==k);
        goodidx2=find(ke==k);
        goodidx=intersect(goodidx,goodidx2);
        if ~isempty(goodidx)
            figure(k+maxFigNum); clf;
            for k2=1:length(goodidx);
                fname=strcat('Test',num2str(goodidx(k2)));
                plot(WaveTuning.RandWave.Tests.(fname).time,...
                    WaveTuning.RandWave.Tests.(fname).wg6);
                if k2==1;
                    xlabel('Time (s)')
                    ylabel('WG6 Displacement (m)')
                    grid on;
                    hold on;
                    tstr= ['Rand. Wave', ' ','Hs=',num2str(Hsun(k)), ' ',...
                        'Tp=',num2str(Tpun(k)), '  ',...
                        'Tests=',num2str(length(goodidx))];
                    title(tstr);
                end
            end
        end
    end
end

