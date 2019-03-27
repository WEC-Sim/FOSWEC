% calculates period and per-unit-mass damping of surge free decay of FOSWEC

close all
clear

% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% directory information
homeDir=pwd;
load_folder = './inter/SurgeDecay';
output_variable_name = 'SurgeDecay.mat';
output_folder = './final/SurgeDecay';

% Test information (see test log)
trials = [8 10:20];
deltax = [repmat([0.15 0.20 0.10 .07],1,3)];                                % displacements (m)
dt=0.02;                                                                    % sample interval
ramp=0.5;                                                                   % retained time (s) before release

process_inter=0;                                                            % load raw .txt data files into structure
process_final=0;                                                            % process structure data
plot_data=1;                                                                % plot processes results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROCESS INTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if process_inter==1;
% Create output folder if it doesn't exist
if exist(output_folder) == 0
    mkdir(output_folder)
end

 % Create structure variable SurgeDecay
    for i = 1:length(trials)
        trial_str=sprintf('%02d',trials(i));
        TrialName{i}=strcat('Trial',trial_str);
        SurgeDecay.inter.(TrialName{i}).del_x = deltax(i);                  % initial displacement
        cd(load_folder)
        cd(['./' TrialName{i}])
        contents = dir( '*.txt' );
        for j = 1:length(contents)
            sensor_name = contents(j).name;
            sensor_name_notxt = strrep(sensor_name,'.txt','');
            SurgeDecay.inter.(TrialName{i}).(sensor_name_notxt) = load(sensor_name);  % load data here    
        end
        cd '../../..'
        date_utc = datetime(SurgeDecay.inter.(TrialName{i}).time,'convertfrom','datenum');
        SurgeDecay.inter.(TrialName{i}).datetime_utc = datetime(date_utc,'TimeZone','UTC');
        SurgeDecay.inter.(TrialName{i}).datetime_local = datetime(SurgeDecay.inter.(TrialName{i}).datetime_utc,'TimeZone','America/Los_Angeles');
        SurgeDecay.inter.(TrialName{i}).time_sec = (SurgeDecay.inter.(TrialName{i}).time - SurgeDecay.inter.(TrialName{i}).time(1))*60*60*24;
    end
       
cd(homeDir)  
cd(output_folder)    
save(output_variable_name,'SurgeDecay') % save structure variable to final folder
cd '../..'

else
    cd(homeDir)
    cd(output_folder)
    load(output_variable_name);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROCESS FINAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if process_final==1;
    
    filtNum=0.125*ones(8,1);                                                % FIR filter coefficients, a 12.5 Hz cut-off for 100 Hz sampling
    
    deltaxU = unique(deltax);
    cases = length(deltaxU);
    
     x_lim = [5000 5300; 5000 5300; 4800 5300;3600 3900; 4000 4300;4600 4800; ... % bound the approximate locations of platform release
        3600 4000;2400 2600; 2250 2500; 4100 4400;6150 6400;2900 3200];
    
    for i = 1:length(trials)
        trial_str=sprintf('%02d',trials(i));
        TrialName{i}=strcat('Trial',trial_str);
        der = gradient(SurgeDecay.inter.(TrialName{i}).platPosx,SurgeDecay.inter.(TrialName{i}).time_sec);    % calculate flap velocity
        der = filtFB(filtNum,1,der,[],2);                                   % filter velocity
        der2 = gradient(der,SurgeDecay.inter.(TrialName{i}).time_sec);      % calculate acceleration
        der2 = filtFB(filtNum,1,der2,[],2);                                 % filter acceleration
        [slope_min,i_slope_min] =  max(der2(x_lim(i,1):x_lim(i,2)));        % this finds the release point of the platform
        
        % correct time stamps relative to release point and mean position
        t_offset(i) = SurgeDecay.inter.(TrialName{i}).time_sec(x_lim(i)+i_slope_min);
        z_offset(i) = mean(SurgeDecay.inter.(TrialName{i}).platPosx(7000:8000)); % mean final (restored) position, after all motion has settled
        i_win_min = find(SurgeDecay.inter.(TrialName{i}).time_sec - t_offset(i) >= -0.5001);
        i_win_min = i_win_min(1);
        i_win_max = find(SurgeDecay.inter.(TrialName{i}).time_sec - t_offset(i) <= 3.5001);
        i_win_max = i_win_max(end);
        i_window_range(i,:) = [i_win_min i_win_max];                        % window over same relative time period
        t_corrected{i}(:,1) = SurgeDecay.inter.(TrialName{i}).time_sec(i_win_min:i_win_max) - t_offset(i); % shift time stamps to considered window
        z_corrected{i}(:,1) = (SurgeDecay.inter.(TrialName{i}).platPosx(i_win_min:i_win_max) - z_offset(i) )... % normalize position
            ./mean(SurgeDecay.inter.(TrialName{i}).platPosx(i_win_min:i_win_min+19) - z_offset(i)); % the average initial position
        
        lengthmin(i)=length(t_corrected{i}(:,1));
        
    end
    cutoff=min(lengthmin);                                                  % consider time window of common length for all trials
    
    for i=1:length(trials)
        t_corrected2{i}(:,1) = t_corrected{i}(1:cutoff,1);
        z_corrected2{i}(:,1) = z_corrected{i}(1:cutoff,1);
    end
    
    % convert t, z data format to matrix, now that lengths are consistent
    t_corrected2=cell2mat(t_corrected2);
    z_corrected2=cell2mat(z_corrected2);
    t_corrected=t_corrected2;
    z_corrected=z_corrected2;
    clear t_corrected2 z_corrected2;
    
   SurgeDecay.final.t = t_corrected(:,1);                                   % the common time stamps for all processed cases
    
    % check plot; plots each trial
    figure; clf;
    for i = 1:length(trials)
        plot(t_corrected(:,i),z_corrected(:,i)), hold on
        leg{i} = ['Trial ' num2str(trials(i),2)];
    end
    set(gcf, 'Color', 'w');
    grid on
    xlabel('t [s]')
    ylabel('Z/Z_o')
    axis([-0.5 8 -.75 1])
    legend(leg,'location','northeast')
    box off
    grid on
    title('Check Plot: Normalized Heave Decay trials')
    % check plot
    
    for i = 1:cases
        fname=strcat('X',num2str(deltaxU(i)));
        fname=strrep(fname,'.','pt');
        i_trial = trials(i);
        i_delz = find(deltax== deltaxU(i));
        cmean = zeros(i_window_range(1,2) - i_window_range(1,1),1);         % initialize mean position matrix
        cstd = zeros(i_window_range(1,2) - i_window_range(1,1),1);          % initialize postion std. dev matrix
        N(i) = length(i_delz);
        
        for j = 1:N(i)                                                      % for each repeated trial of given del_z, find average position, std of positions
            c3(:,j) = z_corrected(:,i_delz(j));
        end
        
        Cstd(:,1) = std(c3,0,2);                                            % measure of run-to-run variance
        Cmean(:,1) = mean(c3,2);                                            % mean of run-to-run position
        Cmean(1:(ramp/dt),1)=1;                                             % set all retained values prior to release of flap to 1
        
        
        % Log normalized data output data structure
        SurgeDecay.final.(fname).x_mean = Cmean;                            % mean position for all trials with given initial displacement
        SurgeDecay.final.(fname).x_std = Cstd;                              % std dev of position for all trials with given initial displacement
        SurgeDecay.final.(fname).tDistCoeff = tinvDom(0.95,N(i));           % the 95% confidence interval based on t-distribution of sample size N(i)
        
        % estimate of damping ratio and natural frequency. NOTE: this case
        % is over-damped, and thus the percent overshoot and the log
        % decrement method are invalid estimators. The canonical form of
        % damped oscillator is used for fitting.
        
        [params,fitError(i)]=fminsearch(@(x)fitFunc(x,SurgeDecay.final.t(ramp/dt:end),Cmean(ramp/dt : end)),[1 1]);
        wn(i)=params(1);
        zeta(i)=params(2);
        Tn(i)=2*pi/wn(i);
        wd(i)=wn(i)*sqrt(1-(zeta(i))^2);
        Td(i)=1/wd(i);
        
        %calculate damping-per-unit-mass based on wn and zeta
        c(i)=2*zeta(i)*wn(i);    
   % Save damping values and frequency estimates in output data
        % structure
        SurgeDecay.final.(fname).c=c(i);                                    % per unit mass/inertia damping
        SurgeDecay.final.(fname).zeta=zeta(i);                              % damping ratio
        SurgeDecay.final.(fname).Td=Td(i);                                  % damped period
        SurgeDecay.final.(fname).Tn=Tn(i);                                  % natural period
        
    end
    cd(homeDir)
    cd(output_folder)
    save(output_variable_name,'SurgeDecay') % save structure variable to final folder
    cd '../..'
    
else
    cd(homeDir)
    cd(output_folder)
    load(output_variable_name);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
fnames=fieldnames(SurgeDecay.final);                                        % fieldnames of final structure
x_vec=unique(deltax);                                                       % unique displacement values (m)

% plot average decay trajectory for each z displacement
for k=2:length(fnames)                                                      % the first field is time, all others are unique x
    figure(k)
    plot(SurgeDecay.final.t,SurgeDecay.final.(fnames{k}).x_mean,'b',...     % plot mean
        'LineWidth',1.2)
    hold on; grid on
    plot(SurgeDecay.final.t,SurgeDecay.final.(fnames{k}).x_mean+...         % plot error bars, 95% confidence intervals
        SurgeDecay.final.(fnames{k}).tDistCoeff.*SurgeDecay.final.(fnames{k}).x_std,'--b','LineWidth',1.2)
    plot(SurgeDecay.final.t,SurgeDecay.final.(fnames{k}).x_mean-...    
        SurgeDecay.final.(fnames{k}).tDistCoeff.*SurgeDecay.final.(fnames{k}).x_std,'--b','LineWidth',1.2)
    xlabel('Time (s)');
    ylabel('Normalized Flap Position');
    titStr=strrep(fnames{k},'_','.');
    title(titStr)
end
    
% plot trend in per-unit-mass damping, period with amplitude
figure; 
for k=2:length(fnames)
    subplot(2,1,1)
    scatter(x_vec(k-1),SurgeDecay.final.(fnames{k}).c,[],'b')               % per unit mass damping vs. initial displacement amplitude
    if k==2;
        ylabel('Per-unit-inertia Damping (N-s/kg')
        hold on; grid on;
    end
    subplot(2,1,2)
    scatter(x_vec(k-1),SurgeDecay.final.(fnames{k}).Tn,[],'b')              % natural period vs. initial displacement amplitude
    if k==2;
        ylabel('Damped Period (s)')
        xlabel('\Delta z (m)')
        hold on; grid on
    end
end
    
    
    



        
        