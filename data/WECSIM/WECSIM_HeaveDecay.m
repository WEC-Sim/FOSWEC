% calculates per unit mass damping and period for FOSWEC free decay in
% heave

close all
clear 

% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% directory information
homeDir=pwd;
load_folder = './inter/HeaveDecay';
output_variable_name = 'HeaveDecay.mat';
output_folder = './final/HeaveDecay';

% Test information (see test log)
trials = [1:3 5:16];                                                        % Good trial numbers 
deltaz = [repmat([0.03 0.05 0.07 0.10 0.15],1,3)];                          % initial displacements
dt=0.02;                                                                    % sampling interval
ramp=0.5;                                                                   % retained seconds before release

% User inputs
process_inter=1;                                                            % load raw .txt data files into '.mat' structure
process_final=1;                                                            % process structure data
plot_data=1;                                                                % plot processed results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROCESS INTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if process_inter==1;
    
    % Create output folder if it doesn't exist
    if exist(output_folder) == 0
        mkdir(output_folder)
    end

 % Create structure variable HeaveDecay
    for i = 1:length(trials)
        trial_str=sprintf('%02d',trials(i));
        TrialName{i}=strcat('Trial',trial_str);
        HeaveDecay.inter.(TrialName{i}).del_z = deltaz(i);                  % initial displacement
        cd(load_folder)
        cd(['./' TrialName{i}])
        contents = dir( '*.txt' );
        for j = 1:length(contents)
            sensor_name = contents(j).name;
            sensor_name_notxt = strrep(sensor_name,'.txt','');
            HeaveDecay.inter.(TrialName{i}).(sensor_name_notxt) = load(sensor_name);  % load data here    
        end
        cd '../../..'
        date_utc = datetime(HeaveDecay.inter.(TrialName{i}).time,'convertfrom','datenum');
        HeaveDecay.inter.(TrialName{i}).datetime_utc = datetime(date_utc,'TimeZone','UTC');
        HeaveDecay.inter.(TrialName{i}).datetime_local = datetime(HeaveDecay.inter.(TrialName{i}).datetime_utc,'TimeZone','America/Los_Angeles');
        HeaveDecay.inter.(TrialName{i}).time_sec = (HeaveDecay.inter.(TrialName{i}).time - HeaveDecay.inter.(TrialName{i}).time(1))*60*60*24;
    end
       
cd(homeDir);  
cd(output_folder)    
save(output_variable_name,'HeaveDecay') % save structure variable to final folder
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
    
    deltazU = unique(deltaz);
    cases = length(deltazU);
    
    for i = 1:length(trials)
        trial_str=sprintf('%02d',trials(i));
        TrialName{i}=strcat('Trial',trial_str);
        der = gradient(HeaveDecay.inter.(TrialName{i}).platPosz,HeaveDecay.inter.(TrialName{i}).time_sec);    % calculate flap velocity
        der = filtFB(filtNum,1,der,[],2);                                   % filter velocity
        der2 = gradient(der,HeaveDecay.inter.(TrialName{i}).time_sec);      % calculate acceleration
        der2 = filtFB(filtNum,1,der2,[],2);                                 % filter acceleration
        [slope_min,i_slope_min] =  min(der2);                               % this finds the release point of the platform
        
        % correct time stamps relative to release point and mean position
        t_offset(i) = HeaveDecay.inter.(TrialName{i}).time_sec(i_slope_min);
        z_offset(i) = mean(HeaveDecay.inter.(TrialName{i}).platPosz(i_slope_min+200:i_slope_min+375)); % mean final (restored) position, after all motion has settled
        i_win_min = find(HeaveDecay.inter.(TrialName{i}).time_sec - t_offset(i) >= -0.5001);
        i_win_min = i_win_min(1);
        i_win_max = find(HeaveDecay.inter.(TrialName{i}).time_sec - t_offset(i) <= 3.5001);
        i_win_max = i_win_max(end);
        i_window_range(i,:) = [i_win_min i_win_max];                        % window over same relative time period
        t_corrected{i}(:,1) = HeaveDecay.inter.(TrialName{i}).time_sec(i_win_min:i_win_max) - t_offset(i); % shift time stamps to considered window
        z_corrected{i}(:,1) = (HeaveDecay.inter.(TrialName{i}).platPosz(i_win_min:i_win_max) - z_offset(i) )... % normalize position
            ./mean(HeaveDecay.inter.(TrialName{i}).platPosz(i_win_min:i_win_min+19) - z_offset(i)); % the average initial position
        
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
    
    HeaveDecay.final.t = t_corrected(:,1);                                  % the common time stamps for all processed cases
    
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
        fname=strcat('Z',num2str(deltazU(i)));
        fname=strrep(fname,'.','pt');
        i_trial = trials(i);
        i_delz = find(deltaz== deltazU(i));
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
        HeaveDecay.final.(fname).z_mean = Cmean;                            % mean position for all trials with given initial displacement
        HeaveDecay.final.(fname).z_std = Cstd;                              % std dev of position for all trials with given initial displacement
        HeaveDecay.final.(fname).tDistCoeff = tinvDom(0.95,N(i));           % the 95% confidence interval based on t-distribution of sample size N(i)
        % find pos/neg peaks and calculate per-unit-mass damping using log
        % decrement method
        
        [loc_max maxval]=peakfinder(Cmean, 0, 0, 1, true, false);
        [loc_min minval]=peakfinder(Cmean, 0, 0, -1, true, false);
        
        zeta(i)=1/sqrt(1+(2*pi/(log(maxval(1)/maxval(2))))^2);              % approx for damping ratio, appropriate for oscillating cases
        % (if overdamped or high damping, bad estimate)
        
        if zeta(i) > 0.5                                                    % for high damping, uses fractional overshoot method
            OS=(maxval(1)-maxval(2))/maxval(2);
            zeta(i)=1/sqrt(1+(pi/log(OS))^2);
        end
        Td(i) = loc_max(2)*dt-loc_max(1)*dt-ramp;                           % damped period
        wd(i)= 2*pi/Td(i);                                                  % damped frequency
        wn(i) = wd(i)/sqrt(1-(zeta(i))^2);                                  % natural frequency
        Tn(i)= 2*pi/wn(i);                                                  % natural period
        
        %calculate damping-per-unit-mass based on first log decrement wn and zeta
        c(i) = 2*zeta(i)*wn(i);
        
        % Save damping values and frequency estimates in output data
        % structure
        HeaveDecay.final.(fname).c=c(i);                                    % per unit mass/inertia damping
        HeaveDecay.final.(fname).zeta=zeta(i);                              % damping ratio
        HeaveDecay.final.(fname).Td=Td(i);                                  % damped period
        HeaveDecay.final.(fname).Tn=Tn(i);                                  % natural period
        
    end
    cd(homeDir)
    cd(output_folder)
    save(output_variable_name,'HeaveDecay') % save structure variable to final folder
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
fnames=fieldnames(HeaveDecay.final);                                        % fieldnames of final structure
z_vec=unique(deltaz);                                                       % unique displacement values (m or dg)

% plot average decay trajectory for each z displacement
for k=2:length(fnames)                                                      % the first field is time, all others are unique z
    figure(k)
    plot(HeaveDecay.final.t,HeaveDecay.final.(fnames{k}).z_mean,'b',...     % plot mean
        'LineWidth',1.2)
    hold on; grid on
    plot(HeaveDecay.final.t,HeaveDecay.final.(fnames{k}).z_mean+...         % plot error bars, 95% confidence intervals
        HeaveDecay.final.(fnames{k}).tDistCoeff.*HeaveDecay.final.(fnames{k}).z_std,'--b','LineWidth',1.2)
    plot(HeaveDecay.final.t,HeaveDecay.final.(fnames{k}).z_mean-...    
        HeaveDecay.final.(fnames{k}).tDistCoeff.*HeaveDecay.final.(fnames{k}).z_std,'--b','LineWidth',1.2)
    xlabel('Time (s)');
    ylabel('Normalized Flap Position');
    titStr=strrep(fnames{k},'_','.');
    title(titStr)
end
    
% plot trend in per-unit-mass damping, period with amplitude
figure; 
for k=2:length(fnames)
    subplot(2,1,1)
    scatter(z_vec(k-1),HeaveDecay.final.(fnames{k}).c,[],'b')               % per unit mass damping vs. initial displacement amplitude
    if k==2;
        ylabel('Per-unit-inertia Damping (N-s/kg')
        hold on; grid on;
    end
    subplot(2,1,2)
    scatter(z_vec(k-1),HeaveDecay.final.(fnames{k}).Td,[],'b')              % damped period vs. initial displacement amplitude
    if k==2;
        ylabel('Damped Period (s)')
        xlabel('\Delta z (m)')
        hold on; grid on
    end
end
    
    
    


