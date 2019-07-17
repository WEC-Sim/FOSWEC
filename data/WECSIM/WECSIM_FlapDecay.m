% calculates per-unit-mass damping and period for free decay of the FOSWEC
% flap

close all
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
process_inter = 0;   % 1:processes inter data, 0:loads inter *.mat
process_final =0;   % 1:processes final data, 0:loads final *.mat
plot_data=1;        % plot processed results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Directory Info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
homeDir=pwd;
inter_folder = './inter/FlapDecay1';
final_folder = './final/FlapDecay1';
final_file = 'FlapDecay1.mat';

% Test information (see test log)
trials = [1:5 7:17];                                                        % successful trial numbers (excluding static offset case)
delta_theta = [3 repmat([5 7 10 15 20],1,3)];                               % nominal initial theta displacement (degrees)
theta_i = [0.91 0.4 1.4 2.4 1.5 0.71 0.8 0.7 0.9 0.9 0.9 0.7...             % initial flap angle before displacement (degrees). Should be nominally zero.
    0.9 0.7 0.7 0.7];
dt=0.02;                                                                    % Sampling interval of FOSWEC mounted sensors
ramp=0.5;                                                                   % by inspection, the time prior to t=0 (in corrected time series) when flap motion begins

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROCESS INTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if process_inter==1;
    
    % Create output folder if it doesn't exist
    if exist(final_folder) == 0
        mkdir(final_folder)
    end
    
    % Create structure vairable Flap1Decay
    for i = 1:length(trials)
        trial_str=sprintf('%02d',trials(i));
        TrialName{i}=strcat('Trial',trial_str);
        Flap1Decay.inter.(TrialName{i}).del_theta = delta_theta(i);         % initial theta displacement
        Flap1Decay.inter.(TrialName{i}).theta_i = theta_i(i);               % steady-state theta (after decay motion has dissipated)
        cd(inter_folder)
        cd(['./' TrialName{i}])
        contents = dir( '*.txt' );
        for j = 1:length(contents)
            sensor_name = contents(j).name;
            sensor_name_notxt = strrep(sensor_name,'.txt','');
            Flap1Decay.inter.(TrialName{i}).(sensor_name_notxt) = load(sensor_name);  % load data here
            
        end
        cd '../../..'
        date_utc = datetime(Flap1Decay.inter.(TrialName{i}).time,'convertfrom','datenum');
        Flap1Decay.inter.(TrialName{i}).datetime_utc = datetime(date_utc,'TimeZone','UTC');
        Flap1Decay.inter.(TrialName{i}).datetime_local = datetime(Flap1Decay.inter.(TrialName{i}).datetime_utc,'TimeZone','America/Los_Angeles');
        Flap1Decay.inter.(TrialName{i}).time_sec = (Flap1Decay.inter.(TrialName{i}).time - Flap1Decay.inter.(TrialName{i}).time(1))*60*60*24;
    end    
    
    cd(final_folder)
    save(final_file,'Flap1Decay') % save structure variable to final folder
    cd '../..'
    
else
    cd(homeDir)
    cd(final_folder);
    load(final_file);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROCESS FINAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if process_final==1;
    
    filtNum=0.125*ones(8,1);                                                  % FIR filter coefficients, a 12.5 Hz cut-off for 100 Hz sampling
    
    delta_thetaU = unique(delta_theta);
    cases = length(delta_thetaU);
    
    x_lim = [1500 1700; 2075 2150; 1850 1950; 2650 2750; 2600 2650; 2510 2550;... % bounds sample numbers to the area around the release point
        1000 1100; 1700 1850; 1800 1900; 1900 2050; 1150 1250; 1800 1900; ...
        1450 1550; 1500 1600; 1700 1900; 1200 1300 ];
    
    
    for i = 1:length(trials)
        trial_str=sprintf('%02d',trials(i));
        TrialName{i}=strcat('Trial',trial_str);
        der = gradient(Flap1Decay.inter.(TrialName{i}).flapPosF1,Flap1Decay.inter.(TrialName{i}).time_sec);    % calculate flap velocity
        der = filtFB(filtNum,1,der,[],2);                                   % filter flap velocity
        der2 = gradient(der,Flap1Decay.inter.(TrialName{i}).time_sec);      % calculate flap acceleration
        der2 = filtFB(filtNum,1,der2,[],2);                                 % filter flap acceleration
        [slope_min,i_slope_min] =  max(der2(x_lim(i,1):x_lim(i,2)));        % this finds the release point of the flap
        
        % correct time stamps relative to release point and mean position
        t_offset(i) = Flap1Decay.inter.(TrialName{i}).time_sec(x_lim(i)+i_slope_min);
        z_offset(i) = mean(Flap1Decay.inter.(TrialName{i}).flapPosF1(2940:2950));             % mean final (restored) position, after all motion has settled
        i_win_min = find(Flap1Decay.inter.(TrialName{i}).time_sec - t_offset(i) >= -0.5001);
        i_win_min = i_win_min(1);
        i_win_max = find(Flap1Decay.inter.(TrialName{i}).time_sec - t_offset(i) <= 6.301);
        i_win_max = i_win_max(end);
        i_window_range(i,:) = [i_win_min i_win_max];                       % window over same relative time period
        t_corrected{i}(:,1) = Flap1Decay.inter.(TrialName{i}).time_sec(i_win_min:i_win_max) - t_offset(i); % shift time stamps to considered window
        z_corrected{i}(:,1) = (Flap1Decay.inter.(TrialName{i}).flapPosF1(i_win_min:i_win_max) - z_offset(i) )... % normalize flap position
            ./mean(Flap1Decay.inter.(TrialName{i}).flapPosF1(i_win_min:i_win_min+19) - z_offset(i)); % the average initial position
        
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
    
    Flap1Decay.final.t = t_corrected(:,1);                                  % the common time stamps for all processed cases
    
    % check plot; plots each trial
    figure; clf;
    for i = 1:length(trials)
        plot(t_corrected(:,i),z_corrected(:,i)), hold on
        leg{i} = ['Trial ' num2str(trials(i),2)];
    end
    set(gcf, 'Color', 'w');
    grid on
    xlabel('t [s]')
    ylabel('\theta/\theta_o')
    axis([-0.5 8 -.75 1])
    legend(leg,'location','northeast')
    box off
    grid on
    title('Check Plot: Normalized Flap Decay trials')
    % check plot
    
    for i = 1:cases
        fname=strcat('Theta',num2str(delta_thetaU(i)));
        i_trial = trials(i);
        i_delz = find(delta_theta== delta_thetaU(i));
        cmean = zeros(i_window_range(1,2) - i_window_range(1,1),1);         % initialize mean position matrix
        cstd = zeros(i_window_range(1,2) - i_window_range(1,1),1);          % initialize postion std. dev matrix
        N(i) = length(i_delz);
        
        for j = 1:N(i)                                                      % for each repeated trial of given del_theta, find average position, std of positions
            c3(:,j) = z_corrected(:,i_delz(j));
        end
        Cstd(:,1) = std(c3,0,2);                                            % measure of run-to-run variance
        Cmean(:,1) = mean(c3,2);                                            % mean of run-to-run position
        
        Cmean(1:(ramp/dt),1)=1;                                             % set all retained values prior to release of flap to 1
        
        
        % Log normalized data output data structure
        Flap1Decay.final.(fname).theta_mean = Cmean;                        % mean position for all trials with given initial displacement
        Flap1Decay.final.(fname).theta_std = Cstd;                          % std dev of position for all trials with given initial displacement
        Flap1Decay.final.(fname).tDistCoeff = tinvDom(0.95,N(i));           % the 95% confidence interval based on t-distribution of sample size N(i)
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
        Flap1Decay.final.(fname).c=c(i);                                    % per unit mass damping
        Flap1Decay.final.(fname).zeta=zeta(i);                              % damping ratio
        Flap1Decay.final.(fname).Td=Td(i);                                  % damped period
        Flap1Decay.final.(fname).Tn=Tn(i);                                  % natural period
        
    end
    cd(homeDir)
    cd(final_folder)
    save(final_file,'Flap1Decay') % save structure variable to final folder
    cd '../..'
    
else
    cd(homeDir)
    cd(final_folder)
    load(final_file);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
fnames=fieldnames(Flap1Decay.final);                                        % fieldnames of final structure
theta_vec=unique(delta_theta);                                              % unique theta displacement values (dg)

% plot average decay trajectory for each theta displacement
for k=2:length(fnames)                                                      % the first field is time, all others are unique Theta
    figure(k)
    plot(Flap1Decay.final.t,Flap1Decay.final.(fnames{k}).theta_mean,'b',... % plot mean
        'LineWidth',1.2)
    hold on; grid on
    plot(Flap1Decay.final.t,Flap1Decay.final.(fnames{k}).theta_mean+...     % plot error bars, 95% confidence intervals
        Flap1Decay.final.(fnames{k}).tDistCoeff.*Flap1Decay.final.(fnames{k}).theta_std,'--b','LineWidth',1.2)
    plot(Flap1Decay.final.t,Flap1Decay.final.(fnames{k}).theta_mean-...    
        Flap1Decay.final.(fnames{k}).tDistCoeff.*Flap1Decay.final.(fnames{k}).theta_std,'--b','LineWidth',1.2)
    xlabel('Time (s)');
    ylabel('Normalized Flap Position');
    title(fnames{k})
end
    
% plot trend in per-unit-mass damping, period with amplitude
figure; 
for k=2:length(fnames)
    subplot(2,1,1)
    scatter(theta_vec(k-1),Flap1Decay.final.(fnames{k}).c,[],'b')             % per unit mass damping vs. initial displacement amplitude
    if k==2;
        ylabel('Per-unit-inertia Damping (N-m-s/kg-m^2')
        hold on; grid on;
    end
    subplot(2,1,2)
    scatter(theta_vec(k-1),Flap1Decay.final.(fnames{k}).Td,[],'b')            % damped period vs. initial displacement amplitude
    if k==2;
        ylabel('Damped Period (s)')
        xlabel('\Delta \theta (dg)')
        hold on; grid on
    end
end
    
    
    
    
    
    
    