% calculates the period and per-unit mass damping of FOSWEC in pitch free
% decay

close all
clear

% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% directory information
homeDir=pwd;
load_folder = './inter/PitchDecay';
output_variable_name = 'PitchDecay.mat';
output_folder = './final/PitchDecay';

% test information (see test log)
trials = [1:15];                                                            % good trial numbers             
delta_theta = [repmat([2 3 5 7 8.4],1,3)];                                  % initial displacements
dt=0.02;                                                                    % sampling interval
ramp=0.44;                                                                   % retained seconds before release

% User inputs
process_inter=1;                                                            % puts raw .txt files into data structure
process_final=1;                                                            % processes data structure
plot_data=1;                                                                % plots results of processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROCESS INTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if process_inter==1    
    
    % Create output folder if it doesn't exist
    if exist(output_folder) == 0
        mkdir(output_folder)
    end

 % Create structure variable PitchDecay
    for i = 1:length(trials)
        trial_str=sprintf('%02d',trials(i));
        TrialName{i}=strcat('Trial',trial_str);
        PitchDecay.inter.(TrialName{i}).del_theta = delta_theta(i);                  % initial theta displacement
        cd(load_folder)
        cd(['./' TrialName{i}])
        contents = dir( '*.txt' );
        for j = 1:length(contents)
            sensor_name = contents(j).name;
            sensor_name_notxt = strrep(sensor_name,'.txt','');
            PitchDecay.inter.(TrialName{i}).(sensor_name_notxt) = load(sensor_name);  % load data here    
        end
        cd '../../..'
        date_utc = datetime(PitchDecay.inter.(TrialName{i}).time,'convertfrom','datenum');
        PitchDecay.inter.(TrialName{i}).datetime_utc = datetime(date_utc,'TimeZone','UTC');
        PitchDecay.inter.(TrialName{i}).datetime_local = datetime(PitchDecay.inter.(TrialName{i}).datetime_utc,'TimeZone','America/Los_Angeles');
        PitchDecay.inter.(TrialName{i}).time_sec = (PitchDecay.inter.(TrialName{i}).time - PitchDecay.inter.(TrialName{i}).time(1))*60*60*24;
    end
       
cd(homeDir)  
cd(output_folder)    
save(output_variable_name,'PitchDecay') % save structure variable to final folder
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
    
    delta_thetaU = unique(delta_theta);
    cases = length(delta_thetaU);
    
    x_lim = [4600 4700; 3600 3700; 3500 3600;3800 3900;4200 4300;5500 5600; 2600 2700;2600 2800;4100 4200; 3980 4060; 2240 2280;...
               2920 2980; 3800 3900 ; 2600 2650;2900 2950 ];                % bounds sample numbers to the area around the release point
    
    for i = 1:length(trials)
        trial_str=sprintf('%02d',trials(i));
        TrialName{i}=strcat('Trial',trial_str);
        der = gradient(PitchDecay.inter.(TrialName{i}).psPlatPosRy)./dt;    % calculate velocity
        der = filtFB(filtNum,1,der,[],2);                                   % filter velocity
        der2 = gradient(der)./dt;                                           % calculate acceleration
        der2 = filtFB(filtNum,1,der2,[],2);                                 % filter acceleration
        [slope_min,i_slope_min] =  max(der2(x_lim(i,1):x_lim(i,2)));        % this finds the release point of the flap
        
        % correct time stamps relative to release point and mean position
        t_offset(i) = PitchDecay.inter.(TrialName{i}).time_sec(x_lim(i,1) + i_slope_min);
        z_offset(i) = mean(PitchDecay.inter.(TrialName{i}).psPlatPosRy(8450:8550)); % mean final (restored) position, after all motion has settled
        i_win_min = find(PitchDecay.inter.(TrialName{i}).time_sec - t_offset(i) >= -0.5001);
        i_win_min = i_win_min(1);
        i_win_max = find(PitchDecay.inter.(TrialName{i}).time_sec - t_offset(i) <= 8.001);
        i_win_max = i_win_max(end);
        i_window_range(i,:) = [i_win_min i_win_max];                        % window over same relative time period
        t_corrected{i}(:,1) = PitchDecay.inter.(TrialName{i}).time_sec(i_win_min:i_win_max) - t_offset(i); % shift time stamps to considered window
        z_corrected{i}(:,1) = (PitchDecay.inter.(TrialName{i}).psPlatPosRy(i_win_min:i_win_max) - z_offset(i) )... % normalize position
            ./mean(PitchDecay.inter.(TrialName{i}).psPlatPosRy(i_win_min:i_win_min+19) - z_offset(i)); % the average initial position
        
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
    
    PitchDecay.final.t = t_corrected(:,1);                                  % the common time stamps for all processed cases
    
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
    title('Check Plot: Normalized Pitch Decay trials')
    % check plot

 for i = 1:cases
        fname=strcat('Theta',num2str(delta_thetaU(i)));
        fname=strrep(fname,'.','pt');                                       % replace any decimals in field name with 'pt'
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
        
        Cmean(1:(ramp/dt),1)=1;                                             % set all retained values prior to release to 1
        
        % Log normalized data output data structure
        PitchDecay.final.(fname).theta_mean = Cmean;                        % mean position for all trials with given initial displacement
        PitchDecay.final.(fname).theta_std = Cstd;                          % std dev of position for all trials with given initial displacement
        PitchDecay.final.(fname).tDistCoeff = tinvDom(0.95,N(i));           % the 95% confidence interval based on t-distribution of sample size N(i)
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
        PitchDecay.final.(fname).c=c(i);                                    % per unit mass damping
        PitchDecay.final.(fname).zeta=zeta(i);                              % damping ratio
        PitchDecay.final.(fname).Td=Td(i);                                  % damped period
        PitchDecay.final.(fname).Tn=Tn(i);                                  % natural period
        
 end
    cd(homeDir)
    cd(output_folder)
    save(output_variable_name,'PitchDecay') % save structure variable to final folder
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
fnames=fieldnames(PitchDecay.final);                                        % fieldnames of final structure
theta_vec=unique(delta_theta);                                              % unique theta displacement values (dg)

% plot average decay trajectory for each theta displacement
for k=2:length(fnames)                                                      % the first field is time, all others are unique Theta
    figure(k)
    plot(PitchDecay.final.t,PitchDecay.final.(fnames{k}).theta_mean,'b',... % plot mean
        'LineWidth',1.2)
    hold on; grid on
    plot(PitchDecay.final.t,PitchDecay.final.(fnames{k}).theta_mean+...     % plot error bars, 95% confidence intervals
        PitchDecay.final.(fnames{k}).tDistCoeff.*PitchDecay.final.(fnames{k}).theta_std,'--b','LineWidth',1.2)
    plot(PitchDecay.final.t,PitchDecay.final.(fnames{k}).theta_mean-...    
        PitchDecay.final.(fnames{k}).tDistCoeff.*PitchDecay.final.(fnames{k}).theta_std,'--b','LineWidth',1.2)
    xlabel('Time (s)');
    ylabel('Normalized Flap Position');
    title(fnames{k})
end
    
% plot trend in per-unit-mass damping, period with amplitude
figure; 
for k=2:length(fnames)
    subplot(2,1,1)
    scatter(theta_vec(k-1),PitchDecay.final.(fnames{k}).c,[],'b')             % per unit mass damping vs. initial displacement amplitude
    if k==2;
        ylabel('Per-unit-inertia Damping (N-m-s/kg-m^2')
        hold on; grid on;
    end
    subplot(2,1,2)
    scatter(theta_vec(k-1),PitchDecay.final.(fnames{k}).Td,[],'b')            % damped period vs. initial displacement amplitude
    if k==2;
        ylabel('Damped Period (s)')
        xlabel('\Delta \theta (dg)')
        hold on; grid on
    end
end