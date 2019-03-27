% Dominic Forbush 4.7.2018
% A script to calculate the hydrostatic coefficient for heave of FOSWEC.

close all; clear;

% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User controls
process_inter=1;                                                            % load raw '.txt' data files into '.mat' structure
process_final=1;                                                            % process '.mat' structure
plotData=1;                                                                 % plot processed results

% directory information
homeDir=pwd;                                                                % home directory (facilitates relative paths)
inter_folder='./inter/HeaveDecay';                                          % location of 'Trial ##' folders containing '.txt' files
load_file = 'HeaveStatic.mat';                                              % '.mat' structure to be created
final_folder = './final/HeaveStatic';                                       % directory to which '.mat' strucutre will be saved
trials = [17]; % this trial number is the static offset run of the decay tests (see test log)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS INTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if process_inter==1;
   
    % Create output folder if it doesn't exist
    if exist(final_folder) == 0
        mkdir(final_folder)
    end
    
    % Create structure variable HeaveStatic
    for i = 1:length(trials)                                                % load in data from Trial##.txt data files
        trial_str = sprintf('%02d',trials(i));
        trial_str = ['Trial' trial_str];
        cd(inter_folder)
        cd(['./' trial_str])
        contents = dir( '*.txt' );
        for j = 1:length(contents)
            sensor_name = contents(j).name;
            sensor_name_notxt = strrep(sensor_name,'.txt','');
            HeaveStatic.inter.(trial_str).(sensor_name_notxt) = load(sensor_name); % load data to structure
            
        end
        cd '../../..'
        date_utc = datetime(HeaveStatic.inter.(trial_str).time,'convertfrom','datenum');
        HeaveStatic.inter.(trial_str).datetime_utc = datetime(date_utc,'TimeZone','UTC');
        HeaveStatic.inter.(trial_str).datetime_local = datetime(HeaveStatic.inter.(trial_str).datetime_utc,'TimeZone','America/Los_Angeles');
        HeaveStatic.inter.(trial_str).time_sec = (HeaveStatic.inter.(trial_str).time - HeaveStatic.inter.(trial_str).time(1))*60*60*24;
        
    end
    cd(homeDir);
    cd(final_folder);
    save(load_file);
else
    cd(homeDir);
    load(fullfile(final_folder,load_file));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS FINAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if process_final==1;
   
    trial_str = sprintf('%02d',trials);
    trial_str = ['Trial' trial_str]; 
    
    % check plot
    figure; clf;
    subplot(2,1,1)
    plot(HeaveStatic.inter.(trial_str).platPosz); hold on
    ylabel('Platform Position (m)')
    grid on
    subplot(2,1,2)
    plot(HeaveStatic.inter.(trial_str).lcCrane)
    xlabel('Sample Number')
    ylabel('Applied Force (N)')
    grid on
    % check plot
    
    % the sample numbers of positions that were "held" were determined via inspection
    % of the check plot
    x_equil=mean(HeaveStatic.inter.(trial_str).platPosz([1:1200]));         % equilibrium position
    x_lim = [2200 2600; 3100 3500; 3900 4300; 4900 5300; 5900 6700; 7400 7600; 8000 8200; 8600 8700];
    
    for i = 1:length(x_lim)
        HeaveStatic.final.mean_z(i) = median(HeaveStatic.inter.(trial_str)...
            .platPosz(x_lim(i,1):x_lim(i,2)),1)-x_equil;                    % median of steady-state displacements
        HeaveStatic.final.mean_F(i) = median(HeaveStatic.inter.(trial_str).lcCrane(x_lim(i,1):x_lim(i,2)),1);  % median of force needed to hold position
        HeaveStatic.final.std_z(i) = std( HeaveStatic.inter.(trial_str).platPosz(x_lim(i,1):x_lim(i,2)),1);   % std dev. of steady-state positions
        HeaveStatic.final.std_F(i) = std( HeaveStatic.inter.(trial_str).lcCrane(x_lim(i,1):x_lim(i,2)),1);    % std dev. of holding force
        
    end
    cd(homeDir);
    cd(final_folder);
    save(load_file);
    
else
    cd(homeDir);
    load(fullfile(final_folder,load_file));
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotData==1;
    
    figure
    plot(HeaveStatic.final.mean_z,HeaveStatic.final.mean_F,'ko');           % plot force vs. displacement. Flip Load cell (if needed) so positive displacement -> positive force
    set(gcf, 'Color', 'w');
    grid on
    xlabel('\Delta z [m]')
    ylabel('F [N]')
    box off
    grid on
    title('Static Heave Test')
    p = polyfit(HeaveStatic.final.mean_z,HeaveStatic.final.mean_F,1);       % linear fit (for linear spring constant determination)
    hold on
    fit = polyval(p,HeaveStatic.final.mean_z);
    plot(HeaveStatic.final.mean_z,fit,'k')                                  % plot linear fit
    text(0.08,240,['k = ' num2str(p(1))])
    Rval=corrcoef(fit,HeaveStatic.final.mean_F);                            % correlation coefficient for linear fit
    Rsq=Rval(1,2)^2;
    
    % constrained quadratic fit (constrained so zero force at zero
    % displacement), to capture non-linearities
    cfitcoeff=mrdivide(HeaveStatic.final.mean_F,HeaveStatic.final.mean_z);
    cfitcoeff=[cfitcoeff, 0];
    cfit=polyval(cfitcoeff,HeaveStatic.final.mean_z);
    cRval=corrcoef(cfit,HeaveStatic.final.mean_F);
    cRsq=Rval(1,2)^2;                                                       % correlation coefficient for constrained quadratic fit 
end




