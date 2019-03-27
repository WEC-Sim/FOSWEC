% calculates spring constant of surge-restoring bungees via FOSWEC static
% offset testing

close all
clear

% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Controls
process_inter=1;
process_final=1;
plotData=1;

% directory information
homeDir=pwd;
load_folder = './inter/SurgeDecay';
load_file = 'SurgeStatic.mat';
output_folder = './final/SurgeStatic';
trials = [21];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if process_inter==1;
    
    % Create output folder if it doesn't exist
    if exist(output_folder) == 0
        mkdir(output_folder)
    end

 % Create structure variable SurgeStatic
    for i = 1:length(trials)                                                % load in data from Trial##.txt data files
        trial_str = sprintf('%02d',trials(i));
        trial_str = ['Trial' trial_str];
        cd(load_folder)
        cd(['./' trial_str])
        contents = dir( '*.txt' );
        for j = 1:length(contents)
            sensor_name = contents(j).name;
            sensor_name_notxt = strrep(sensor_name,'.txt','');
            SurgeStatic.inter.(trial_str).(sensor_name_notxt) = load(sensor_name); % load data to structure
            
        end
        cd '../../..'
        date_utc = datetime(SurgeStatic.inter.(trial_str).time,'convertfrom','datenum');
        SurgeStatic.inter.(trial_str).datetime_utc = datetime(date_utc,'TimeZone','UTC');
        SurgeStatic.inter.(trial_str).datetime_local = datetime(SurgeStatic.inter.(trial_str).datetime_utc,'TimeZone','America/Los_Angeles');
        SurgeStatic.inter.(trial_str).time_sec = (SurgeStatic.inter.(trial_str).time - SurgeStatic.inter.(trial_str).time(1))*60*60*24;
        
    end
    cd(homeDir);
    cd(output_folder);
    save(load_file);
else
    cd(homeDir);
    load(fullfile(output_folder,load_file));
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
    plot(SurgeStatic.inter.(trial_str).platPosx); hold on
    ylabel('Platform Position (m)')
    grid on
    subplot(2,1,2)
    plot(SurgeStatic.inter.(trial_str).lcLadder)
    xlabel('Sample Number')
    ylabel('Applied Force (N)')
    grid on
    % check plot
    
x_equil=mean(SurgeStatic.inter.(trial_str).platPosx([1:1200]));             % equilibrium position        
x_lim = [2200 2600; 3100 3500; 3900 4300; 4900 5300; 5900 6700; 7400 7600; 8000 8200; 8600 8700]; % sample numbers of held positions

  
    for i = 1:length(x_lim)
        SurgeStatic.final.mean_x(i) = median(SurgeStatic.inter.(trial_str)...
            .platPosx(x_lim(i,1):x_lim(i,2)),1)-x_equil;                    % median of steady-state positions
        SurgeStatic.final.mean_F(i) = median(SurgeStatic.inter.(trial_str).lcLadder(x_lim(i,1):x_lim(i,2)),1);  % median of force needed to hold position
        SurgeStatic.final.std_x(i) = std( SurgeStatic.inter.(trial_str).platPosx(x_lim(i,1):x_lim(i,2)),1);   % std dev. of steady-state positions
        SurgeStatic.final.std_F(i) = std( SurgeStatic.inter.(trial_str).lcLadder(x_lim(i,1):x_lim(i,2)),1);    % std dev. of holding force
        
    end
    cd(homeDir);
    cd(output_folder);
    save(load_file);
    
else
    cd(homeDir);
    load(fullfile(output_folder,load_file));
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotData==1;
    
    figure
    plot(-1.*SurgeStatic.final.mean_x,SurgeStatic.final.mean_F,'ko');       % plot force vs. displacement. Flip sign (if needed) so positive displacement -> positive force
    set(gcf, 'Color', 'w');
    grid on
    xlabel('\Delta x [m]')
    ylabel('F [N]')
    box off
    grid on
    title('Static Surge Test')
    p = polyfit(-1.*SurgeStatic.final.mean_x,SurgeStatic.final.mean_F,1);   % linear fit (for linear spring constant determination)
    hold on
    fit = polyval(p,-1.*SurgeStatic.final.mean_x);
    plot(-1.*SurgeStatic.final.mean_x,fit,'k')                              % plot linear fit
    text(0.08,240,['k = ' num2str(p(1))]);
    Rval=corrcoef(fit,SurgeStatic.final.mean_F);                            % correlation coefficient for linear fit
    Rsq=Rval(1,2)^2;
    
    % constrained quadratic fit (constrained so zero force at zero
    % displacement), to capture non-linearities
    cfitcoeff=mrdivide(SurgeStatic.final.mean_F,-1.*SurgeStatic.final.mean_x);
    cfitcoeff=[cfitcoeff, 0];
    cfit=polyval(cfitcoeff,-1.*SurgeStatic.final.mean_x);
    cRval=corrcoef(cfit,SurgeStatic.final.mean_F);
    cRsq=Rval(1,2)^2;                                                       % correlation coefficient for constrained quadratic fit 
end
