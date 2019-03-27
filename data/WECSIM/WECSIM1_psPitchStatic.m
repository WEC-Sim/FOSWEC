% Dominic Forbush 4.7.2018
% A script to calculate the hydrostatic coefficient for pitch of FOSWEC.

close all; clear;

% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User controls
process_inter=1;                                                            % load raw '.txt' data files into '.mat' structure
process_final=1;                                                            % process '.mat' structure    
plotData=1;                                                                 % plot processed results

% directory
homeDir=pwd;                                                                % home directory (facilitates relative paths)
inter_folder='./inter/PitchDecay';                                          % location of 'Trial ##' folders containing '.txt' files
load_file = 'PitchStatic.mat';                                              % '.mat' structure to be created
final_folder = './final/PitchStatic';                                       % directory to which '.mat' strucutre will be saved
trials = [1]; % this trial number is the static offset run of the decay tests (see test log)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS INTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if process_inter==1;
    % Create structure variable PitchStatic
    
    % Create output folder if it doesn't exist
    if exist(final_folder) == 0
        mkdir(final_folder)
    end
    
    for i = 1:length(trials)                                                % load in data from Trial##.txt data files
        trial_str = sprintf('%02d',trials(i));
        trial_str = ['Trial' trial_str];
        cd(inter_folder)
        cd(['./' trial_str])
        contents = dir( '*.txt' );
        for j = 1:length(contents)
            sensor_name = contents(j).name;
            sensor_name_notxt = strrep(sensor_name,'.txt','');
            PitchStatic.inter.(trial_str).(sensor_name_notxt) = load(sensor_name);         % load data to structure
            
        end
        cd '../../..'
        date_utc = datetime(PitchStatic.inter.(trial_str).time,'convertfrom','datenum');
        PitchStatic.inter.(trial_str).datetime_utc = datetime(date_utc,'TimeZone','UTC');
        PitchStatic.inter.(trial_str).datetime_local = datetime(PitchStatic.inter.(trial_str).datetime_utc,'TimeZone','America/Los_Angeles');
        PitchStatic.inter.(trial_str).time_sec = (PitchStatic.inter.(trial_str).time - PitchStatic.inter.(trial_str).time(1))*60*60*24;
        
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
    
    R= 0.7;                                                                 % platform center to point of load application (m), for torque calculation
    trial_str = sprintf('%02d',trials);
    trial_str = ['Trial' trial_str];
      
    % check plot
    figure; clf;
    subplot(2,1,1)
    plot(PitchStatic.inter.(trial_str).psPlatPosRy); hold on
    ylabel('Platform Position (dg)')
    grid on
    subplot(2,1,2)
    plot(PitchStatic.inter.(trial_str).lcCrane)
    xlabel('Sample Number')
    ylabel('Applied Force (N)')
    grid on
    % check plot
    
    % the positions that were "held" were determined via inspection of the
    % check plot
    x_equil=mean(PitchStatic.inter.(trial_str).psPlatPosRy([1:1200]));      % equilibrium position
    x_lim = [2150 2300; 2660 2900; 3300 3500; 3800 4000; 4300 4500; 4750 4950;...
        5300 5500; 5850 6050; 6350 6500; 6750 6950; 8000 9000];
    
   for i = 1:length(x_lim)
        PitchStatic.final.mean_x(i) = median(PitchStatic.inter.(trial_str)...
            .psPlatPosRy(x_lim(i,1):x_lim(i,2)),1)-x_equil;                 % median of steady-state positions (dg)
        PitchStatic.final.mean_F(i) =R .* median(PitchStatic.inter.(trial_str).lcCrane(x_lim(i,1):x_lim(i,2)),1)...% median of torque needed to hold position
            *cosd(PitchStatic.final.mean_x(i));
        PitchStatic.final.std_x(i) = std( PitchStatic.inter.(trial_str).psPlatPosRy(x_lim(i,1):x_lim(i,2)),1);    % std dev. of steady-state positions
        PitchStatic.final.std_F(i) = std( PitchStatic.inter.(trial_str).lcCrane(x_lim(i,1):x_lim(i,2)),1);        % std dev. of holding torque
        
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
if plotData==1    
    
    figure;                                                                 % plot force vs. displacement. Flip sign (if needed) so positive displacement -> positive force
    plot(-1.*PitchStatic.final.mean_x,PitchStatic.final.mean_F,'ko');
    set(gcf, 'Color', 'w');
    grid on
    xlabel('\Delta \theta [^o]')
    ylabel('\tau [N-m]')
    box off
    grid on
    title('Static Platform Pitch Test')
    p = polyfit(-1.*PitchStatic.final.mean_x,PitchStatic.final.mean_F,1);   % linear fit (for linear spring constant determination)
    hold on
    fit = polyval(p,-1.*PitchStatic.final.mean_x);
    plot(-1.*PitchStatic.final.mean_x,fit,'k')                              % plot linear fit
    text(2,10,['k = ' num2str(p(1))])
    Rval=corrcoef(fit,PitchStatic.final.mean_F);                                              
    Rsq=Rval(1,2)^2;                                                        % correlation coefficient for linear fit
    
    % constrained quadratic fit (constrained so zero force at zero
    % displacement), to capture non-linearities
    param=[-1.*(PitchStatic.final.mean_x.^2); -1.*PitchStatic.final.mean_x];        
    cfitcoeff=mrdivide(PitchStatic.final.mean_F,param) ;
    cfitcoeff=[cfitcoeff, 0];
    cfit=polyval(cfitcoeff,-1.*PitchStatic.final.mean_x);
    cRval=corrcoef(cfit,PitchStatic.final.mean_F);
    cRsq=Rval(1,2)^2;                                                       % correlation coefficient for constrained quadratic fit
    
end

    


