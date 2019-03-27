% A script to calculate the hydrostatic coefficient for FOSWEC flap.
close all; clear;

% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User controls
process_inter=1;                                                            % load raw '.txt' data files into '.mat' structure
process_final=1;                                                            % process '.mat' structure
plotData=1;                                                                 % plot processed results

% directory information
homeDir=pwd;                                                                % home directory (facilitates relative paths)
inter_folder='./inter/Flap1Decay';                                          % location of 'Trial##' folders containing '.txt' files
load_file = 'FlapStatic.mat';                                               % '.mat' structure to be created
final_folder = './final/Flap1Static';                                       % directory to which '.mat' structure will be saved
trials = [19]; % this trial number is the static offset run of the decay tests (see test log)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS INTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if process_inter==1;
    % Create structure variable Flap1Static
    for i = 1:length(trials)                                                % load in data from Trial##.txt data files
        trial_str = sprintf('%02d',trials(i));
        trial_str = ['Trial' trial_str];
        cd(inter_folder)
        cd(['./' trial_str])
        contents = dir( '*.txt' );
        for j = 1:length(contents)
            sensor_name = contents(j).name;
            sensor_name_notxt = strrep(sensor_name,'.txt','');
            Flap1Static.inter.(trial_str).(sensor_name_notxt) = load(sensor_name); % load data to structure
            
        end
        cd '../../..'
        date_utc = datetime(Flap1Static.inter.(trial_str).time,'convertfrom','datenum');
        Flap1Static.inter.(trial_str).datetime_utc = datetime(date_utc,'TimeZone','UTC');
        Flap1Static.inter.(trial_str).datetime_local = datetime(Flap1Static.inter.(trial_str).datetime_utc,'TimeZone','America/Los_Angeles');
        Flap1Static.inter.(trial_str).time_sec = (Flap1Static.inter.(trial_str).time - Flap1Static.inter.(trial_str).time(1))*60*60*24;
        
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
    
    R= 0.58;                                                                % flap length (m), for torque calculation
    trial_str = sprintf('%02d',trials);
    trial_str = ['Trial' trial_str];
      
    % check plot
    figure; clf;
    subplot(2,1,1)
    plot(Flap1Static.inter.(trial_str).flapPosF1); hold on
    ylabel('Flap Position (dg)')
    grid on
    subplot(2,1,2)
    plot(Flap1Static.inter.(trial_str).lcLadder)
    xlabel('Sample Number')
    ylabel('Applied Force (N)')
    grid on
    % check plot
    
    % the sample numbers of positions that were "held" were determined via
    % inspection of the check plot
    x_equil=0;                                                              % equilibrium position
    x_lim = [900 970; 1270 1390; 1625 1745; 1990 2080; 2350 2450; 2695 2845; 3075 3175;...
        3400 3525; 3775 3880; 4150 4250];% 4450 4625; 4875 5000; 5275 5500; 5750 5925;...
    % 6170 6300; 6575 6775; 7050 7225; 7450 7650; 7950 8075; 8375 8475]
    
    for i = 1:length(x_lim)
        Flap1Static.final.mean_x(i) = median(Flap1Static.inter.(trial_str)...
            .flapPosF1(x_lim(i,1):x_lim(i,2)),1)-x_equil;                   % median of steady-state displacements
        Flap1Static.final.mean_F(i) = median(Flap1Static.inter.(trial_str).lcLadder(x_lim(i,1):x_lim(i,2)),1)...
            .*R.*cosd(Flap1Static.final.mean_x(i));                                                          % median of torque needed to hold position
        Flap1Static.final.std_x(i) = std( Flap1Static.inter.(trial_str).flapPosF1(x_lim(i,1):x_lim(i,2)),1);    % std dev. of steady-state positions
        Flap1Static.final.std_F(i) = std( Flap1Static.inter.(trial_str).lcLadder(x_lim(i,1):x_lim(i,2)),1);     % std dev. of holding torque
        
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
    plot(-1.*Flap1Static.final.mean_x,Flap1Static.final.mean_F,'ko');
    set(gcf, 'Color', 'w');
    grid on
    xlabel('\Delta \theta [^o]')
    ylabel('\tau [N-m]')
    box off
    grid on
    title('Static Flap Test')
    p = polyfit(-1.*Flap1Static.final.mean_x,Flap1Static.final.mean_F,1);   % linear fit (for linear spring constant determination)
    hold on
    fit = polyval(p,-1.*Flap1Static.final.mean_x);
    plot(-1.*Flap1Static.final.mean_x,fit,'k')                              % plot linear fit
    text(4,15,['k = ' num2str(p(1))])
    Rval=corrcoef(fit,Flap1Static.final.mean_F);                                              
    Rsq=Rval(1,2)^2;                                                        % correlation coefficient for linear fit
    
    % constrained quadratic fit (constrained so zero force at zero
    % displacement), to capture non-linearities
    param=[-1.*(Flap1Static.final.mean_x.^2); -1.*Flap1Static.final.mean_x];        
    cfitcoeff=mrdivide(Flap1Static.final.mean_F,param) ;
    cfitcoeff=[cfitcoeff, 0];
    cfit=polyval(cfitcoeff,-1.*Flap1Static.final.mean_x);
    cRval=corrcoef(cfit,Flap1Static.final.mean_F);
    cRsq=Rval(1,2)^2;                                                       % correlation coefficient for constrained quadratic fit
    
end



