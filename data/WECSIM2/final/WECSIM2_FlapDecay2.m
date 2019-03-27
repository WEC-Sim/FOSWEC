% Determines per-unit-mass damping and natural and damped frequencies for
% the flap using the log decrement method for free-decay tests.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
process_inter = 1;   % 1:processes inter data, 0:loads inter *.mat
process_final = 1;   % 1:processes final data, 0:loads final *.mat
plot_data = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Directory Info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_folder = pwd;
final_folder = './FlapDecay2';
inter_folder = '../inter/FlapDecay2';
log_folder = '../logs';
inter_file = 'FlapDecay2_inter.mat';
final_file = 'FlapDecay2_final.mat';
addpath(genpath(strrep(pwd,'\WECSIM2\final\FlapDecay2','')))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import Test Log and 'inter' Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if process_inter == 1
    % import test log *.xlsx
    cd(log_folder)
    [num,txt,raw] = xlsread('WECSIM2_FlapDecay2.xlsx','Log');
    data.Exp        = 'FlapDecay2';
    data.Header     = txt(5,2:end);
    data.Trial      = num(:,2);
    data.Theta0     = num(:,3);
    data.Notes      = txt(6:end,6);
    data.Flag       = num(:,7);
    clear txt num raw    
    
    FlapDecay2.Theta0   = data.Theta0;    
    FlapDecay2.Flag     = data.Flag;   
    
    % import 'inter' data for all trials          
    numTrials = length(data.Trial);
    cd(inter_folder)        
    for i = 1:numTrials      
        trial_str = sprintf('%02d',data.Trial(i)); 
        Trial = ['Trial' trial_str];    
        file_str = sprintf('%01d',data.Trial(i)); 
        File = ['Trial_' file_str '.mat'];  
        load(File)
        FlapDecay2.inter.(Trial).Theta0     = data.Theta0(i);    
        FlapDecay2.inter.(Trial).Notes      = data.Notes(i);       
        FlapDecay2.inter.(Trial).Flag       = data.Flag(i);  
        FlapDecay2.inter.(Trial).time   	= FlapDecay.time;  
        FlapDecay2.inter.(Trial).flapPosF1 	= FlapDecay.signals(2).values;  
    end
    save(inter_file,'FlapDecay2')                                          % save 'inter' data to *.mat    
    clearvars -except FlapDecay2 final_file final_folder inter_file inter_folder process_final plot_data base_folder
    cd(base_folder)
else 
    load([inter_folder,'\', inter_file])                                   % load 'inter' *.mat    
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
    
    % visual estimate from inter data of oscillation start
    x_lim = [5800 6000;5800 6000;5860 5880; 8150 8350; 8600 8750; 9520 9720; 8350 8550;...
         10550 10750; 9300 9500; 10150 10350; 9300 9500; 9400 9550;... 
         12670 12850; 11200 11400; 9600 9800; 10500 10650; 9900 10750;];
        
    %plot and process inter data
    figure     
    f = gcf;  
    numTrials = length(fieldnames(FlapDecay2.inter));
    for i = 1:numTrials;                                                    %total number of trials (including flagged)  
        ii = i-1;                                                           %raw trial number        
        trial_str = sprintf('%02d',ii); 
        Trial = ['Trial' trial_str]; 
        
        time = FlapDecay2.inter.(Trial).time;
        flapPosF1 = FlapDecay2.inter.(Trial).flapPosF1;
        % filter parameters
        filtNum=0.125*ones(8,1);
        
        % save to final
        if FlapDecay2.inter.(Trial).Flag == 1                         
            FlapDecay2.final.(Trial).Flag = 1; 
            FlapDecay2.final.(Trial).time = FlapDecay2.inter.(Trial).time; % load time here   
            duration = length(FlapDecay2.inter.(Trial).flapPosF1);
            FlapDecay2.final.(Trial).flapPosF1 = zeros(duration,1);        % set data to 0 for flagged cases
        else
            der = gradient(flapPosF1,time);                                % velocity and acceleration calculations
            der = filtFB(filtNum,1,der,[],2);
            der2 = gradient(der, time);
            der2 = filtFB(filtNum,1,der2,[],2);
            [slope_min,i_slope_min] =  max(der2(x_lim(i,1):x_lim(i,2)));                
            t_offset(i) = time(x_lim(i)+i_slope_min);
            theta_ss(i) = mean(flapPosF1(15150:15200));                    % steady-state theta value when oscillation ends
            i_win_min = min(find(time - t_offset(i) >= -0.5001));
            i_win_max = max(find(time - t_offset(i) <= 25.301));
            i_window_range(i,:) = [i_win_min i_win_max];

            % clips to relevant time-period
            time_c(:,i) = time(i_win_min:i_win_max) - t_offset(i);                                 

            % normalizes displacement and sets time<0 and initial
            % displacement equal to one
            theta_norm(:,i) = (flapPosF1(i_win_min:i_win_max) - theta_ss(i) )...
                ./mean(flapPosF1(i_win_min:i_win_min+19) - theta_ss(i));    % normalized flap position           
            dt = 0.01;
            ramp = 0.5;                
            theta_norm(1:ramp/dt,i) = 1;                                    % setting initial norm'd displacement to 1   	            

            % save final (normalized) data
            FlapDecay2.final.(Trial).time = time_c(:,i);
            FlapDecay2.final.(Trial).flapPosF1 = theta_norm(:,i);                    
        end
        
    end   
    
    leg = fieldnames(FlapDecay2.inter);
    legend(leg,'location','northeast')
    title('Flap Decay - Inter Data ')    

    del_theta = FlapDecay2.Theta0;
    theta = nonzeros(unique(FlapDecay2.Theta0));
    for i = 1:length(theta);                                                % number of unique initial displacements
        trialNum = find([del_theta] == theta(i));
        cmean = zeros(i_window_range(1,2) - i_window_range(1,1)+1,1);
        cstd = zeros(i_window_range(1,2) - i_window_range(1,1)+1,1);
        for j = 1:length(trialNum)
            cmean = cmean + theta_norm(:,trialNum(j));                      % summing norm. thetas for each theta setting
            flapPosF1_norm(:,j) = theta_norm(:,trialNum(j));                % normalized flap position
        end
        Cstd(:,i) = std(flapPosF1_norm,0,2);                                % standard deviation between trials common theta
        Cmean(:,i) = cmean./length(trialNum);                               % mean of trials for common theta                        
        
        trialName = fieldnames(FlapDecay2.inter);
        theta_str = sprintf('%01d',theta(i)); 
        Theta = ['Theta_' theta_str]; 
        
        % save calculated values        
        FlapDecay2.FlapDecay.(Theta).trialName = trialName(trialNum);  
        FlapDecay2.FlapDecay.(Theta).time = time_c(:,trialNum(1));
        FlapDecay2.FlapDecay.(Theta).mean = Cmean(:,i); 
        FlapDecay2.FlapDecay.(Theta).stdev = Cstd(:,i);
        
        % find pos/neg peaks and calc damping per-unit-mass using log decrement method
        dt = 0.01;
        ramp = 0.5;
        [loc_max maxval] = peakfinder(FlapDecay2.FlapDecay.(Theta).mean,0,0,...
            1,true,false);  
        loc_max(6:end)=[]; maxval(6:end)=[];
        [loc_min minval] = peakfinder(FlapDecay2.FlapDecay.(Theta).mean,0,0,...
            -1,true,false);
        loc_min(5:end)=[]; minval(5:end)=[];
        for i = 1:length(maxval)-1
            zeta(i) = 1/sqrt(1+(2*pi/(log(maxval(i)/maxval(i+1))))^2);      % damping ratio
            Td(i) = loc_max(i+1)*dt-loc_max(i)*dt-ramp;                     % damped period
            wd(i)= 2*pi/Td(i);                                              % damped frequency 
            wn(i) = wd(i)/sqrt(1-(zeta(i))^2);                              % natural frequency
            Tn(i)= 2*pi/wn(i);                                              % natural period
        end
        
        %calculate damping-per-unit-mass based on first log decrement wn and zeta
        c = 2*zeta(1)*wn(1);
      
        % save calculated values
        FlapDecay2.FlapDecay.(Theta).max_pks = [maxval loc_max];  
        FlapDecay2.FlapDecay.(Theta).min_pks = [minval loc_min];  
        FlapDecay2.FlapDecay.(Theta).zeta = zeta;
        FlapDecay2.FlapDecay.(Theta).Td = Td;
        FlapDecay2.FlapDecay.(Theta).Tn = Tn;
        FlapDecay2.FlapDecay.(Theta).c = c;                 
    end        
    FlapDecay2 = rmfield(FlapDecay2,'inter');
    save(final_file,'FlapDecay2')                    
else
    cd(final_folder)
    load(final_file)   
end       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_data == 1

    theta = nonzeros(unique(FlapDecay2.Theta0));

    %% Plot Final Data by Trial
    
    figure     
        f = gcf;
        numTrials = length(fieldnames(FlapDecay2.final));
        for i = 1:numTrials;                                                %total number of trials (including flagged)       
            ii = i -1;
            trial_str = sprintf('%02d',ii); 
            Trial = ['Trial' trial_str]; 
            plot(FlapDecay2.final.(Trial).time,FlapDecay2.final.(Trial).flapPosF1)       
            hold on       
        end
        set(gcf, 'Color', 'w');
        xlabel('t [s]')
        ylabel('\theta/\theta_o')
        axis([-0.5 25 -.75 1])    
        leg = fieldnames(FlapDecay2.final);
        legend(leg,'location','northeast')         
        box off
        grid on
        title('Flap Decay - Final Data')
   savefig('FlapDecay2_Trials.fig')
    

    %% Plot Trials by initial displacement with error and 
    % fitted exp decay based on log decrement
    
    figure     
        f = gcf;
        f.Position = [1000 918 560 420];
        clear leg
        
        for i = 3:length(theta)
            theta_str = sprintf('%01d',theta(i)); 
            Theta = ['Theta_' theta_str];    
            plot(FlapDecay2.FlapDecay.(Theta).time,FlapDecay2.FlapDecay.(Theta).mean,'LineWidth',1.5)
            hold on
            leg{i} = ['\theta_0 = ' num2str(theta(i),3) '^o'];
        end    

        set(gca,'ColorOrderIndex',1)
        set(gcf, 'Color', 'w');    

        for i = 1:length(theta)
            theta_str = sprintf('%01d',theta(i)); 
            Theta = ['Theta_' theta_str];  

            loc_max = FlapDecay2.FlapDecay.(Theta).max_pks(:,2);
            loc_min = FlapDecay2.FlapDecay.(Theta).min_pks(:,2);
            iter = nonzeros(sort([loc_max;loc_min]));           

            errorbar(FlapDecay2.FlapDecay.(Theta).time(iter),...
               FlapDecay2.FlapDecay.(Theta).mean(iter),...
               tinvDom(0.95,length(FlapDecay2.FlapDecay.(Theta).trialName))...
               *FlapDecay2.FlapDecay.(Theta).stdev(iter),'.','LineWidth',1.5)
            hold on
        end

        grid on
        xlabel('t [s]')
        ylabel('\theta_{norm}')
        axis([-0.5 25 -1 1])
        legend(leg(3:5),'location','northeast')
        box off
        grid on
        title('Flap Decay - Mean Normalized with 95% CI')        
    savefig('FlapDecay2_Theta.fig')        


    %% Plot mean normalized Trials for each trial with error and 
    % fitted exp decay based on log decrement
    
    for i = 1:length(theta)
        figure     
        f = gcf;
        clear leg
    
        %Plot mean norm data                                                % mean normalized flap position for given theta as a function of time
        theta_str = sprintf('%01d',theta(i)); 
        Theta = ['Theta_' theta_str];    
        plot(FlapDecay2.FlapDecay.(Theta).time,FlapDecay2.FlapDecay.(Theta).mean,'k','LineWidth',1.5)
        hold on  
        leg = ['\theta_0 = ' num2str(theta(i),3) '^o'];   

        %Plot error bars
        loc_max = FlapDecay2.FlapDecay.(Theta).max_pks(:,2);
        loc_min = FlapDecay2.FlapDecay.(Theta).min_pks(:,2);
        iter = nonzeros(sort([loc_max;loc_min]));
        
        set(gca,'ColorOrderIndex',1)
        set(gcf, 'Color', 'w');  
        errorbar(FlapDecay2.FlapDecay.(Theta).time(iter),...                % error bar +/- 2 std of the mean
           FlapDecay2.FlapDecay.(Theta).mean(iter),...
           2*FlapDecay2.FlapDecay.(Theta).stdev(iter),'k.');
        hold on

        %Plot exp decay                                                     % the modeled first-order decay envelope, determined from c
        set(gca,'ColorOrderIndex',1)
        set(gcf, 'Color', 'w');     
        plot(FlapDecay2.FlapDecay.(Theta).time,...
            exp(-FlapDecay2.FlapDecay.(Theta).c/2*FlapDecay2.FlapDecay.(Theta).time),'k--');
        hold on

        %Plot trials for same Theta0
        trialNames = FlapDecay2.FlapDecay.(Theta).trialName;
        for ii = 1:length(trialNames)
            trialName = char(trialNames(ii));
            plot(FlapDecay2.final.(trialName).time,FlapDecay2.final.(trialName).flapPosF1)       
            hold on      
        end 

        grid on
        xlabel('t [s]')
        ylabel('\theta_{norm}')
        axis([-0.5 25 -.8 1])
        leg2 = [leg;'95% CI';'Decay';trialNames];
        legend(leg2)
        box off
        grid on
        title(['Flap Decay \theta_0 = ', theta_str,'deg - Trials and Mean Normalized with 95% CI'])        
        savefig(['FlapDecay2_Theta', theta_str,'.fig'])        
    end

%% Plot C by theta    
    figure     
        f = gcf;
        clear leg
        for i = 1:length(theta)                                             % plot damping estimate as a function of theta (non-linearity check)
            theta_str = sprintf('%01d',theta(i)); 
            Theta = ['Theta_' theta_str];    
            c(i) = FlapDecay2.FlapDecay.(Theta).c;
            plot(theta(i),c(i),'sk','LineWidth',1.5)
            hold on
        end
        p = polyfit(theta,c',2);                                            % fit quadratic to (theta,c) data.
        x1 = linspace(0,25);
        y1 = polyval(p,x1);
        plot(x1,y1);

        xlim([0 25])
        ylim([0 0.5])
        title({'mean(c) = ' num2str(mean(c))})
        xlabel('Theta0')
        ylabel('damping per unit mass')
    savefig(['FlapDecay2_damping.fig'])   
    
    cd ..
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
