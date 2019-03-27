% Make power plot for RAO cases for specifed configuration (using flap1 or
% flap 1 and 2, as appropriate). Power is calculated at several points in
% the drivetrain system to estimate component efficiency and bound errors.
% This processing is appropriate for REGULAR waves, for which there are a
% wider range of periods.

% line 225 'SORT' may error for some combinations of configurations and
% Matlab editions. To circumvent error, run one flap at a time (e.g.,
% flap2Check = 1 or flap2Check=2, not flap2Check=[1 2]).

clear; close all; clc;
%% Setup
homeDir=cd;

% select your configuration. Use a single value, chosen from [1:4].
% Config1: only flap 1 unlocked. Config2: flaps 1 and 2 unlocked. 
% Config3: flaps and platform heave unlocked. 
% Config4: flaps and platform heave, surge and pitch are
% unlocked.
configNum=4; 
RAOdir=strcat('./final/Config',num2str(configNum),'Reg');                  % the following processing is appropriate for regular wave cases
filName=strcat('Config',num2str(configNum),'Reg_final.mat');
flap2Check=[1 2];                                                           % flap numbers to check (flap 1, 2)

% mechanical parameters
Ng=4;                                                                       % gear ratio gearbox output shaft to flap shaft (flap moves 1/4 motor speed)
Ngb=71;                                                                     % gearbox ratio: encoder on input side, load cell on output side

% processing parameters
filtNum=0.125.*ones(8,1);                                                   % FIR filter numerator coefficients (100 Hz sampling, 12.5 Hz cut-off)
dt=0.01;                                                                    % sampling rate (s)
sidx=5000;                                                                  % the start index, in samples, of the fully developed oscillation
dur=20;                                                                     % the duration after the stop time to consider, multiples of nominal period (s)

%% load data
f1name=strrep(filName,'_final.mat','');                                  % the first field name of the loaded structure
cd(RAOdir);
data=load(filName);
temp=fieldnames(data.(f1name).inter);
i_min=sidx*ones(length(temp),1);                                            % visual estimate of start times of device response

% process needed intermediate data
for k=1:length(temp)
    if data.(f1name).Flag(k) < 0.5
        
        %% Process flap 1 case
        if length(flap2Check)==1 && abs(flap2Check-1)<0.01
            
            % flap power calculation
            flapVel{k}(:,1)=filtFB(filtNum,1,diff(data.(f1name).inter.(temp{k})... 
                .flapPosF1(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt))./dt,[],2).*(pi/180);
            flapT{k}(:,1)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})...
                .lcTyF1(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            flapP{k}(:,1)=abs(flapT{k}(1:end-1,1).*flapVel{k}(:,1));        % Flap power is assumed agnostic to flap direction
            
            % motor power calculation
            motVel{k}(:,1)=filtFB(filtNum,1,diff(data.(f1name).inter.(temp{k})...
                .motPosF1(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt))./dt,[],2).*(pi/180)./Ngb; % velocity must be stepped through GB to relate to transducer-measured torque
            motT{k}(:,1)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})...
                .ttTrqF1(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            motP{k}(:,1)=abs(motT{k}(1:end-1,1).*motVel{k}(:,1));          % Power is assumed agnostic to direction
            
            % motor command power calculation
            motCmdT{k}(:,1)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})...
                .motCmdTrqF1(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            motCmdP{k}(:,1)=abs(motCmdT{k}(1:end-1,1).*motVel{k}(:,1))*Ngb;% Power is assumed agnostic to direction
            
            % motor electrical power calculation
            elecV{k}(:,1)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})... % supply voltage
                .motVPwrSupF1(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            elecI{k}(:,1)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})... % supply current
                .motILoadF1(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            elecP{k}(:,1)=abs(elecV{k}(:,1).*elecI{k}(:,1));                % Power is assumed agnostic to direction
            
            legendVec={'P_{flap,1}','P_{motor,1}','P_{cmd,1}'};
            
            %% Process flap 2 case
        elseif length(flap2Check)==1 &&  abs(flap2Check-2)<0.01
            % flap power
            flapVel{k}(:,1)=filtFB(filtNum,1,diff(data.(f1name).inter.(temp{k})...
                .flapPosF2(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt))./dt,[],2).*(pi/180);
            flapT{k}(:,1)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})...
                .lcTyF2(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            flapP{k}(:,1)=abs(flapT{k}(1:end-1,1).*flapVel{k}(:,1));        % Power is assumed agnostic to direction
            
            % motor power
            motVel{k}(:,1)=filtFB(filtNum,1,diff(data.(f1name).inter.(temp{k})...
                .motPosF2(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt))./dt,[],2).*(pi/180)./Ngb; % velocity must be stepped through GB to relate to transducer-measured torque
            motT{k}(:,1)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})...
                .ttTrqF2(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            motP{k}(:,1)=abs(motT{k}(1:end-1,1).*motVel{k}(:,1));           % Power is assumed agnostic to direction
            
            % motor command power
            motCmdT{k}(:,1)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})...
                .motCmdTrqF2(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            motCmdP{k}(:,1)=abs(motCmdT{k}(1:end-1,1).*motVel{k}(:,1))*Ngb; % Power is assumed agnostic to direction
            
            % motor electrical power
            elecV{k}(:,1)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})... % supply voltage
                .motVPwrSupF2(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            elecI{k}(:,1)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})... % supply current
                .motILoadF2(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            elecP{k}(:,1)=abs(elecV{k}(:,1).*elecI{k}(:,1));                % Power is assumed agnostic to direction
            
            legendVec={'P_{flap,2}','P_{motor,2}','P_{cmd,2}'};
            
            %% Process both flaps
        else % implies both flaps active
            % Flap 1
            % initialize figures
            % flap power
            flapVel{k}(:,1)=filtFB(filtNum,1,diff(data.(f1name).inter.(temp{k})...
                .flapPosF1(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt))./dt,[],2).*(pi/180);
            flapT{k}(:,1)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})...
                .lcTyF1(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            flapP{k}(:,1)=abs(flapT{k}(1:end-1,1).*flapVel{k}(:,1));        % Power is assumed agnostic to direction
            
            % motor power
            motVel{k}(:,1)=filtFB(filtNum,1,diff(data.(f1name).inter.(temp{k})...
                .motPosF1(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt))./dt,[],2).*(pi/180)./Ngb; % velocity must be stepped through GB to relate to transducer-measured torque
            motT{k}(:,1)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})...
                .ttTrqF1(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            motP{k}(:,1)=abs(motT{k}(1:end-1,1).*motVel{k}(:,1));           % Power is assumed agnostic to direction
            
            % motor command power
            motCmdT{k}(:,1)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})...
                .motCmdTrqF1(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            motCmdP{k}(:,1)=abs(motCmdT{k}(1:end-1,1).*motVel{k}(:,1))*Ngb; % Power is assumed agnostic to direction
            
            % motor electrical power
            elecV{k}(:,1)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})... % supply voltage
                .motVPwrSupF1(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            elecI{k}(:,1)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})... % supply current
                .motILoadF1(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt)./dt,[],2);
            elecP{k}(:,1)=abs(elecV{k}(:,1).*elecI{k}(:,1));                % Power is assumed agnostic to direction
            
            % Flap 2
            % flap power
            flapVel{k}(:,2)=filtFB(filtNum,1,diff(data.(f1name).inter.(temp{k})...
                .flapPosF2(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt).*(pi/180))./dt,[],2);
            flapT{k}(:,2)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})...
                .lcTyF2(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            flapP{k}(:,2)=abs(flapT{k}(1:end-1,2).*flapVel{k}(:,2));        % Power is assumed agnostic to direction
            
            % motor power
            motVel{k}(:,2)=filtFB(filtNum,1,diff(data.(f1name).inter.(temp{k})...
                .motPosF2(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt))./dt,[],2).*(pi/180)./Ngb; % velocity must be stepped through GB to relate to transducer-measured torque
            motT{k}(:,2)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})...
                .ttTrqF2(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            motP{k}(:,2)=abs(motT{k}(1:end-1,2).*motVel{k}(:,2));           % Power is assumed agnostic to direction
            
            % motor command power
            motCmdT{k}(:,2)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})...
                .motCmdTrqF2(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            motCmdP{k}(:,2)=abs(motCmdT{k}(1:end-1,2).*motVel{k}(:,2))*Ngb;
            
            % motor electrical power
            elecV{k}(:,2)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})... % supply voltage
                .motVPwrSupF2(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            elecI{k}(:,2)=filtFB(filtNum,1,data.(f1name).inter.(temp{k})... % supply current
                .motILoadF2(i_min(k):i_min(k)+ dur*data.(f1name).T(k)./dt),[],2);
            elecP{k}(:,2)=abs(elecV{k}(:,2).*elecI{k}(:,2));                % Power is assumed agnostic to direction
            
            % pre-format plot legend
            legendVec={'P_{flap,1}','P_{motor,1}','P_{cmd,1}','P_{flap,2}','P_{motor,2}','P_{cmd,2}'};
            
        end
    else % catch case flagged errors
        if length(flap2Check)==2
            flapVel{k}=[NaN,NaN; NaN, NaN];
            flapT{k}=[NaN,NaN; NaN, NaN];
            flapP{k}=[NaN,NaN; NaN, NaN];
            motVel{k}=[NaN,NaN; NaN, NaN];
            motT{k}=[NaN,NaN; NaN, NaN];
            motP{k}=[NaN,NaN; NaN, NaN];
            motCmdT{k}=[NaN,NaN; NaN, NaN];
            motCmdP{k}=[NaN,NaN; NaN, NaN];
            elecV{k}=[NaN,NaN; NaN, NaN];
            elecI{k}=[NaN,NaN; NaN, NaN];
            elecP{k}=[NaN,NaN; NaN, NaN];
        else
            flapVel{k}=NaN;
            flapT{k}=NaN;
            flapP{k}=NaN;
            motVel{k}=NaN;
            motT{k}=NaN;
            motP{k}=NaN;
            motCmdT{k}=NaN;
            motCmdP{k}=NaN;
            elecV{k}=NaN;
            elecI{k}=NaN;
            elecP{k}=NaN;
        end
    end
end

%% Plotting functions
[Hunq,ia,ic]=unique(data.(f1name).H);

% initialize figures and plotting functions
if length(flap2Check)==1 && flap2Check==1                                   % formats title vectors based on flap selection
    titlePre='Flap 1';
elseif length(flap2Check)==1 && flap2Check==2
    titlePre='Flap 2';
else
    titlePre='Flap 1+2';
end

for k=1:length(Hunq);                                                       % a separate figure will be made for each wave height
    figure(k); clf
    goodidx=find(ic==k);
    titleSuff=strcat('H=',num2str(data.(f1name).H(goodidx(1))));
    for k2=1:length(goodidx)                                                % for all indices at selected wave height
        flapPmed(k2,:)=median(flapP{goodidx(k2)});                          % take median by column
        motPmed(k2,:)=median(motP{goodidx(k2)});                            
        motCmdPmed(k2,:)=median(motCmdP{goodidx(k2)});
        elecPmed(k2,:)=median(elecP{goodidx(k2)});
    end
    [Tsort,sortIdx]=sort(data.(f1name).T(goodidx))                          % sort by increasing wave period
    plot(Tsort,flapPmed(sortIdx,1),'-ob','MarkerFaceColor','b')             % plot power measured at flap
    hold on; grid on;
    plot(Tsort,motPmed(sortIdx,1),'-or','MarkerFaceColor','r')              % plot power measured at motor
    plot(Tsort,motCmdPmed(sortIdx,1),'-og','MarkerFaceColor','g')           % plot power command to motor

    if length(flap2Check)==2                                                % overlay flap 2 data if meaningful 
        [Tsort2,sortidx2]=sort(data.(f1name).T(goodidx));
        plot(Tsort2,flapPmed(sortidx2,2),'--^b','MarkerFaceColor','b') % plot power measured at flap
        plot(Tsort2,motPmed(sortidx2,2),'--^r','MarkerFaceColor','r')  % plot power measured at motor
        plot(Tsort2,motCmdPmed(sortidx2,2),'--^g','MarkerFaceColor','g')  % plot power command to motor
    end
    legend(legendVec);
    title({titlePre titleSuff});
    xlabel('Wave Period (s)')
    ylabel('Power (W)')
    clear flapPmed motPmed motCmdPmed elecPmed sortidx sortidx2 Tsort Tsort2
end