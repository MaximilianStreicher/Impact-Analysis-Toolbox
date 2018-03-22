%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date: 11/05/2017 
% Author: Maximilian Streicher, Phd student Ghent University.
% Email: Maximilian.Streicher@UGent.be (str.max1@gmx.de)
%
% Script: Peak detection method for integrated pressures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%% MANUAL USER INPUT 

clc; clear all; close all;

dir = '.\Output\Pressure\';                                                 % input directory
file = '#bi_01_5_DR00x15-1kHz.asc';                                         % name of original input file or test name, for naming purposes
skiprows = 1;                                                               % rows to skip at the beginning of the input file
skipcol = 0;                                                                % colums to skip at the beginning of input file
del =  '\t';                                                                % delimiter in input file. space: ' ', tab: 't'
fs = 1000;                                                                  % sampling frequency of pressure sensor signals [Hz]
sel = 0.06;                                                                 % minimum magnitude difference between two consecutive detected peaks [kN/m]
hpt = 0.06;                                                                 % high pass threshold [kN/m]
dt = 2;                                                                     % minimum time between impacts [s] 
addpath('C:\Program Files\Matlab\ResearchR2015b\toolbox\peakfinder');       % add toolboxpath

%% LOAD DATA

Sens = dlmread([dir '\Integrated_pressures.txt'],del,skiprows,skipcol);     % read in integrated pressure data from file


%% PEAK DETECTION

[peakLoc peakMag] = peakfinder(Sens(:,2), sel, hpt, 1, 1);                  % use Matlab peakfinder toolbox to find peaks. 1 stands for maxima and second 1 for including the endpoints

cnt = 1;
X = zeros(length(peakLoc),1);
Y = zeros(length(peakLoc),1);

for i = 1:length(peakLoc) 
    
    window = peakLoc(i)-dt*1000:1:peakLoc(i);                               % filter values which are too close in time
    n = find(ismember(peakLoc,window)==1);
    
    if length(n) == 1                                                       % check if only 1 peak was detected in the range of dt
       X(cnt,1)=peakLoc(n);
       Y(cnt,1)=peakMag(n);
       cnt = cnt+1;
       
    else
       [mx I] = max(peakMag(n));                                            % check backwards if the last entrance in X and Y is part of the same impact and maybe overwrite
    
       if peakMag(n(I)) >= Y(cnt-1)
            X(cnt-1,1)=peakLoc(n(I));
            Y(cnt-1,1)=peakMag(n(I));
        else
        end
    end
    
end

Sens_events = [Sens(X(Y>=hpt),1) Y(Y>=hpt)];                                % high pass threshold


%% VISUALISATION

f1 = figure('units','normalized','outerposition',[0 0 1 1]); 
hold on; grid minor;
plot(Sens(:,1),Sens(:,2),'-','Color',[0 0 0]);
plot(Sens_events(:,1),Sens_events(:,2),'ob','MarkerSize',15);
set(gca,'FontSize',24);
xlabel('Time [s]'); ylabel('Force [kN/m]'); 
legend({'Integrated pressure','Peak event'},...
    'FontSize',24,'Location','NorthWest');


%% SAVE OUTPUT

fileID = fopen([dir 'Events.txt'],'w');                                     % save peak events
fprintf(fileID,'%s\t%s\n','Time [s]','Events [kN/m]');
fprintf(fileID,'%.4f\t%.4f\n',cat(1,Sens_events(:,1)',Sens_events(:,2)'));
fclose(fileID);

fileID = fopen([dir 'Metadata.txt'],'w');                                   % save meta data
fprintf(fileID,'%s\t%s\n','original file:',file);
fprintf(fileID,'%s\t%s\n','high pass threshold [kN]:',num2str(hpt));
fprintf(fileID,'%s\t%s\n','min. difference in magnitude between 2 consecutive impacts [kN/m]:',num2str(sel));
fprintf(fileID,'%s\t%s\n','min. time difference between 2 detected impacts [s]:', num2str(dt));
fprintf(fileID,'%s\t%s\n','force unit:','kN/m');
fclose(fileID);
 
savefig(f1,[dir 'Events.fig']);                                             % save figure


