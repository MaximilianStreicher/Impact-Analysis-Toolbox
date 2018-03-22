%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date: 11/05/2017 
% Author: Maximilian Streicher, Phd student Ghent University.
% Email: Maximilian.Streicher@UGent.be (str.max1@gmx.de)
%
% Script: Integration of pressure sensors to obtain a force value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%% MANUAL USER INPUT 

clc; clear all; close all;

dir = '.\Output\Pressure\';                                                 % input directory
skiprows = 1;                                                               % rows to skip at the beginning of the input file
skipcol = 1;                                                                % colums to skip at the beginning of input file
del =  '\t';                                                                % delimiter in input file. space: ' ', tab: 't'
fs = 1000;                                                                  % sampling frequency of pressure sensor signals [Hz]
Loc = [0.02 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.12 0.12 0.12 0.12];   % distance between pressure sensors (starting from bottom)

%% LOAD DATA

Sens1 = dlmread([dir 'Sensor_1_filtered.txt'],del,skiprows,skipcol);        % preallocation of matrix
Sens = zeros(length(Sens1),length(Loc));

for i = 1:1:length(Loc)
Sens(:,i) = dlmread([dir 'Sensor_' num2str(i) '_filtered.txt'],...
    '\t',skiprows,skipcol);                                                 % read in data from file
end

time = (0:1/fs:(length(Sens)-1)/fs)';                                       % generate time vector

%% RECTANGULAR INTEGRATION

Int = Sens(:,1).*(Loc(1)+Loc(2)/2);                                         % first sensor above bottom 

for i = 2:length(Loc)-1
    tmp = Sens(:,i).*(Loc(i)/2+Loc(i+1)/2);                                 % integration of pressure sensors 
    Int = Int + tmp;
end

tmp = Sens(:,end).*(Loc(i+1)/2);                                            % upper sensor
Int = (Int + tmp).*100.*1;                                                  % conversion to kN/m (1 Bar = 100 kN/m)
    
%% VISUALISATION

f1 = figure('units','normalized','outerposition',[0 0 1 1]); 
grid minor; hold on;
plot(time,Int,'-k');
set(gca,'FontSize',24);
xlabel('Time [s]'); ylabel('Force [kN/m]'); 
legend('Integrated pressure','FontSize',24);
hold off;

%% SAVE OUTPUT

fileID = fopen([dir '\Integrated_pressures.txt'],'w');                      % save time series after integration
fprintf(fileID,'%s\t%s\n','Time [s]','Integrated Pressure [kN/m]');
fprintf(fileID,'%.4f\t%.4f\n',cat(1,time',Int'));
fclose(fileID);

savefig(f1,[dir 'Integrated_pressures.fig']);                               % save figure
