%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date: 09/05/2017 
% Author: Maximilian Streicher, Phd student Ghent University.
% Email: Maximilian.Streicher@UGent.be (str.max1@gmx.de)
%
% Script: calculate the sum of load cells and force per unit width
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%% MANUAL USER INPUT 

clc; clear all; close all;

dir = '.\Output\Force\';                                                    % input directory
skiprows = 1;                                                               % rows to skip at the beginning of the input file
skipcol = 1;                                                                % colums to skip at the beginning of input file
del =  '\t';                                                                % delimiter in input file. space: ' ', tab: 't'
width = 0.2;                                                                % width of the force unit measurement plate for each load cell [m]
fs = 1000;                                                                  % sample frequency
nr = 2;                                                                     % number of load cells attached to the same measurement plate

%% LOAD DATA

Sens1 = dlmread([dir 'Loadcell_1_filtered.txt'],del,skiprows,skipcol);      % preallocation of matrix
Sens = zeros(length(Sens1),nr);

for i = 1:1:nr
Sens(:,i) = dlmread([dir 'Loadcell_' num2str(i) '_filtered.txt'],...
    '\t',skiprows,skipcol);                                                 % read in data from file
end

time = (0:1/fs:(length(Sens)-1)/fs)';                                       % generate time vector

%% CALCULATION

Sens_sum = sum(Sens,2)./0.2;                                                % calculate the sum of all load cell measurements and convert number to meter unit width [kN/m]    


%% VISUALISATION

f1 = figure('units','normalized','outerposition',[0 0 1 1]); 
grid minor; hold on;
plot(time,Sens_sum,'-k');
set(gca,'FontSize',24);
xlabel('Time [s]'); ylabel('Force [kN/m]'); 
legend('Total force','FontSize',24);
hold off;


%% SAVE OUTPUT

fileID = fopen([dir '\Total_force.txt'],'w');                               % save time series of total impact force
fprintf(fileID,'%s\t%s\n','Time [s]','Total force [kN/m]');
fprintf(fileID,'%.4f\t%.4f\n',cat(1,time',Sens_sum'));
fclose(fileID);

savefig(f1,[dir 'Total_force.fig']);                                        % save figure

