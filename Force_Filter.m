%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date: 11/05/2017 
% Author: Maximilian Streicher, Phd student Ghent University.
% Email: Maximilian.Streicher@UGent.be (str.max1@gmx.de)
%
% Script: Filtering of force signals obtained with load cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%% MANUAL USER INPUT 

clc; clear all; close all;

dir = '.\Input\Force\';                                                     % input directory
file = '#bi_01_5_Krm2x-1kHz.asc';                                           % input .txt file with load cell signals (nr. of colums equals the number of load cells)
skiprows = 1;                                                               % rows to skip at the beginning of the input file
skipcol = 0;                                                                % colums to skip at the beginning of input file
del =  ' ';                                                                 % delimiter in input file. space: ' ', tab: 't'
addpath('.\wafo'); initwafo('full');                                        % add path to wafo toolbox and initialize (http://www.maths.lth.se/matstat/wafo/)
fs = 1000;                                                                  % sampling frequency of load cell signals [Hz]
hpf = 0;                                                                    % high pass filter [Hz]. Set to 0 if no high pass filter is used
lpf = 25;                                                                   % low pass filter [Hz]. Set to 0 if no low pass filter is used
lbpf = [1.32 4.01];                                                         % lower boundaries of bandpass filters [Hz]. Set to 0 if no bandwidth filter is used
ubpf = [1.37 4.12];                                                         % upper boundaries of bandpass filters [Hz]. Set to 0 if no bandwidth filter is used
noise = [0.02 -0.01];                                                          % upper noise boundary for each load cell signal. All values below the upper noise boundary will be used to calculate the polynomial fit through the data in order to do the zero correction of the signal


%% LOAD DATA 

in = dlmread([dir file],del,skiprows,skipcol);                              % read in data from file
time = (0:1/fs:(length(in)-1)/fs)';                                         % generate time vector


%% FILTERING


for i = 1:size(in,2)                                                        % go over each load cell individually
    
% REMOVE DRIFT WITH DETREND   
detrY = detrend(in(:,i));                                                   % removes the best straight-line fit from the signal

% REMOVE DRIFT WITH POLYNOMIAL BEST-FIT
idx = find(detrY<=noise(i));                                                % find values below the upper noise level to generat polynomial best-fit line, in order to correct for non-linear drift in the signal
yint = interp1(time(idx),detrY(idx),time);
[p,s,mu] = polyfit(time(idx),detrY(idx),20);                                % fit polynomial best-fit line through selected data 
f_y = polyval(p,time,[],mu);                                                % calculate values for best-fit line
noDrift = detrY-f_y;                                                        % remove the polynomia best-fit line from the signal 

% COMPUTE SPECTRUM FOR COMPARISON
Sest = dat2spec(cat(2,time,noDrift),length(noDrift)/2,'f');                 % compute spectrum using Wafo toolbox for visualisation purposes
[P,F] = pwelch(noDrift,fs*100,(fs/2)*100,fs*100,fs);                        % compute spectrum using pwelch function in Matlab for visualisation purposes

% FILTERING IN FREQUENCY DOMAIN
X_mags = fft(noDrift);                                                      % Fourier transformation
num_bins = length(X_mags);                                                  % number of bins

bin_vals = 0 : num_bins-1;                                                  % x-vector
fax_Hz = (bin_vals*fs/num_bins);                                            % convert bins into frequencies
N = floor(num_bins/2);                                                      % half of the spectrum only

psdx = (1/(fs*num_bins)) * abs(X_mags).^2;                                  % caculate power spectral density
NY = fs/2;                                                                  % nyquist frequency, half the sampling frequency


% FILTERING
FIL = noDrift;
      
% low pass filter
if lpf ~= 0
     [b a] = butter(10,lpf/NY,'low');                                       % low pass butterworth filter design
     R1 = freqz(b,a,N);                                                     % frequency response butterworth filter
     FIL = filter(b,a,FIL);                                                 % filtering
     clear a b
end

% high pass filter
if hpf ~= 0
     [b a] = butter(2,hpf/NY,'high');                                       % high pass butterworth filter design
     R2 = freqz(b,a,N);                                                     % frequency response butterworth filter
     FIL = filter(b,a,FIL);                                                 % filtering
     clear a b
end

% bandwidth filter
R3 = [];
if ubpf ~= 0
    for cnt = 1:length(ubpf)
    [b a] = butter(1,[lbpf(cnt) ubpf(cnt)]/NY,'stop');                      % bandwidth butterworth filter design
    tmp = freqz(b,a,N);                                                     % frequency response butterworth filter
    FIL = filter(b,a,FIL);                                                  % filtering
    R3 = cat(2,R3,tmp);
    clear a b tmp
    end
end


%% VISUALIZE OUTPUT

% plot of time series from measured pressure sensor signal
f1 = figure('units','normalized','outerposition',[0 0 1 1]); 
hold on; grid on; 
plot(time,in(:,i),'-k');                                                    % plot raw data signal
set(gca,'FontSize',24);
xlabel('Time [s]');ylabel('Force [kN]');
title('Raw data: Force signal');
legend('Time series measured force','FontSize',24);
hold off;

% plot of frequency spectrum
f2 = figure('units','normalized','outerposition',[0 0 1 1]); 
hold on; grid on;
leg(1) = plot(fax_Hz(1:N),psdx(1:N)./max(psdx(1:N)),'-b');                  % plot spectrum computed with fft function
txt{1} = 'FFT spectrum';
leg(2) = plot(F,P./max(P),'-r');                                            % plot spectrum computed with Pwelch function
txt{2} = 'PWelch spectrum';
leg(3) = plot(Sest.f,Sest.S./max(Sest.S),'-g');                             % plot spectrum computed with wafo toolbox
txt{3} = 'Wafo spectrum';
ct = 4;                                                                     % initialise counter
if lpf ~= 0                                                                 % low pass filter
    leg(ct) = plot(fax_Hz(1:N),abs(R1),'-k');                               % plot frequency response 
    txt{ct} = 'Low-pass filter';
    ct = ct + 1;
end
if hpf ~= 0                                                                 % high pass filter
    leg(ct) = plot(fax_Hz(1:N),abs(R2),'.-k');                              % plot frequency response 
    txt{ct} = 'High-pass filter';
    ct = ct + 1;
end
if ubpf ~= 0
     for cnt = 1:length(ubpf)
          leg(ct) = plot(fax_Hz(1:N),abs(R3(:,cnt)),':k');                  % plot frequency response 
          txt{ct} = 'Bandwidth filter';
     end
end
set(gca,'FontSize',24);
xlabel('Frequency [Hz]');ylabel('Normalized magnitude [-]');
title('Frequency domain');
legend(leg,txt,'FontSize',24);
axis([0 lpf+30 0 1.2]);
hold off;

% plot of time series of signal after each filter step
f3 = figure('units','normalized','outerposition',[0 0 1 1]); 
hold on; grid minor; 
plot(time,in(:,i),'-','Color',[0.5 0.5 0.5]);                               % original signal                    
plot(time,detrY,'-','Color',[1 0 0]);                                       % after applying Matlab detrend function
plot(time,noDrift,'-','Color',[0 1 0]);                                     % after removing polynomial best-fit ine
plot(time,f_y,'-','Color',[0.5 0 0]);                                       % best-fit polynom
plot(time,FIL,'-','Color',[0 0 1]);                                         % after FFT filtering
set(gca,'FontSize',24);
xlabel('Time [s]'); ylabel('Force [kN]'); 
title('Time series after each post-processing step'); 
legend('1. Measured signal', '2. After removing best-fit straight line',...
    '3. After removing polynomial best-fit line', '4. Polynomia best-fit line',...
    '5. After FFT-filtering','FontSize',24); 
hold off;
 
%% SAVE OUTPUT

fileID = fopen(['.\Output\Force\Loadcell_' num2str(i)...
    '_filtered.txt'],'w');                                                  % save filtered time series
fprintf(fileID,'%s\t%s\n','Time [s]',['Loadcell_' num2str(i) ' [kN]']);
fprintf(fileID,'%.4f\t%.4f\n',cat(1,time',FIL'));
fclose(fileID);

fileID = fopen(['.\Output\Force\Loadcell_' num2str(i)...
    '_Metadata.txt'],'w');                                                  % save meta data
fprintf(fileID,'%s\t%s\n','original file:',file);
fprintf(fileID,'%s\t%s\n','sample frequency [Hz]:',num2str(fs));
fprintf(fileID,'%s\t%s\n','high-pass filter [Hz]:',num2str(hpf));
fprintf(fileID,'%s\t%s\n','low-pass filter [Hz]:',num2str(lpf));
fprintf(fileID,'%s\t%s\n','bandwidth filter [Hz]:',...
    [num2str(lbpf) '-' num2str(ubpf)]);
fprintf(fileID,'%s\t%s\n','filter design:','butterworth');
fprintf(fileID,'%s\t%s\n','force unit:','kN');
fclose(fileID);

savefig(f1,['.\Output\Force\Loadcell_' num2str(i) '_measured.fig']);        % save figures
savefig(f2,['.\Output\Force\Loadcell_' num2str(i) '_spectrum.fig']);               
savefig(f3,['.\Output\Force\Loadcell_' num2str(i)...
    '_postprocessing steps.fig']);

close all;

end





