clear all;
close all;
clc;
listOfFiles = dir("Data\Subject0*_1.mat")
Fs = 500;
signal1_1 = load("Data\Subject01_1.mat");
signal1_2 = load("Data\Subject01_2.mat");
load("filter_fir_teta_1.mat");
%subjectsXchannels = cell(length(listOfFiles), length(signal1_1));
% 
% for
% 
% end

%% band pass teta1
EEG = signal1_2.P3;
EEG = filtfilt(filter_fir_teta_1,1,EEG); % filtering

N = length(EEG); 
t = linspace(0,N,N); % time definition

order = 20;
[PSD, f] = pyulear(EEG, order, (N-1)*2, Fs, "onesided");
fprintf('calcolo: %d \n', sum(PSD(3*500:7*500)));
figure
plot(f,PSD)
hold on
xlabel('Frequency [Hz]')
ylabel('PSD [V^2*Hz^{-1}]')


%% baseline

EEG = signal1_1.P3;
EEG = EEG(31000:31000*2);
EEG = filtfilt(filter_fir_teta_1,1,EEG); % filtering

N = length(EEG); 
t = linspace(0,N,N); % time definition

order = 20;
[PSD, f] = pyulear(EEG, order, (N-1)*2, Fs, "onesided");
fprintf('riposo: %d \n', sum(PSD(3*500:7*500)));
plot(f,PSD)
% xlabel('Frequency [Hz]')
% ylabel('PSD [V^2*Hz^{-1}]')
xlim([3, 7])
legend("Calcolo", "Riposo")



% [magnitude,freq]=periodogram(EEG,rectwin(N),N,Fs);
% 
% figure
% plot(t,EEG)
% xlabel('Time [s]')
% ylabel('\mu V')
% 
% figure
% semilogy(freq,magnitude)
% xlabel('Frequency [Hz]')
% ylabel('\mu log V^{2} /Hz')
