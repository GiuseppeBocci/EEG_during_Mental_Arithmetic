clear; close all; clc

% inizialization of variables
fs = 500;
sub_1 = load("Data\Subject01_1.mat");

% band definition
banda_delta = [ 1,  4]; % 1 to filter noise
banda_theta = [ 4,  8];
banda_alpha = [ 8, 13];
banda_beta  = [13, 30];

% creation of the coefficients of filter
[~,filter_delta] = bandpass(sub_1(1).C3, banda_delta, fs);
[~,filter_theta] = bandpass(sub_1(1).C3, banda_theta, fs);
[~,filter_alpha] = bandpass(sub_1(1).C3, banda_alpha, fs);
[~,filter_beta] = bandpass(sub_1(1).C3, banda_beta, fs);

%save('filters.mat', 'filter_beta', 'filter_alpha', 'filter_theta', 'filter_delta');

%% Theta1,2 Beta1,2

banda_theta1 = [  4.1, 5.8];
banda_theta2 = [  5.9, 7.4];
banda_beta1  = [ 13,  19.9];
banda_beta2  = [ 20,  25  ];
fs = 500;

tmp = load("Data\baseline_resampled.mat");
EEG_baseline_resampled = struct2cell(tmp);
EEG_baseline_resampled = EEG_baseline_resampled{1};

% creation of the coefficients of filter (all fir)
[~,filter_theta1] = bandpass(EEG_baseline_resampled(1).F8, banda_theta1, fs);
[~,filter_theta2] = bandpass(EEG_baseline_resampled(1).F8, banda_theta2, fs);
[~,filter_beta1]  = bandpass(EEG_baseline_resampled(1).F8, banda_beta1, fs);
[~,filter_beta2]  = bandpass(EEG_baseline_resampled(1).F8, banda_beta2, fs);

save('filters_sub_bands.mat', 'filter_theta1', 'filter_theta2', 'filter_beta1', 'filter_beta2');
