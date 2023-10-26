clear; close all; clc

% inizialization of variables
fs = 140;
tmp = load("Data\baseline_resampled.mat");
EEG_baseline_resampled = struct2cell(tmp);
EEG_baseline_resampled = EEG_baseline_resampled{1};

% band definition
banda_delta = [ 1,  4]; % 1 to filter noise
banda_theta = [ 4,  8];
banda_alpha = [ 8, 13];
banda_beta  = [13, 30];

% creation of the coefficients of filter
[~,filter_delta] = bandpass(EEG_baseline_resampled, banda_delta, fs, 'ImpulseResponse', 'fir');
[~,filter_theta] = bandpass(EEG_baseline_resampled, banda_theta, fs, 'ImpulseResponse', 'fir');
[~,filter_alpha] = bandpass(EEG_baseline_resampled, banda_alpha, fs, 'ImpulseResponse' ,'fir');
[~,filter_beta] = bandpass(EEG_baseline_resampled, banda_beta, fs, 'ImpulseResponse', 'fir');

save('filters.mat', 'filter_beta', 'filter_alpha', 'filter_theta', 'filter_delta');

%% Theta1,2 Beta1,2

banda_theta1 = [  4.1, 5.8];
banda_theta2 = [  5.9, 7.4];
banda_beta1  = [ 13,  19.9];
banda_beta2  = [ 20,  25  ];

% creation of the coefficients of filter (all fir)
[~,filter_theta1] = bandpass(EEG_baseline_resampled, banda_theta1, fs, 'ImpulseResponse','fir');
[~,filter_theta2] = bandpass(EEG_baseline_resampled, banda_theta2, fs, 'ImpulseResponse','fir');
[~,filter_beta1]  = bandpass(EEG_baseline_resampled, banda_beta1, fs, 'ImpulseResponse','fir');
[~,filter_beta2]  = bandpass(EEG_baseline_resampled, banda_beta2, fs, 'ImpulseResponse','fir');

save('filters_sub_bands.mat', 'filter_theta1', 'filter_theta2', 'filter_beta1', 'filter_beta2');
