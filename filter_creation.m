clear; close all; clc

% inizialization of variables
fs = 500;
sub_1 = load("Data\Subject01_1.mat");

% band definition
banda_delta = [0.5, 4];
banda_theta = [ 4,  8];
banda_alpha = [ 8, 13];
banda_beta  = [ 13, 30];

% creation of the coefficients of filter
[~,filter_delta] = bandpass(sub_1(1).C3, banda_delta, fs);
[~,filter_theta] = bandpass(sub_1(1).C3, banda_theta, fs);
[~,filter_alpha] = bandpass(sub_1(1).C3, banda_alpha, fs);
[~,filter_beta] = bandpass(sub_1(1).C3, banda_beta, fs);

%save('filters.mat', 'filter_beta', 'filter_alpha', 'filter_theta', 'filter_delta');