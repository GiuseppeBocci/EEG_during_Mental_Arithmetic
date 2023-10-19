clear; close; clc

datafolder='Data\';
flist_1=dir([datafolder,'Subject0*_1.mat']); %91000
flist_2=dir([datafolder,'Subject0*_2.mat']); %31000

for i=1:numel(flist_1)
    fname_1=flist_1(i).name;
    fpath_1=flist_1(i).folder;
    sub_1(i)=load([fpath_1,'/',fname_1]);

    fname_2=flist_2(i).name;
    fpath_2=flist_2(i).folder;
    sub_2(i)=load([fpath_2,'/',fname_2]);
end

channels = fieldnames(sub_2);
N_chan = numel(channels);
N_sub = numel(sub_2);

%% View signal and fft
%Pre-processing
%High-pass filter 0.5 Hz --> al posto di detrend
%Low-pass filter 45 Hz
%Notch filter 50 Hz

fs = 500;
res = 0.3; %time constant of the amplification tract in seconds
N_1 = length(sub_1(1).C3);
N_2 = length(sub_2(1).C3);

%N_min = fs*res; %fs/N aumenta N, res in Hz diminuisce, res in s aumenta

t_1 = linspace(0,fs,N_1);
t_2 = linspace(0,fs,N_2);

% for s = 1 : N_sub
%     figure(s)
%     for c = 1 : N_chan
%         subplot(5,4,c)
%         plot(t_1, sub_1(s).(channels{c}))
%         title('EEG rest')
%     end
%         
% end
% 
% for s = 1 : N_sub
%     figure(s+N_sub)
%     for c = 1 : N_chan
%         subplot(5,4,c)
%         plot(t_2, sub_2(s).(channels{c}))
%         title('EEG mental arith')
%     end
%         
% end
% 
% for s = 1 : N_sub
%     figure(s+N_sub*2)
%     for c = 1 : N_chan
%         EEG_1_fft=fft(sub_1(s).(channels{c}));
%         subplot(5,4,c)
%         plot(t_1, abs(EEG_1_fft))
%         %plot(t_1,(1/(fs*N_1)) * abs(EEG_1_fft).^2)
%         title('FFT EEG')
%         xlim([0,fs/2])
%     end
% end

%% PSD all signal
lim = [0,fs/2];
noverlap = 0;
l_wind = fs * 20; %window of 20 s
rest = [1, 3, 5, 7 , 9, 11];
aritm = [2, 4, 6, 8, 10, 12];

for c = 1 : N_chan
    for s = 1 : N_sub

        [PSDp_1(s).(channels{c}),fp_1] = pwelch(sub_1(s).(channels{c}),rectwin(l_wind),noverlap,[],fs);
        [PSDp_2(s).(channels{c}),fp_2] = pwelch(sub_2(s).(channels{c}),rectwin(l_wind),noverlap,[],fs);

%         figure(c)
%         subplot(6,2,rest(s))
%         plot(fp_1,PSDp_1(s).(channels{c}))
%         xlim(lim)
%         title('PSD rest')
% 
%         subplot(6,2,aritm(s))
%         plot(fp_2,PSDp_2(s).(channels{c}))
%         xlim(lim)
%         title('PSD mental arith')
    end
end

%% Resample the signal
fcut = 70;
fs_new = 140;
EEG_1_rs = struct();
EEG_2_rs = struct();

for s = 1 : N_sub
    for c = 1 : N_chan
        EEG_1_rs(s).(channels{c}) = resample(sub_1(s).(channels{c}),fs_new,fs); %basically, take one sample every fs/fs_new
        t_1_rs = resample(t_1,fs_new,fs);
        N_1_rs = length(EEG_1_rs);
        EEG_2_rs(s).(channels{c}) = resample(sub_2(s).(channels{c}),fs_new,fs);
        t_2_rs = resample(t_2,fs_new,fs);
        N_2_rs = length(EEG_2_rs);

%         figure(s)
%         subplot(5,4,c)
%         plot(t_1,sub_1(s).(channels{c}))
%         hold on
%         plot(t_1_rs,EEG_1_rs(s).(channels{c}))
%         %legend('original EEG', 'resampled EEG')
    end
end

%% PSD all resampled signal
lim = [0, fcut];
l_wind_rs = fs_new * 10;

for c = 1 : N_chan
    for s = 1 : N_sub

        [PSDp_1_rs(s).(channels{c}),fp_1_rs] = pwelch(EEG_1_rs(s).(channels{c}),rectwin(l_wind_rs),noverlap,[],fs_new);
        [PSDp_2_rs(s).(channels{c}),fp_2_rs] = pwelch(EEG_2_rs(s).(channels{c}),rectwin(l_wind_rs),noverlap,[],fs_new);

%         figure(c)
%         subplot(6,2,rest(s))
%         plot(fp_1_rs,PSDp_1_rs(s).(channels{c}))
%         xlim(lim)
%         title('PSD rest')
% 
%         subplot(6,2,aritm(s))
%         plot(fp_2_rs,PSDp_2_rs(s).(channels{c}))
%         xlim(lim)
%         title('PSD mental arith')
    end
end

%% Band_subdivision
load('filters_resample.mat')
coefficienti_Bande= struct('delta',  filter_delta.Coefficients, ...
    'theta',  filter_theta.Coefficients,...
    'alpha',  filter_alpha.Coefficients,...
    'beta',  filter_beta.Coefficients);

[resample_band_1, resample_band_2] = Band_sub(coefficienti_Bande, fs_new, fs);

%% PDS band
lim_delta = [0.5, 4];
lim_theta = [4,  8];
lim_alpha = [8, 13];
lim_beta  = [13, 30];
overlap = 0.5;
folder = 'Band\';

for c = 1 : N_chan
    for s = 1 : N_sub
        %Delta
        [PSD_1_delta(s).(channels{c}),fp_1_delta] = pwelch(resample_band_1.(channels{c}){s,1},rectwin(l_wind_rs),overlap,[],fs_new);
        [PSD_2_delta(s).(channels{c}),fp_2_delta] = pwelch(resample_band_2.(channels{c}){s,1},rectwin(l_wind_rs),overlap,[],fs_new);
        %Theta
        [PSD_1_theta(s).(channels{c}),fp_1_theta] = pwelch(resample_band_1.(channels{c}){s,2},rectwin(l_wind_rs),overlap,[],fs_new);
        [PSD_2_theta(s).(channels{c}),fp_2_theta] = pwelch(resample_band_2.(channels{c}){s,2},rectwin(l_wind_rs),overlap,[],fs_new);
        %Alpha
        [PSD_1_alpha(s).(channels{c}),fp_1_alpha] = pwelch(resample_band_1.(channels{c}){s,3},rectwin(l_wind_rs),overlap,[],fs_new);
        [PSD_2_alpha(s).(channels{c}),fp_2_alpha] = pwelch(resample_band_2.(channels{c}){s,3},rectwin(l_wind_rs),overlap,[],fs_new);
        %Beta
        [PSD_1_beta(s).(channels{c}),fp_1_beta] = pwelch(resample_band_1.(channels{c}){s,4},rectwin(l_wind_rs),overlap,[],fs_new);
        [PSD_2_beta(s).(channels{c}),fp_2_beta] = pwelch(resample_band_2.(channels{c}){s,4},rectwin(l_wind_rs),overlap,[],fs_new);

%         %Delta
%         figure(c)
%         subplot(6,2,rest(s))
%         plot(fp_1_delta,PSD_1_delta(s).(channels{c}))
%         xlim(lim_delta)
%         title('PSD rest')
%         subplot(6,2,aritm(s))
%         plot(fp_2_delta,PSD_2_delta(s).(channels{c}))
%         xlim(lim_delta)
%         title('PSD mental arith')
% 
%         %Theta
%         figure(c+N_sub)
%         subplot(6,2,rest(s))
%         plot(fp_1_theta,PSD_1_theta(s).(channels{c}))
%         xlim(lim_theta)
%         title('PSD rest')
%         subplot(6,2,aritm(s))
%         plot(fp_2_theta,PSD_2_theta(s).(channels{c}))
%         xlim(lim_theta)
%         title('PSD mental arith')
% 
%         %Alpha
%         figure(c+N_sub*2)
%         subplot(6,2,rest(s))
%         plot(fp_1_alpha,PSD_1_alpha(s).(channels{c}))
%         xlim(lim_alpha)
%         title('PSD rest')
%         subplot(6,2,aritm(s))
%         plot(fp_2_alpha,PSD_2_alpha(s).(channels{c}))
%         xlim(lim_alpha)
%         title('PSD mental arith')
% 
%         %Beta
%         figure(c+N_sub*3)
%         subplot(6,2,rest(s))
%         plot(fp_1_beta,PSD_1_beta(s).(channels{c}))
%         xlim(lim_beta)
%         title('PSD rest')
%         subplot(6,2,aritm(s))
%         plot(fp_2_beta,PSD_2_beta(s).(channels{c}))
%         xlim(lim_beta)
%         title('PSD mental arith')
    end
%          %Salvare figura
%          hfig=figure((c));
%          titolo_figura=sprintf('%s_delta', channels{c});
%          fullFileName=fullfile(folder,titolo_figura);
%          saveas(hfig,fullFileName)
%  
%          hfig=figure((c+N_sub));
%          titolo_figura=sprintf('%s_theta', channels{c});
%          fullFileName=fullfile(folder,titolo_figura);
%          saveas(hfig,fullFileName)
%  
%          hfig=figure((c+N_sub*2));
%          titolo_figura=sprintf('%s_alpha', channels{c});
%          fullFileName=fullfile(folder,titolo_figura);
%          saveas(hfig,fullFileName)
%  
%          hfig=figure((c+N_sub*3));
%          titolo_figura=sprintf('%s_beta', channels{c});
%          fullFileName=fullfile(folder,titolo_figura);
%          saveas(hfig,fullFileName)
end


%% Band and channels selected
load('filters_sub_bands.mat')
coefficienti_Bande= struct('theta1',  filter_theta1.Coefficients, ...
    'theta2',  filter_theta2.Coefficients,...
    'beta1',  filter_beta1.Coefficients,...
    'beta2',  filter_beta2.Coefficients);

[signal_sub_band_1, signal_sub_band_2] = Band_sub(coefficienti_Bande, fs_new, fs);

%% PDS band
lim_theta1 = [4.1, 5.8];
lim_theta2 = [5.9,  7.4];
lim_beta1 = [13, 19.9];
lim_beta2 = [20, 25];
folder = 'Band_selected\';
chan_choosen = [4, 7, 8, 12]; % F3, F8, Fp1, O2

PSD_1_theta1 = struct();

for c = 1 : numel(chan_choosen)
    for s = 1 : N_sub
        %Theta1
        [PSD_1_theta1(s).(channels{chan_choosen(c)}),fp_1_theta1] = pwelch(signal_sub_band_1.(channels{chan_choosen(c)}){s,1},rectwin(l_wind_rs),overlap,[],fs_new);
        [PSD_2_theta1(s).(channels{chan_choosen(c)}),fp_2_theta1] = pwelch(signal_sub_band_2.(channels{chan_choosen(c)}){s,1},rectwin(l_wind_rs),overlap,[],fs_new);
        %Theta2
        [PSD_1_theta2(s).(channels{chan_choosen(c)}),fp_1_theta2] = pwelch(signal_sub_band_1.(channels{chan_choosen(c)}){s,2},rectwin(l_wind_rs),overlap,[],fs_new);
        [PSD_2_theta2(s).(channels{chan_choosen(c)}),fp_2_theta2] = pwelch(signal_sub_band_2.(channels{chan_choosen(c)}){s,2},rectwin(l_wind_rs),overlap,[],fs_new);
        %Beta1
        [PSD_1_beta1(s).(channels{chan_choosen(c)}),fp_1_beta1] = pwelch(signal_sub_band_1.(channels{chan_choosen(c)}){s,3},rectwin(l_wind_rs),overlap,[],fs_new);
        [PSD_2_beta1(s).(channels{chan_choosen(c)}),fp_2_beta1] = pwelch(signal_sub_band_2.(channels{chan_choosen(c)}){s,3},rectwin(l_wind_rs),overlap,[],fs_new);
        %Beta2
        [PSD_1_beta2(s).(channels{chan_choosen(c)}),fp_1_beta2] = pwelch(signal_sub_band_1.(channels{chan_choosen(c)}){s,4},rectwin(l_wind_rs),overlap,[],fs_new);
        [PSD_2_beta2(s).(channels{chan_choosen(c)}),fp_2_beta2] = pwelch(signal_sub_band_2.(channels{chan_choosen(c)}){s,4},rectwin(l_wind_rs),overlap,[],fs_new);

%         %Theta1
%         figure(c)
%         subplot(6,2,rest(s))
%         plot(fp_1_theta1,PSD_1_theta1(s).(channels{chan_choosen}))
%         xlim(lim_delta)
%         title('PSD rest')
%         subplot(6,2,aritm(s))
%         plot(fp_2_theta1,PSD_2_theta1(s).(channels{chan_choosen}))
%         xlim(lim_delta)
%         title('PSD mental arith')
% 
%         %Theta2
%         figure(c+N_sub)
%         subplot(6,2,rest(s))
%         plot(fp_1_theta2,PSD_1_theta2(s).(channels{chan_choosen}))
%         xlim(lim_theta)
%         title('PSD rest')
%         subplot(6,2,aritm(s))
%         plot(fp_2_theta2,PSD_2_theta2(s).(channels{chan_choosen}))
%         xlim(lim_theta)
%         title('PSD mental arith')
% 
%         %Beta1
%         figure(c+N_sub*2)
%         subplot(6,2,rest(s))
%         plot(fp_1_beta1,PSD_1_beta1(s).(channels{chan_choosen(c)}))
%         xlim(lim_alpha)
%         title('PSD rest')
%         subplot(6,2,aritm(s))
%         plot(fp_2_beta1,PSD_2_beta1(s).(channels{chan_choosen(c)}))
%         xlim(lim_alpha)
%         title('PSD mental arith')
% 
%         %Beta2
%         figure(c+N_sub*3)
%         subplot(6,2,rest(s))
%         plot(fp_1_beta2,PSD_1_beta2(s).(channels{chan_choosen(c)}))
%         xlim(lim_beta)
%         title('PSD rest')
%         subplot(6,2,aritm(s))
%         plot(fp_2_beta2,PSD_2_beta2(s).(channels{chan_choosen(c)}))
%         xlim(lim_beta)
%         title('PSD mental arith')
    end
%          %Salvare figura
%          hfig=figure((c));
%          titolo_figura=sprintf('%s_delta', channels{chan_choosen});
%          fullFileName=fullfile(folder,titolo_figura);
%          saveas(hfig,fullFileName)
%  
%          hfig=figure((c+N_sub));
%          titolo_figura=sprintf('%s_theta', channels{chan_choosen});
%          fullFileName=fullfile(folder,titolo_figura);
%          saveas(hfig,fullFileName)
%  
%          hfig=figure((c+N_sub*2));
%          titolo_figura=sprintf('%s_alpha', channels{chan_choosen});
%          fullFileName=fullfile(folder,titolo_figura);
%          saveas(hfig,fullFileName)
%  
%          hfig=figure((c+N_sub*3));
%          titolo_figura=sprintf('%s_beta', channels{chan_choosen});
%          fullFileName=fullfile(folder,titolo_figura);
%          saveas(hfig,fullFileName)
end


%% Area
N_chan = 4;

for s = 1 : N_sub
    for c = 1 : N_chan
        area_1{s,c,1} = sum(PSD_1_theta1(s).(channels{chan_choosen(c)}));
        area_1{s,c,2} = sum(PSD_1_theta2(s).(channels{chan_choosen(c)}));
        area_2{s,c,1} = sum(PSD_2_theta1(s).(channels{chan_choosen(c)}));
        area_2{s,c,2} = sum(PSD_2_theta2(s).(channels{chan_choosen(c)}));
        area_1{s,c,3} = sum(PSD_1_beta1(s).(channels{chan_choosen(c)}));
        area_1{s,c,4} = sum(PSD_1_beta2(s).(channels{chan_choosen(c)}));
        area_2{s,c,3} = sum(PSD_2_beta1(s).(channels{chan_choosen(c)}));
        area_2{s,c,4} = sum(PSD_2_beta2(s).(channels{chan_choosen(c)}));
    end
end
disp('a')
%%

s=2;
can=4;
figure
scatter(ones(1,6),cell2mat(area_1(:,1,can)))
hold on
scatter(1.2*ones(1,6), cell2mat(area_2(:,1,can)))
hold on
scatter(1.5*ones(1,6),cell2mat(area_1(:,2,can)))
hold on
scatter(1.7*ones(1,6), cell2mat(area_2(:,2,can)))
hold on
scatter(2*ones(1,6),cell2mat(area_1(:,3,can)))
hold on
scatter(2.2*ones(1,6), cell2mat(area_2(:,3,can)))
hold on
scatter(2.5*ones(1,6),cell2mat(area_1(:,4,can)))
hold on
scatter(2.7*ones(1,6), cell2mat(area_2(:,4,can)))
hold on
plot([1, 1.2], [cell2mat(area_1(s,1,can)), cell2mat(area_2(s,1,can))])
hold on
plot([1.5, 1.7], [cell2mat(area_1(s,2,can)), cell2mat(area_2(s,2,can))])
hold on
plot([2, 2.2], [cell2mat(area_1(s,3,can)), cell2mat(area_2(s,3,can))])
hold on
plot([2.5, 2.7], [cell2mat(area_1(s,4,can)), cell2mat(area_2(s,4,can))])