%% Initialization and loading EEGs

clear; close; clc;
fig=0; %set 1 for true
datafolder='Data\';
flist_1 = dir([datafolder,'Subject0*_1.mat']); % baseline has 91000 samples 
flist_2 = dir([datafolder,'Subject0*_2.mat']); % task has 31000 samples

N_sub = numel(flist_1); % Number of subjects

% NOTE: not always the fastest code is the prittiest :), but :_)

% external loading of the fist subject to than copy the structure
subj_rest_1 = importdata([datafolder,flist_1(1).name]);
subj_task_1 = load([datafolder,flist_2(1).name]);

% definition arries of subjects during rest(baseline) and task(mental arithmetic)
subj_rest = repmat(subj_rest_1, 1, N_sub);
subj_task = repmat(subj_task_1, 1, N_sub);
clear subj_rest_1;
clear subj_task_1;

% loading subjects EEGs
for i = 2:N_sub
    subj_rest(i) = load([datafolder,flist_1(i).name]);
    subj_task(i) = load([datafolder,flist_2(i).name]);
end

channels = fieldnames(subj_task); % available channels 
N_chan = numel(channels); % Number of available channels


%% View signal and fft
% pre-processing already performed:
% High-pass filter 0.5 Hz --> insead of detrend
% Notch filter 50 Hz --> remove power network noise

fs = 500;
fn = fs/2;

% Number of samples
N_rest = length(subj_rest(1).C3); % 91000
N_task = length(subj_task(1).C3); % 31000

t_rest = linspace(0,fs,N_rest);
t_task = linspace(0,fs,N_task);

% for s = 1 : N_sub
%     figure(s)
%     for c = 1 : N_chan
%         subplot(5,4,c)
%         plot(t_1, subj_rest(s).(channels{c}))
%         title('EEG rest')
%     end
%         
% end
% 
% for s = 1 : N_sub
%     figure(s+N_sub)
%     for c = 1 : N_chan
%         subplot(5,4,c)
%         plot(t_2, subj_task(s).(channels{c}))
%         title('EEG mental arith')
%     end
%         
% end
% 

%TODO: cambiare da questa visualizzazione a quella con la mappa del
%cervello e al click dell'elettrodo esce la fft associata, comunque 6
%finestre per ogni soggetto, o no? da vedere
if fig==1
    for s = 1:N_sub
        figure(s+N_sub*2)
        for c = 1 : N_chan
            EEG_rest_fft = fft(subj_rest(s).(channels{c}));
            subplot(5,4,c)
            plot(t_rest, abs(EEG_rest_fft))
            % plot(t_1,(1/(fs*N_1)) * abs(EEG_rest_fft).^2) PSD
            title({'FFT EEG rest ', channels{c}})
            xlim([0,fn])
        end
    end
end
%% PSD all signal

noverlap = 0;
l_wind = fs * 20; % window of 20 s

for ch = 1 : N_chan
    [PSDp_rest_1.(channels{ch}),fp_1] = pwelch(subj_rest(1).(channels{ch}), rectwin(l_wind), noverlap, [], fs);
    [PSDp_task_1.(channels{ch}),fp_2] = pwelch(subj_task(1).(channels{ch}), rectwin(l_wind), noverlap, [], fs);
end

PSDp_rest = repmat(PSDp_rest_1, 1, N_sub);
PSDp_task = repmat(PSDp_task_1, 1, N_sub);
clear PSDp_rest_1;
clear PSDp_task_1;

for ch = 1 : N_chan
    for s = 2 : N_sub

        [PSDp_rest(s).(channels{ch}),fp_1] = pwelch(subj_rest(s).(channels{ch}), rectwin(l_wind), noverlap, [], fs);
        [PSDp_task(s).(channels{ch}),fp_2] = pwelch(subj_task(s).(channels{ch}), rectwin(l_wind), noverlap, [], fs);
        
        % figure(ch)
        % subplot(6,2,s*2-1)
        % plot(fp_1,PSDp_rest(s).(channels{ch}))
        % xlim([0,fn])
        % title('PSD rest')
        % 
        % subplot(6,2,s*2)
        % plot(fp_2,PSDp_task(s).(channels{ch}))
        % xlim([0,fn])
        % title('PSD mental arith')
    end
end

%% Resample the signal

fcut = 70;
fs_new = 140;
EEG_rest_rs = struct();
EEG_task_rs = struct();

for s = 1 : N_sub
    for c = 1 : N_chan
        EEG_rest_rs(s).(channels{c}) = resample(subj_rest(s).(channels{c}),fs_new,fs); %basically, take one sample every fs/fs_new
        t_1_rs = resample(t_rest,fs_new,fs);
        N_1_rs = length(EEG_rest_rs);
        EEG_task_rs(s).(channels{c}) = resample(subj_task(s).(channels{c}),fs_new,fs);
        t_2_rs = resample(t_task,fs_new,fs);
        N_2_rs = length(EEG_task_rs);

%         figure(s)
%         subplot(5,4,c)
%         plot(t_1,subj_rest(s).(channels{c}))
%         hold on
%         plot(t_1_rs,EEG_rest_rs(s).(channels{c}))
%         %legend('original EEG', 'resampled EEG')
    end
end
tmp = EEG_rest_rs(1).C3;
save("Data\baseline_resampled.mat", 'tmp');

%% PSD all resampled signal
lim = [0, fcut];
l_wind_rs = fs_new * 10;

for c = 1 : N_chan
    [PSDp_rest_rs_1.(channels{c}),~] = pwelch(EEG_rest_rs(s).(channels{c}),rectwin(l_wind_rs),noverlap,[],fs_new);
    [PSDp_task_rs_1.(channels{c}),~] = pwelch(EEG_task_rs(s).(channels{c}),rectwin(l_wind_rs),noverlap,[],fs_new);
end

PSDp_rest_rs = repmat(PSDp_rest_rs_1, 1, N_sub);
PSDp_task_rs = repmat(PSDp_task_rs_1, 1, N_sub);
clear PSDp_rest_rs_1;
clear PSDp_task_rs_1;

for c = 1 : N_chan
    for s = 1 : N_sub

        [PSDp_rest_rs(s).(channels{c}),~] = pwelch(EEG_rest_rs(s).(channels{c}),rectwin(l_wind_rs),noverlap,[],fs_new);
        [PSDp_task_rs(s).(channels{c}),~] = pwelch(EEG_task_rs(s).(channels{c}),rectwin(l_wind_rs),noverlap,[],fs_new);

%         figure(c)
%         subplot(6,2,rest(s))
%         plot(fp_1_rs,PSDp_rest_rs(s).(channels{c}))
%         xlim(lim)
%         title('PSD rest')
% 
%         subplot(6,2,aritm(s))
%         plot(fp_2_rs,PSDp_task_rs(s).(channels{c}))
%         xlim(lim)
%         title('PSD mental arith')
    end
end

%% Band_subdivision

load('filters.mat');

band_coefficients= struct('delta',  filter_delta.Coefficients, ...
    'theta',  filter_theta.Coefficients,...
    'alpha',  filter_alpha.Coefficients,...
    'beta',  filter_beta.Coefficients);

[resample_band_1, resample_band_2] = Band_sub(band_coefficients, fs_new, fs);

%% PDS band
bands = ["delta", "theta", "alpha", "beta"];
lim_delta = [0.5, 4];
lim_theta = [4,  8];
lim_alpha = [8, 13];
lim_beta  = [13, 30];
overlap = floor(l_wind_rs*0.5);

%TODO: si puo fare in modo diverso?
[tmp ,fp] = pwelch(resample_band_1.(channels{1}){s,1},rectwin(l_wind_rs),overlap,[],fs_new);
len_welch = numel(tmp);

for c = 1 : N_chan
    PSD_rest.(channels{c}) = zeros(N_sub, numel(bands), len_welch);
    PSD_task.(channels{c}) = zeros(N_sub, numel(bands), len_welch);
end

for c = 1 : N_chan
    for s = 1 : N_sub
        for b = 1 : numel(bands)
            [PSD_rest.(channels{c})(s,b,:),~] = pwelch(resample_band_1.(channels{c}){s,b},rectwin(l_wind_rs),overlap,[],fs_new);
            [PSD_task.(channels{c})(s,b,:),~] = pwelch(resample_band_2.(channels{c}){s,b},rectwin(l_wind_rs),overlap,[],fs_new);
        end
    end
end


%% Band and channels selected
load('filters_sub_bands.mat')
band_coefficients= struct('theta1',  filter_theta1.Coefficients, ...
    'theta2',  filter_theta2.Coefficients,...
    'beta1',  filter_beta1.Coefficients,...
    'beta2',  filter_beta2.Coefficients);

[signal_sub_band_1, signal_sub_band_2] = Band_sub(band_coefficients, fs_new, fs);

%% PDS band
lim_theta1 = [4.1, 5.8];
lim_theta2 = [5.9,  7.4];
lim_beta1 = [13, 19.9];
lim_beta2 = [20, 25];

for c = 1 : N_chan
    PSD_rest_theta.(channels{c}) = zeros(N_sub, 2, len_welch);
    PSD_rest_beta.(channels{c}) = zeros(N_sub, 2, len_welch);

    PSD_task_theta.(channels{c}) = zeros(N_sub, 2, len_welch);
    PSD_task_beta.(channels{c}) = zeros(N_sub, 2, len_welch);
end

for c = 1 : N_chan
    for s = 1 : N_sub
        for i = 1:2 %sub bands
            %Theta
            [PSD_rest_theta.(channels{c})(s,i,:), ~] = pwelch(signal_sub_band_1.(channels{c}){s,i},rectwin(l_wind_rs),overlap,[],fs_new);
            [PSD_task_theta.(channels{c})(s,i,:), ~] = pwelch(signal_sub_band_2.(channels{c}){s,i},rectwin(l_wind_rs),overlap,[],fs_new);
            %Beta
            [PSD_rest_beta.(channels{c})(s,i,:), ~] = pwelch(signal_sub_band_1.(channels{c}){s,i+2},rectwin(l_wind_rs),overlap,[],fs_new);
            [PSD_task_beta.(channels{c})(s,i,:), ~] = pwelch(signal_sub_band_2.(channels{c}){s,i+2},rectwin(l_wind_rs),overlap,[],fs_new);
        end
    end
end


%% Area
N_chan = 4;
chan_choosen = [4, 7, 8, 12]; % F3, F8, Fp1, O2
area_rest = zeros(6, 4, 4); %TODO
area_task = zeros(6, 4, 4); %TODO
for c = 1 : N_chan
    area_rest(:,c,1:2) = sum(PSD_rest_theta.(channels{chan_choosen(c)}), 3); %TODO
    area_rest(:,c,3:4) = sum(PSD_rest_beta.(channels{chan_choosen(c)}), 3); %TODO

    area_task(:,c,1:2) = sum(PSD_task_theta.(channels{chan_choosen(c)}), 3); %TODO
    area_task(:,c,3:4) = sum(PSD_task_beta.(channels{chan_choosen(c)}), 3); %TODO
end

%% Plot results

band_N = numel(fieldnames(band_coefficients)); %TODO
inc = 0.2;
mul = 0.5;
%inc_b = 1/band_N;
y_lim = 350;
X_ticks = ones(1,band_N);
X_ticks(1) = X_ticks(1) + inc/2;
for i = 2:band_N
    X_ticks(i) = X_ticks(i) + (i-1)*mul + inc/2;
end
if fig==1
    for c = 1:N_chan
        mul = 0.5;
        figure('Name', "Scatter channel "+channels{chan_choosen(c)})
        hold on
        o_flag = false; % outliers flag
        for i = 1:band_N
            % for j = 1:N_sub
            %     rgb(j,1) = rgb(j,1) + 0.1;
            %     rgb(j,2) = rgb(j,2) + 0.23;
            % end
            mul = mul + 0.5;
            scatter(mul*ones(1,6), area_rest(:,c,i),"blue", 'filled'); %TODO 6=N_SUB
            outliers_rest = find(area_rest(:,c,i) > y_lim);
    
            scatter((mul+inc)*ones(1,6), area_task(:,c,i),"red", 'filled');
            outliers_task = find(area_task(:,c,i) > y_lim);
    
            for k = 1:numel(outliers_rest)
                o_flag = true;
                o = scatter(mul, y_lim-k*2, "black");
                o.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Actual Y', area_rest(outliers_rest(k),c,i));
                o.DataTipTemplate.DataTipRows(4) = '';
            end
            for k = 1:numel(outliers_task)
                o_flag = true;
                o = scatter(mul+inc, y_lim-k*2, "black");
                o.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Actual Y', area_task(outliers_task(k),c,i));
                o.DataTipTemplate.DataTipRows(4) = '';
            end
        end
        xticks(X_ticks);
        xticklabels({'\theta1', ...
            '\theta2','\beta1', ...
            '\beta2'});
        xlim([X_ticks(1)-inc, X_ticks(numel(X_ticks))+inc]);
        ylim([0, y_lim]);
        if o_flag
            legend("rest", "task", "outliers", "Location", "eastoutside");
        else
            legend("rest", "task", "Location", "eastoutside");
        end
    end
end
%% Median
median_rest = median(area_rest,1);
median_task = median(area_task,1);


arraytable = zeros(size(median_rest, 3),size(median_rest, 3)*2);
for b=1:size(median_rest, 3)
    arraytable(:,b*2-1) = median_rest(:,:,b);
    arraytable(:,b*2) = median_task(:,:,b);
end


variableNames = {'Theta1_rest', 'Theta1_task', 'Theta2_rest', 'Theta2_task',...
                'Beta1_rest', 'Beta1_task', 'Beta2_rest', 'Beta2_task'};
tableData = array2table(arraytable, 'VariableNames', variableNames);

writetable(tableData, 'Tabella_Mediane.xlsx')

%% 
area_rest = zeros(N_sub,numel(channels), numel(bands));
area_task = zeros(N_sub,numel(channels), numel(bands));

for c = 1 : numel(channels)
    area_rest(:,c,1:2) = sum(PSD_rest_theta.(channels{c}), 3); %TODO
    area_rest(:,c,3:4) = sum(PSD_rest_beta.(channels{c}), 3); %TODO

    area_task(:,c,1:2) = sum(PSD_task_theta.(channels{c}), 3); %TODO
    area_task(:,c,3:4) = sum(PSD_task_beta.(channels{c}), 3); %TODO
end

median_rest = median(area_rest,1);
median_task = median(area_task,1);
%%
load("Data\chanlocs.mat")
names = {chanlocs(1:19).labels};
[~, idx] = sort(names);
chanlocs = chanlocs(idx);

%%
cd("helper_code")
figure('Name', "Theta1")
subplot(1,2,1)
title("Rest")
topoplot(median_rest(:,:,1), chanlocs, 'electrodes', 'labels', 'style', 'fill')
subplot(1,2,2)
title("Task")
topoplot(median_task(:,:,1), chanlocs, 'electrodes', 'labels', 'style', 'fill')

figure('Name', "Theta2")
subplot(1,2,1)
title("Rest")
topoplot(median_rest(:,:,2), chanlocs, 'electrodes', 'labels', 'style', 'fill')
subplot(1,2,2)
title("Task")
topoplot(median_task(:,:,2), chanlocs, 'electrodes', 'labels', 'style', 'fill')

figure('Name', "Beta1")
subplot(1,2,1)
title("Rest")
topoplot(median_rest(:,:,3), chanlocs, 'electrodes', 'labels', 'style', 'fill')
subplot(1,2,2)
title("Task")
topoplot(median_task(:,:,3), chanlocs, 'electrodes', 'labels', 'style', 'fill')

figure('Name', "Beta2")
subplot(1,2,1)
title("Rest")
topoplot(median_rest(:,:,4), chanlocs, 'electrodes', 'labels', 'style', 'fill')
subplot(1,2,2)
title("Task")
topoplot(median_task(:,:,4), chanlocs, 'electrodes', 'labels', 'style', 'fill')

cd('..')
%%
Ratio = abs(((median_task./median_rest)-1).*100);
bnd=4;
tmp=Ratio(:,:,bnd);
Ratio(:,isoutlier(tmp),bnd)
find(isoutlier(tmp))
%%
bnd=4;
tmp=Ratio(:,:,bnd);
Q3 = quantile(tmp, 0.75);
Q1 = quantile(tmp, 0.25);
position_high = ((tmp >= Q3));
position_low = ((tmp <= Q1));
% Ratio(:,position,bnd)
find(position_high)
find(position_low)
%%
