%% Initialization and loading EEGs

clear; close all; clc;

show_fig = false; % set to true to show all the figures
subj_show = 1; % subject to show when show_fig = true
chan_show = "FP1"; % channel to show when show_fig = true;
datafolder='Data\';
% Palace here your files
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

if show_fig
    c_show = chan_show; % redefine chan_show here for this figure
    figure
    hold on
    plot(t_rest, abs(fft(subj_rest(subj_show).(c_show))));
    plot(t_task, abs(fft(subj_task(subj_show).(c_show))));
    title('FFT EEG '+c_show);
    xlabel("Frequency [Hz]");
    ylabel("Amplitude [\muV^2/Hz]"); %TODO
    xlim([0,fn])

    legend("rest","task",'ItemHitFcn',@show_hide_graph);

    helpdlg("Click the legend to hide/show a graph");
end

%% PSD all signal

noverlap = 0;
l_wind = fs * 20; % window of 20 s
rect_w = rectwin(l_wind);

for ch = 1 : N_chan
    [PSDp_rest_1.(channels{ch}),fp_rest] = pwelch(subj_rest(1).(channels{ch}), rect_w, noverlap, [], fs);
    [PSDp_task_1.(channels{ch}),fp_task] = pwelch(subj_task(1).(channels{ch}), rect_w, noverlap, [], fs);
    
end

PSDp_rest = repmat(PSDp_rest_1, 1, N_sub);
PSDp_task = repmat(PSDp_task_1, 1, N_sub);
clear PSDp_rest_1;
clear PSDp_task_1;

for ch = 1 : N_chan
    for s = 2 : N_sub

        [PSDp_rest(s).(channels{ch}),~] = pwelch(subj_rest(s).(channels{ch}), rect_w, noverlap, [], fs);
        [PSDp_task(s).(channels{ch}),~] = pwelch(subj_task(s).(channels{ch}), rect_w, noverlap, [], fs);
       
    end
end

if show_fig
    c_show = chan_show; % redefine chan_show here for this figure
    
    figure
    hold on
    xi = repmat((1:N_chan)', 1,numel(fp_rest))*3;
    yi= repmat(fp_rest', N_chan, 1);
    
    r = waterfall(xi, yi, (reshape(cell2mat(struct2cell(PSDp_rest(subj_show))), numel(fp_rest), numel(channels))'));
    t = waterfall(xi-0.5, yi, (reshape(cell2mat(struct2cell(PSDp_task(subj_show))), numel(fp_task), numel(channels))'));
    
    rc = [0.6350 0.0780 0.1840];
    tc =  [0 0.4470 0.7410];
    set(r, 'FaceColor', rc);
    set(t, 'FaceColor', tc);
    set(r, 'FaceAlpha', 0.6)
    set(t, 'FaceAlpha', 0.6)
    set(r, 'EdgeColor', 'none');
    set(t, 'EdgeColor', 'none');

    ylim([0, 30]);
    zlim([0, 300]);
    xticks((1:N_chan)*3-0.25);
    xticklabels(channels);
    xlabel('Channels');
    ylabel('Frequency [Hz]');
    zlabel('Power [\muV^2/Hz]'); %TODO 
    title("PSD");
    legend("rest","task",'ItemHitFcn',@show_hide_graph);
    view(37, 69);

    clear xi yi;
end

%% Resample the signal

fcut = 70;
fs_new = 140;
EEG_rest_rs = struct();
EEG_task_rs = struct();

for s = 1 : N_sub
    for c = 1 : N_chan
        EEG_rest_rs(s).(channels{c}) = resample(subj_rest(s).(channels{c}),fs_new,fs); %basically, take one sample every fs/fs_new
        t_rest_rs = resample(t_rest,fs_new,fs);
        EEG_task_rs(s).(channels{c}) = resample(subj_task(s).(channels{c}),fs_new,fs);
        t_task_rs = resample(t_task,fs_new,fs);        
    end
end

if show_fig
    c_show = chan_show; % redefine chan_show here for this figure

    figure
    hold on
    plot(t_rest, subj_rest(subj_show).(c_show))
    plot(t_rest_rs, EEG_rest_rs(subj_show).(c_show))
    title('EEG '+c_show);
    xlabel("Time");
    xlim([0, max([t_task, t_task_rs])])
    ylabel("Amplitude [\muV]"); %TODO

    legend("original", "resampled",'ItemHitFcn',@show_hide_graph);
end

% uncomment to resave baseline
% tmp = EEG_rest_rs(1).C3;
% save("Data\baseline_resampled.mat", 'tmp');

% No needed anymore
clear subj_rest subj_task t_rest t_task;

%% PSD all resampled signal
lim = [0, fcut];
l_wind_rs = fs_new * 10;
rect_w = rectwin(l_wind_rs);

for c = 1 : N_chan
    [PSDp_rest_rs_1.(channels{c}),fp_rest_rs] = pwelch(EEG_rest_rs(s).(channels{c}),rect_w,noverlap,[],fs_new);
    [PSDp_task_rs_1.(channels{c}),fp_task_rs] = pwelch(EEG_task_rs(s).(channels{c}),rect_w,noverlap,[],fs_new);
end

PSDp_rest_rs = repmat(PSDp_rest_rs_1, 1, N_sub);
PSDp_task_rs = repmat(PSDp_task_rs_1, 1, N_sub);
clear PSDp_rest_rs_1;
clear PSDp_task_rs_1;

for c = 1 : N_chan
    for s = 1 : N_sub
        [PSDp_rest_rs(s).(channels{c}),~] = pwelch(EEG_rest_rs(s).(channels{c}),rect_w,noverlap,[],fs_new);
        [PSDp_task_rs(s).(channels{c}),~] = pwelch(EEG_task_rs(s).(channels{c}),rect_w,noverlap,[],fs_new);
    end
end

if show_fig
    c_show = chan_show; % redefine chan_show here for this figure
    
    figure
    hold on
    xi = repmat((1:N_chan)', 1,numel(fp_rest_rs))*3;
    yi= repmat(fp_task_rs', N_chan, 1);
    
    r = waterfall(xi, yi, (reshape(cell2mat(struct2cell(PSDp_rest_rs(subj_show))), numel(fp_rest_rs), numel(channels))'));
    t = waterfall(xi-0.5, yi, (reshape(cell2mat(struct2cell(PSDp_task_rs(subj_show))), numel(fp_task_rs), numel(channels))'));
    
    rc = [0.6350 0.0780 0.1840];
    tc =  [0 0.4470 0.7410];
    set(r, 'FaceColor', rc);
    set(t, 'FaceColor', tc);
    set(r, 'FaceAlpha', 0.6)
    set(t, 'FaceAlpha', 0.6)
    set(r, 'EdgeColor', 'none');
    set(t, 'EdgeColor', 'none');

    ylim([0, 30]);
    zlim([0, 300]);
    xticks((1:N_chan)*3-0.25);
    xticklabels(channels);
    xlabel('Channels');
    ylabel('Frequency [Hz]');
    zlabel('Power [\muV^2/Hz]'); %TODO
    legend("rest","task",'ItemHitFcn',@show_hide_graph);
    title("PSD of resampled")
    view(37, 69);

    clear xi yi;
end

%% standard bands subdivision

load('filters.mat');

band_coefficients= struct('delta',  filter_delta.Coefficients, ...
    'theta',  filter_theta.Coefficients,...
    'alpha',  filter_alpha.Coefficients,...
    'beta',  filter_beta.Coefficients);

[resample_band_1, resample_band_2] = Band_sub(band_coefficients, fs_new, fs);

%% PDS standard bands
bands = ["delta", "theta", "alpha", "beta"];
lim_delta = [0.5, 4];
lim_theta = [4,  8];
lim_alpha = [8, 13];
lim_beta  = [13, 30];
overlap = floor(l_wind_rs*0.5);
nfft = 2^ceil(log2(l_wind_rs));
len_welch = nfft/2+1;

for c = 1 : N_chan
    PSD_rest.(channels{c}) = zeros(N_sub, numel(bands), len_welch);
    PSD_task.(channels{c}) = zeros(N_sub, numel(bands), len_welch);
end

for c = 1 : N_chan
    for s = 1 : N_sub
        for b = 1 : numel(bands)
            [PSD_rest.(channels{c})(s,b,:),f_rest] = pwelch(resample_band_1.(channels{c}){s,b},rect_w,overlap,[],fs_new);
            [PSD_task.(channels{c})(s,b,:),f_task] = pwelch(resample_band_2.(channels{c}){s,b},rect_w,overlap,[],fs_new);
        end
    end
end

if show_fig
    c_show = chan_show; % redefine chan_show here for this figure

    figure('Name', "PSD "+c_show);
    hold on
    N_b = numel(bands);
    for b = 1 : N_b
        subplot(N_b, 1, b);
        plot(f_rest, squeeze(PSD_rest.(c_show)(subj_show,b,:)));
        hold on
        plot(f_task, squeeze(PSD_task.(c_show)(subj_show,b,:)));
        title('PSD '+bands(b));
        xlabel("Frequency [Hz]");
        ylabel("Amplitude");
        legend("rest", "task",'ItemHitFcn',@show_hide_graph);
    end
    clear N_b;
end


%% selected bands subdivision
load('filters_sub_bands.mat')
band_coefficients= struct('theta1',  filter_theta1.Coefficients, ...
    'theta2',  filter_theta2.Coefficients,...
    'beta1',  filter_beta1.Coefficients,...
    'beta2',  filter_beta2.Coefficients);

[signal_sub_band_1, signal_sub_band_2] = Band_sub(band_coefficients, fs_new, fs);

%% PDS selected bands
sub_bands = ["theta1", "theta2", "beta1", "beta2"];
N_sub_bands= numel(sub_bands);
% lim_theta1 = [4.1, 5.8];
% lim_theta2 = [5.9,  7.4];
% lim_beta1 = [13, 19.9];
% lim_beta2 = [20, 25];

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
            [PSD_rest_theta.(channels{c})(s,i,:), f_t_rest] = pwelch(signal_sub_band_1.(channels{c}){s,i},rect_w,overlap,[],fs_new);
            [PSD_task_theta.(channels{c})(s,i,:), f_t_task] = pwelch(signal_sub_band_2.(channels{c}){s,i},rect_w,overlap,[],fs_new);
            %Beta
            [PSD_rest_beta.(channels{c})(s,i,:), f_b_rest] = pwelch(signal_sub_band_1.(channels{c}){s,i+2},rect_w,overlap,[],fs_new);
            [PSD_task_beta.(channels{c})(s,i,:), f_b_task] = pwelch(signal_sub_band_2.(channels{c}){s,i+2},rect_w,overlap,[],fs_new);
        end
    end
end

if show_fig
    c_show = chan_show; % redefine chan_show here for this figure

    figure('Name', "PSD "+c_show);
    hold on
    for b = 1 : 2
        subplot(N_sub_bands, 1, b);
        plot(f_t_rest, squeeze(PSD_rest_theta.(c_show)(subj_show,b,:)));
        hold on
        plot(f_t_task, squeeze(PSD_task_theta.(c_show)(subj_show,b,:)));
        title('PSD '+sub_bands(b));
        xlabel("Frequency [Hz]");
        ylabel("Amplitude");
        legend("rest", "task",'ItemHitFcn',@show_hide_graph);
    end
    for b = 1 : 2
        subplot(N_sub_bands, 1, b+2);
        plot(f_b_rest, squeeze(PSD_rest_theta.(c_show)(subj_show,b,:)));
        hold on
        plot(f_b_task, squeeze(PSD_task_theta.(c_show)(subj_show,b,:)));
        title('PSD '+sub_bands(b+2));
        xlabel("Frequency [Hz]");
        ylabel("Amplitude");
        legend("rest", "task",'ItemHitFcn',@show_hide_graph);
    end
    clear N_sub_bands;
end

%% Area PSD by channels and selected bands

area_rest = zeros(N_sub, N_chan, N_sub_bands);
area_task = zeros(N_sub, N_chan, N_sub_bands);
for c = 1 : N_chan
    area_rest(:,c,1:2) = sum(PSD_rest_theta.(channels{c}), 3);
    area_rest(:,c,3:4) = sum(PSD_rest_beta.(channels{c}), 3);

    area_task(:,c,1:2) = sum(PSD_task_theta.(channels{c}), 3);
    area_task(:,c,3:4) = sum(PSD_task_beta.(channels{c}), 3);
end

% Median by channels
median_rest = median(area_rest,1);
median_task = median(area_task,1);

%% Topoplot of median area in selected channels and band

% Ordering channels in case are shuffled
[~, idx] = sort(channels);
chanlocs = chanlocs(idx);

for b = 1:N_sub_bands
    figure('Name', sub_bands(b));
    subplot(1,2,1)
    title("Rest")
    topoplot(median_rest(:,:,b), chanlocs, 'electrodes', 'labels', 'style', 'fill');
    subplot(1,2,2)
    title("Task")
    topoplot(median_task(:,:,b), chanlocs, 'electrodes', 'labels', 'style', 'fill');
end


%% Finding variation to choose channels

for b = 1:N_sub_bands
    fprintf("\n\nBand name: %s", sub_bands(b));
    Ratio = abs(((median_task./median_rest)-1).*100);
    tmp=Ratio(:,:,b);
    Ratio(:,isoutlier(tmp),b);
    find(isoutlier(tmp));
    
    Q3 = quantile(tmp, 0.75);
    Q1 = quantile(tmp, 0.25);
    position_high = ((tmp >= Q3));
    position_low = ((tmp <= Q1));
    cph = channels(position_high);
    cpl = channels(position_low);
    fprintf("\nposition high:");
    for i = 1:numel(cph)
       fprintf(" %s", cph{i}); 
    end
    fprintf("\nposition low:");
    for i = 1:numel(cpl)
       fprintf(" %s", cpl{i}); 
    end
    fprintf("\n");
end

%% Selected channels

chan_choosen_names = ["F3", "F8", "FP1", "O2"];
N_chan_ch = numel(chan_choosen_names);
chan_choosen = ones(N_chan_ch)*-1;

for i=1:N_chan_ch
    j = 1;
    while chan_choosen(i) == -1
        if strcmp(channels{j}, chan_choosen_names(i))
           chan_choosen(i) = j; 
        end
        j = j + 1;
    end
end

%% Medians of choosed channels
median_rest_ch = zeros(1, N_chan_ch, N_sub_bands);
median_task_ch = zeros(1, N_chan_ch, N_sub_bands);
for c = 1:N_chan_ch
    median_rest_ch(:, c,:) = median_rest(:,chan_choosen(c),:);
    median_task_ch(:, c,:) = median_task(:,chan_choosen(c),:);
end

%% Save median into a tabel
arraytable = zeros(size(median_rest_ch, 3),size(median_rest_ch, 3)*2);
for b=1:size(median_rest_ch, 3)
    arraytable(:,b*2-1) = median_rest_ch(:,:,b);
    arraytable(:,b*2) = median_task_ch(:,:,b);
end

variablesNames = reshape([sub_bands; sub_bands], 1, 2*N_sub_bands)+repmat({'_rest', '_task'}, 1, N_sub_bands);
tableData = array2table(arraytable, 'VariableNames', variableNames);

writetable(tableData, 'Medians_Table.xlsx')
%% Plot results in scatter

inc = 0.2;
mul = 0.5;
y_lim = 350;
X_ticks = ones(1,N_sub_bands);
X_ticks(1) = X_ticks(1) + inc/2;
for i = 2:N_sub_bands
    X_ticks(i) = X_ticks(i) + (i-1)*mul + inc/2;
end
if show_fig==1
    for c = 1:N_chan_ch
        mul = 0.5;
        ch = chan_choosen(c);
        figure('Name', "Scatter channel "+chan_choosen_names(c))
        hold on
        o_flag = false; % outliers flag 
        for i = 1:N_sub_bands %TODO: figura con i filtri ottennuti(non qui, in generale)
            mul = mul + 0.5;
            r = scatter(mul*ones(1,N_sub), area_rest(:,ch,i),"blue", 'filled');
            outliers_rest = find(area_rest(:,ch,i) > y_lim);
    
            t = scatter((mul+inc)*ones(1,N_sub), area_task(:,ch,i),"red", 'filled');
            outliers_task = find(area_task(:,ch,i) > y_lim); % TODO: va bene calcolati cosi?
            
            m = scatter([mul, mul+inc], [median_rest_ch(:, c,i),median_task_ch(:, c,i)], "d", "MarkerEdgeColor", "magenta");

            for k = 1:numel(outliers_rest)
                o_flag = true;
                o = scatter(mul, y_lim-k*2, "black");
                o.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Actual Y', area_rest(outliers_rest(k),ch,i));
                o.DataTipTemplate.DataTipRows(4) = ' ';
            end
            for k = 1:numel(outliers_task)
                o_flag = true;
                o = scatter(mul+inc, y_lim-k*2, "black");
                o.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Actual Y', area_task(outliers_task(k),ch,i));
                o.DataTipTemplate.DataTipRows(4) = ' ';
            end
        end
        xticks(X_ticks);
        xticklabels({'\theta1', ...
            '\theta2','\beta1', ...
            '\beta2'});
        xlim([X_ticks(1)-inc, X_ticks(numel(X_ticks))+inc]);
        ylim([0, y_lim]);
        if o_flag
            legend([r,t,m,o], "rest", "task", "median", "outliers", "Location", "eastoutside");
        else
            legend([r,t,m], "rest", "task", "median",  "Location", "eastoutside");
        end
    end
end

%% Functions

function show_hide_graph(~,evt)
    if strcmp(evt.Peer.Visible,'on')
        evt.Peer.Visible = 'off';
    else 
        evt.Peer.Visible = 'on';
    end
end
