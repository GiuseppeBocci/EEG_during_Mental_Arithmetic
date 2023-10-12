clear; close all; clc;

% import data
path_out = "bands";
listofsubjects = dir("Data");
lsubj_1 = find(contains({listofsubjects.name}, '_1'));
listofsubjects_1 = listofsubjects(lsubj_1,:);
lsubj_2 = find(contains({listofsubjects.name}, '_2'));
listofsubjects_2 = listofsubjects(lsubj_2,:);

load("filters.mat");

% initialization of variables
coefficienti_Bande= struct('delta',  filter_delta.Coefficients, ...
    'theta',  filter_theta.Coefficients,...
    'alpha',  filter_alpha.Coefficients,...
    'beta',  filter_beta.Coefficients);
fields_bande= fieldnames(coefficienti_Bande);

bande_EEG_1 = struct();
bande_EEG_2 = struct();


% filter every signal and channels in each band

% --- subjects in rest ---
for s = 1:numel(listofsubjects_1)
    
    EEG = load(['Data/' listofsubjects_1(s).name]);
            Channels= fieldnames(EEG);

    for ch = 1:numel(Channels)
        %BANDS EXTRACTION 
        for banda=1:numel(fields_bande)
            % every row has a subject
            EEG_temp{s,banda}= filter(coefficienti_Bande.(fields_bande{banda}),1,EEG.(Channels{ch}));
            % TODO: filtfilt
        end
        bande_EEG_1.(Channels{ch})=EEG_temp;
    end
end

save('original_band_1.mat',"bande_EEG_1");

% --- subjects during task ---
for s = 1:numel(listofsubjects_2)
    
    EEG = load(['Data/' listofsubjects_2(s).name]);
            Channels= fieldnames(EEG);

    for ch = 1:numel(Channels)
        %BANDS EXTRACTION 
        for banda=1:numel(fields_bande)
            % every row has a subject
            EEG_temp{s,banda}= filter(coefficienti_Bande.(fields_bande{banda}),1,EEG.(Channels{ch}));
            % TODO: filtfilt
        end
        bande_EEG_2.(Channels{ch})=EEG_temp;
    end
end

save('original_band_2.mat',"bande_EEG_2");

