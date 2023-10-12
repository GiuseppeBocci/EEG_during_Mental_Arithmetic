clear; close all; clc;

% import data
path_out = "Risultati\Bande di frequenza";
listofsubjects = dir("Data");
lsubj_1 = find(contains({listofsubjects.name}, '_1'));
listofsubjects_1 = listofsubjects(lsubj_1,:);
lsubj_2 = find(contains({listofsubjects.name}, '_2'));
listofsubjects_2 = listofsubjects(lsubj_2,:);

load("filters.mat");

coefficienti_Bande= struct('delta',  filter_delta.Coefficients, ...
    'theta',  filter_theta.Coefficients,...
    'alpha',  filter_alpha.Coefficients,...
    'beta',  filter_beta.Coefficients);
fields_Bande= fieldnames(coefficienti_Bande);

bande_EEG= struct();


% filter every signal and channels in each band
% --- subjects in rest ---
for s= 1:numel(listofsubjects_1)
    EEG= load(['Data/' listofsubjects(s).name]);
            Channels= fieldnames(EEG);

    for ch=1:numel(Channels)
        %BANDS EXTRACTION 
        for banda=1:numel(fields_Bande)
            % every row has a subject
            EEG_temp{s,banda}= filter(coefficienti_Bande.(fields_Bande{banda}),1,EEG.(Channels{ch}));
            % TODO: filtfilt
        end
        bande_EEG.(Channels{ch})=EEG_temp;
    end
    % save(strcat(path_out, '\subject_', num2str(s), '_Bande_baseline.mat'), 'Bande_EEG')
end

