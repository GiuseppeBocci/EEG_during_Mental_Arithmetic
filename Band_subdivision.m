clear all;
close all;
clc;
path_out="Risultati\Bande di frequenza";
listofsubjects= dir("Data");
lsubj= find(contains({listofsubjects.name}, '_2'));
listofsubjects= listofsubjects(lsubj,:);
load("filters.mat");
% Coefficienti_Bande= struct('theta1',  filter_fir_teta_1, 'theta2',  1:89, 'delta1', 4:20);
coefficienti_Bande= struct('theta',  filter_theta.Coefficients,...
    'alpha',  filter_alpha.Coefficients,...
    'beta',  filter_beta.Coefficients,...
    'delta',  filter_delta.Coefficients);
Fields_Bande= fieldnames(coefficienti_Bande);
Bande_EEG= struct();
%%
for s= 1:numel(listofsubjects)
    EEG= load(['Data/' listofsubjects(s).name]);
            Channels= fieldnames(EEG);
    for ch=1:numel(Channels)

        %BANDS EXTRACTION 
        for banda=1:numel(Fields_Bande)
            % every row has a subject
            EEG_temp{s,banda}= filter(coefficienti_Bande.(Fields_Bande{banda}),1,EEG.(Channels{ch}));
            % TODO: filtfilt
        end
        Bande_EEG.(Channels{ch})=EEG_temp;
    end
    % save(strcat(path_out, '\subject_', num2str(s), '_Bande_baseline.mat'), 'Bande_EEG')
end

