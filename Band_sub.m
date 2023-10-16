function [bande_EEG_1,bande_EEG_2] = Band_sub(coefficienti_Bande)
    % import data
    listofsubjects = dir("Data");
    lsubj_1 = find(contains({listofsubjects.name}, '_1'));
    listofsubjects_1 = listofsubjects(lsubj_1,:);
    lsubj_2 = find(contains({listofsubjects.name}, '_2'));
    listofsubjects_2 = listofsubjects(lsubj_2,:);

    % initialization of variables
    bande_EEG_1 = struct();
    bande_EEG_2 = struct();
    fields_bande= fieldnames(coefficienti_Bande);
    
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
end