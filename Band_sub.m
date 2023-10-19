function [bande_EEG_1,bande_EEG_2] = Band_sub(coefficienti_Bande, new_fs, old_fs)
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
            EEG.(Channels{ch})=resample(EEG.(Channels{ch}),new_fs,old_fs);
            for banda=1:numel(fields_bande)
                % every row has a subject
                [dim1,~] = size(coefficienti_Bande.(fields_bande{banda}));
                if dim1 == 1
                    EEG_temp = filter(coefficienti_Bande.(fields_bande{banda}),1,EEG.(Channels{ch}));
                    % TODO: filtfilt
                else
                    EEG_temp = sosfilt(coefficienti_Bande.(fields_bande{banda}),EEG.(Channels{ch}));
                end
                bande_EEG_1.(Channels{ch}){s,banda} = EEG_temp;
            end
        end
    end

    % --- subjects during task ---
    for s = 1:numel(listofsubjects_2)
        
        EEG = load(['Data/' listofsubjects_2(s).name]);
        Channels= fieldnames(EEG);

        for ch = 1:numel(Channels)
            %BANDS EXTRACTION 
            EEG.(Channels{ch})=resample(EEG.(Channels{ch}),new_fs,old_fs);

            for banda=1:numel(fields_bande)
                % every row has a subject
                [dim1,~] = size(coefficienti_Bande.(fields_bande{banda}));
                if dim1 == 1
                    EEG_temp = filter(coefficienti_Bande.(fields_bande{banda}),1,EEG.(Channels{ch}));
                    % TODO: filtfilt
                else
                    EEG_temp = sosfilt(coefficienti_Bande.(fields_bande{banda}),EEG.(Channels{ch}));
                end
                bande_EEG_2.(Channels{ch}){s,banda} = EEG_temp;
            end
  
        end
    end
end
