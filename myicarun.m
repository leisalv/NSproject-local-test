function EEG = myicarun(EEG, goodchannels, params)
%% EEG = icarun(EEG,goodchannels,params)
%
% Remove bad components

pcan = params.applyPCAtoICA; % PCA before ICA
resm = params.resampleICA; % Frequency for resampling

if params.lowPassFilterICA>0 % Low-pass
    EEG = pop_eegfiltnew(EEG, [], params.lowPassFilterICA, [], 0, [], 0);  % Low pass filter
end

if params.resampleICA>0 % Resample
    EEG = pop_resample(EEG, resm); % Resample data if requested
end 

if params.highPassFilterICA>0  % High pass filter - comment out if doing highpass filtering before ASR
    EEG = pop_eegfiltnew(EEG, params.highPassFilterICA, [], [], 0, [], 0); 
end

if pcan==0
    EEG = pop_runica(EEG,'icatype','runica',...
        'chanind',goodchannels,'rndreset','no');
else
    EEG = pop_runica(EEG,'icatype','runica',...
        'chanind',goodchannels,'pca', pcan,'rndreset','no');
end

EEG = eeg_checkset(EEG);  % Check data structure