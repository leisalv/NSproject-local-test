function EEG = icarun(EEG,goodchannels,params)
%% EEG = icarun(EEG,goodchannels,params)
%
% Remove bad components
%

pcan = params.applyPCAtoICA; % PCA before ICA
resm = params.resampleICA; % Frequency for resampling

if params.lowPassFilter>0 % Low-pass
    EEG = pop_eegfiltnew(EEG, [], params.lowPassFilter, [], 0, [], 0);  % Low pass filter
end

if params.resampleICA>0 % Resample
    EEG = pop_resample(EEG, resm); % Resample data if requested
end 

if params.highPassFilter>0 % High-pass
    EEG = pop_eegfiltnew(EEG, params.highPassFilter, [], [], 0, [], 0); % High pass filter
end

%EEG = pop_epoch( EEG, {'s30'}, [-params.epochStart, 0], 'epochinfo', 'yes');

if pcan==0
    EEG = pop_runica(EEG,'icatype','runica',...
        'chanind',goodchannels,'rndreset','no');
else
    EEG = pop_runica(EEG,'icatype','runica',...
        'chanind',goodchannels,'pca',pcan,'rndreset','no');
end

EEG = eeg_checkset(EEG);  % Check data structure