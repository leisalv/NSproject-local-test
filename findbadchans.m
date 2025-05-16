function badchans = findbadchans(EEG, params)
%% [badchans, report] = findbadchans(EEG, params)
%
% Find bad channels
%
% Defaulf values for PREP pipeline
% noisyIn: (fields are filled in on input if not present and propagated to output)
%     robustDeviationThreshold (5) - z score cutoff for robust channel deviation
%     highFrequencyNoiseThreshold (5) -  z score cutoff for SNR (signal above 50 Hz)
%     correlationWindowSeconds (1) - correlation window size in seconds (default = 1 sec)
%     correlationThreshold (0.4) - correlation below which window is bad (default = 0.4)
%     badTimeThreshold (0.01) - cutoff fraction of bad corr windows (default = 0.01)
%     ransacSampleSize (50) - samples for computing ransac (default = 50)
%     ransacChannelFraction (0.25) - fraction of channels for robust reconstruction (default = 0.25)
%     ransacCorrelationThreshold (0.75) - cutoff correlation for abnormal wrt neighbors(default = 0.75)
%     ransacUnbrokenTime (0.4) - cutoff fraction of time channel can have poor ransac predictability (default = 0.4)
%     ransacWindowSeconds (5) - correlation window for ransac (default = 5 sec)
%

fprintf('Bad channel detection...\n');

% Assign parameters vs. defaults
if ~isempty(params.badchandet.robustDeviationThreshold)
    noisyIn.highFrequencyNoiseThreshold = params.badchandet.robustDeviationThreshold;
end
if ~isempty(params.badchandet.highFrequencyNoiseThreshold)
    noisyIn.highFrequencyNoiseThreshold = params.badchandet.highFrequencyNoiseThreshold;
end
if ~isempty(params.badchandet.correlationWindowSeconds)
    noisyIn.correlationWindowSeconds = params.badchandet.correlationWindowSeconds;
end
if ~isempty(params.badchandet.correlationThreshold)
    noisyIn.correlationThreshold = params.badchandet.correlationThreshold;
end
if ~isempty(params.badchandet.badTimeThreshold)
    noisyIn.badTimeThreshold = params.badchandet.badTimeThreshold;
end
if ~isempty(params.badchandet.ransacSampleSize)
    noisyIn.ransacSampleSize = params.badchandet.ransacSampleSize;
end
if ~isempty(params.badchandet.ransacChannelFraction)
    noisyIn.ransacChannelFraction = params.badchandet.ransacChannelFraction;
end
if ~isempty(params.badchandet.ransacCorrelationThreshold)
    noisyIn.ransacCorrelationThreshold = params.badchandet.ransacCorrelationThreshold;
end
if ~isempty(params.badchandet.ransacUnbrokenTime)
    noisyIn.ransacUnbrokenTime = params.badchandet.ransacUnbrokenTime;
end
if ~isempty(params.badchandet.ransacWindowSeconds)
    noisyIn.ransacWindowSeconds = params.badchandet.ransacWindowSeconds;
end

% Re-epoch data to remove period after alarm sound
epochs = length(EEG.epoch);

% Prepare output structure
badchans.all_marked = cell(1,epochs);
badchans.bc_nodata  = cell(1,epochs);
badchans.bc_lowsnr  = cell(1,epochs);
badchans.bc_hfnoise = cell(1,epochs);
badchans.bc_spacorr = cell(1,epochs);
badchans.bc_deviat  = cell(1,epochs);
badchans.bc_ransac  = cell(1,epochs);
badchans.bc_dropout = cell(1,epochs);

if params.badchanthrICA > 0
    badchans.bc_ica     = []; % Add field for bad channels detected with ICA
end

for ne = 1

    signal = pop_select(EEG, 'trial', ne); % Select epoch
    signal = removeTrend(signal); % High pass filter
    %signal = pop_epoch(signal, {'s30'}, [-params.epochStart, 0], 'epochinfo', 'yes');
    signal = pop_reref(signal, []); % Average reference

    %noisyIn.evaluationChannels = params.evaluationChannels;
    %id_ref = find(noisyIn.evaluationChannels==params.recreference);
    %id_veog = find(noisyIn.evaluationChannels==params.veogchannel);
    %id_heog = find(noisyIn.evaluationChannels==params.heogchannel);
    % if strcmpi(EEG.chanlocs(id_veog).labels,'VEOGR')==0
    %    error('Mismatch in chanlocs'); 
    % end
    %noisyIn.evaluationChannels([id_eog, id_ref]) = [];
    noisyStatistics = findNoisyChannels(signal,noisyIn); % Find noisy channels
 
    % Prepare output structure
    badchans.all_marked{ne} = sort(cat(1,noisyStatistics.noisyChannels.all));
    badchans.bc_nodata{ne}  = noisyStatistics.noisyChannels.badChannelsFromNoData;
    badchans.bc_lowsnr{ne}  = noisyStatistics.noisyChannels.badChannelsFromLowSNR;
    badchans.bc_hfnoise{ne} = noisyStatistics.noisyChannels.badChannelsFromHFNoise;
    badchans.bc_spacorr{ne} = noisyStatistics.noisyChannels.badChannelsFromCorrelation;
    badchans.bc_deviat{ne}  = noisyStatistics.noisyChannels.badChannelsFromDeviation;
    badchans.bc_ransac{ne}  = noisyStatistics.noisyChannels.badChannelsFromRansac;
    badchans.bc_dropout{ne} = noisyStatistics.noisyChannels.badChannelsFromDropOuts;

    clear signal;

end; clear ne;

fprintf('Bad channel detection completed\n');
