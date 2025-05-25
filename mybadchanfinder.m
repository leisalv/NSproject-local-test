function badchans = mybadchanfinder(EEG, params)
%% [badchans, report] = findbadchans(EEG, params)
%
% Find bad channels / adapted to coontinuous, non-epoched data
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

% Prepare output structure
badchans.all_marked = cell(1);
badchans.bc_nodata  = cell(1);
badchans.bc_lowsnr  = cell(1);
badchans.bc_hfnoise = cell(1);
badchans.bc_spacorr = cell(1);
badchans.bc_deviat  = cell(1);
badchans.bc_ransac  = cell(1);
badchans.bc_dropout = cell(1);

if params.badchanthrICA > 0
    badchans.bc_ica     = []; % Add field for bad channels detected with ICA
end

signal = EEG; % Select eeg data
signal = removeTrend(signal); % High pass filter
signal = pop_reref(signal, []); % Average reference

noisyIn.evaluationChannels = params.evaluationChannels;

% Remove reference and auxiliary channels from evaluation channels if present
channelFields = {'recref', 'veogchannel', 'heogchannel', 'emgchannel'}; %'M1channel', 'M2channel', 

for i = 1:length(channelFields)
    chan = params.(channelFields{i});
    if ~isempty(chan) && ismember(chan, noisyIn.evaluationChannels)
        noisyIn.evaluationChannels(noisyIn.evaluationChannels == chan) = [];
    end
end

noisyStatistics = findNoisyChannels(signal, noisyIn); % Find noisy channels
 
% Prepare output structure
badchans.all_marked = sort(cat(1,noisyStatistics.noisyChannels.all));
badchans.bc_nodata  = noisyStatistics.noisyChannels.badChannelsFromNoData;
badchans.bc_lowsnr  = noisyStatistics.noisyChannels.badChannelsFromLowSNR;
badchans.bc_hfnoise = noisyStatistics.noisyChannels.badChannelsFromHFNoise;
badchans.bc_spacorr = noisyStatistics.noisyChannels.badChannelsFromCorrelation;
badchans.bc_deviat  = noisyStatistics.noisyChannels.badChannelsFromDeviation;
badchans.bc_ransac  = noisyStatistics.noisyChannels.badChannelsFromRansac;
badchans.bc_dropout = noisyStatistics.noisyChannels.badChannelsFromDropOuts;

clear signal;

fprintf('Bad channel detection completed\n');
