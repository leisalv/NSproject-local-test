function badchans = mybadchannelfinder(EEG, params)
%% [badchans, report] = findbadchans(EEG, params)
%
% Find bad channels / adapted to coontinuous, non-epoched data / keeping only default params
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
% signal = removeTrend(signal); % High pass filter - commented out if FIR applied before
signal = pop_reref(signal, []); % Average reference

noisyStatistics = findNoisyChannels(signal); % Find noisy channels
 
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
