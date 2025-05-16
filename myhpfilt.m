function [EEG, fakeEEG] = myhpfilt(EEG, hp_freq, mirror_sec)
% MYHIGHPASSFILT performs high-pass filtering with mirror padding to reduce edge artifacts.
%
%   EEG = MYHIGHPASSFILT(EEG, hp_freq, mirror_sec)
%
%   Inputs:
%       EEG         - EEGLAB EEG structure (continuous)
%       hp_freq     - High-pass filter cutoff frequency in Hz (e.g., 0.5)
%       mirror_sec  - Duration in seconds to mirror at beginning (e.g., 5) - end of recording very noisy, will be deleted --> no need to mirror
%
%   Output:
%       EEG         - Filtered EEG structure with mirror region removed

    % Sanity check
    if nargin < 3
        error('Usage: EEG = myhighpassfilt(EEG, hp_freq, mirror_sec)');
    end

    if mirror_sec > 0
        fs = EEG.srate;
        mirror_pts = round(mirror_sec * fs);
    
        % Mirror data
        mirrored_start = fliplr(EEG.data(:, 1:mirror_pts));
        fakeEEG        = EEG;
        fakeEEG.data   = [mirrored_start, EEG.data];
        fakeEEG.pnts   = size(fakeEEG.data, 2);
        fakeEEG.times  = linspace(fakeEEG.xmin, fakeEEG.xmin + (fakeEEG.pnts - 1) / fs, fakeEEG.pnts);
    
        % High-pass filter
        fakeEEG = pop_eegfiltnew(fakeEEG, hp_freq, []);
    
        % Remove mirrored start data
        EEG.data = fakeEEG.data(:, mirror_pts+1:end);
    else
        disp('No timewindow selected for mirroring. Continuing without filtering.');
    end
    EEG = eeg_checkset(EEG);
end