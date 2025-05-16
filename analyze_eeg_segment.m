function [ORIGINAL, OSCILLATORY, SLOPES, INTERCEPTS] = analyze_eeg_segment(EEG_seg, freq_bands, slope_ranges)

    % Convert to FieldTrip structure
    eegdata = eeglab2fieldtrip(EEG_seg, 'preprocessing', 'none');

    % Detrend
    cfg = []; cfg.detrend = 'yes';
    eegdata = ft_preprocessing(cfg, eegdata);

    % IRASA decomposition
    cfg = [];
    cfg.length = 5; cfg.overlap = 0;
    eegdata_epc = ft_redefinetrial(cfg, eegdata);

    cfg = [];
    cfg.method = 'irasa'; cfg.foilim = [0.5 50*1.9]; 
    cfg.pad = 'nextpow2'; cfg.taper = 'hanning';
    cfg.keeptrials = 'yes';

    cfg.output = 'original'; original = ft_freqanalysis(cfg, eegdata_epc);
    cfg.output = 'fractal'; fractal = ft_freqanalysis(cfg, eegdata_epc);

    cfg_math = []; cfg_math.parameter='powspctrm'; cfg_math.operation='x2 ./ x1';
    oscillatory = ft_math(cfg_math, fractal, original);

    % Extract powers and slopes
    % bandpower = NaN(nchans, size(freq_bands,1));
    % slopes = NaN(nchans,size(slope_ranges,1));
    
    ntrls = size(original.powspctrm,1);
    nchans = size(original.powspctrm,2);
    freqs = original.freq;

    ORIGINAL = NaN(ntrls, nchans, size(freq_bands,1));
    OSCILLATORY = NaN(ntrls, nchans, size(freq_bands,1));
    SLOPES = NaN(ntrls, nchans, size(slope_ranges,1));
    INTERCEPTS = NaN(ntrls, nchans, size(slope_ranges,1));

    for t = 1:ntrls
        
        fprintf('Loading trial %d / %d\n', t, trls);

        osc_pow = squeeze(oscillatory.powspctrm(t,:,:));
        orig_pow = squeeze(original.powspctrm(t,:,:));
        frac_pow = squeeze(fractal.powspctrm(t,:,:));

        for ch = 1:nchans

            fprintf('Channel #%d', ch);

            for b = 1:size(freq_bands,1)
                
                fprintf('Frequency band #%d', freq_bands(b,:));

                b_idx = freqs >= freq_bands(b,1) & freqs < freq_bands(b,2);
                ORIGINAL(t,ch,b) = trapz(freqs(b_idx), orig_pow(ch, b_idx));
                OSCILLATORY(t,ch,b) = trapz(freqs(b_idx), osc_pow(ch, b_idx));

            end
    
            for s = 1:size(slope_ranges,1)

                fprintf('Slope range #%d', slope_ranges(s,:));

                % slope estimation
                s_idx = freqs >= slope_ranges(s,1) & freqs <= slope_ranges(s,2);
                [SLOPES(t,ch,s), INTERCEPTS(t,ch,s)] = logfit(freqs(s_idx), frac_pow(ch, s_idx), 'loglog');
            end

        end
    end
end