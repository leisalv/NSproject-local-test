function EEG_segs = extract_eeg_segments(EEG, segments)
    EEG_segs = cell(size(segments,1),1);
    for i = 1:size(segments,1)
        EEG_segs{i} = pop_select(EEG, 'point', segments(i,:));
        EEG_segs{i} = eeg_checkset(EEG_segs{i});
    end
end