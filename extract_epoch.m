function EEG_eps = extract_epoch(EEG)
    if ~isempty(EEG.epoch)
        EEG_eps = cell(length(EEG.epoch),1);
        for i = 1:length(EEG.epoch)
            EEG_eps{i} = pop_select(EEG, 'trial', i);
        end
    else EEG_eps = cell()
end