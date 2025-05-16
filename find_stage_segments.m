function segments = find_stage_segments(stage_vec, stage_label)
    is_stage = stage_vec == stage_label;
    d = diff([0, is_stage, 0]);
    start_idx = find(d == 1);
    end_idx = find(d == -1) - 1;
    segments = [start_idx', end_idx'];
end