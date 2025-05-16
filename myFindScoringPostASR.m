% Step 1: Copy original_labels to a new variable
cleaned_labels = original_labels;

% Step 2: Extract boundary events only
boundary_events = myeeg_prep.event(strcmp({myeeg_prep.event.type}, 'boundary'));

% Step 3: Initialize variable to track total removed points so far
total_removed = 0;

% Step 4: Loop through each boundary sequentially
for i = 1:length(boundary_events)
    % Get latency and duration
    latency = boundary_events(i).latency;
    duration = boundary_events(i).duration;

    % Adjust latency for previously removed samples
    adjusted_start = ceil(latency - total_removed);
    adjusted_end = adjusted_start + duration - 1;

    % Remove the specified segment
    cleaned_labels(adjusted_start:adjusted_end) = [];

    % Update how many total points have been removed
    total_removed = total_removed + (adjusted_end - adjusted_start + 1);
end
