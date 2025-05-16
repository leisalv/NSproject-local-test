% EEGLAB history file generated on the 29-Apr-2025
% ------------------------------------------------
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename','myeeg_superclean.set','filepath','C:\\Users\\Leila\\Local\\NSproject\\data\\derivatives\\NS_nap_preproc_v5_2025-04-29_19-52-09\\prep_data\\');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
pop_saveh( EEG.history, 'eeglabhistory.m', 'C:\Users\Leila\OneDrive - Vrije Universiteit Amsterdam\Desktop\');
eeglab redraw;
