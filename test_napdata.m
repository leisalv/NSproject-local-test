eeglab;
EEG = pop_loadset('C:\Users\Leila\Local\CessationStudy\Raw\Nap\DA_fullNAP_RAW.set');
pop_saveset(EEG, 'filename', 'napadata_converted.set', 'version', '7.3');