function [ep_evt_str, ep_evt_end, ep_evt_szs] = eventcsv2mat(evtfile, srate, evtype)

events = table2cell(evtfile);  
events = events(strcmp(events(:,7), evtype), :);
ep_evt_str = ([events{:,2}].*srate)+1;  % beginning of each artefact
ep_evt_end = ([events{:,3}].*srate);    % termination of each artefact
ep_evt_szs = (ep_evt_end-ep_evt_str+1);

clear events;