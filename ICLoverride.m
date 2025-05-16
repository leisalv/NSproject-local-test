function [reject, keep] = ICLoverride(reject, keep, outpath)
%% [reject, keep] = ICLoverride(reject, keep, outsets, outpath)

ic_keep=[]; ic_rejc=[];
if exist([outpath,'ICLoverride_keep.txt'],'file')>0
    ic_keep = csvread([outsets,'ICLoverride_keep.txt']);
    warning('Overriding ICLabel results: re-adding rejected ICs!');
    [~,id] = intersect(reject, ic_keep);
    keep = sort(unique([keep; ic_keep]));
    reject(id) = []; clear id;
end
if exist([outpath,'ICLoverride_rejc.txt'],'file')>0
    ic_rejc = csvread([outpath,'ICLoverride_rejc.txt']);
    warning('Overriding ICLabel results: rejecting included ICs!');
    [~,id] = intersect(keep, ic_rejc);
    reject = sort(unique([reject; ic_rejc]));
    keep(id) = []; clear id;
end

if ~isempty(ic_keep) || ~isempty(ic_rejc)
    writemod = sort([ic_keep; ic_rejc]);
    filename=fullfile(outpath,'Rep_ICLoverride.txt');
    fid=fopen(filename,'w');
    for nx = 1:length(writemod)
        fprintf(fid, [num2str(writemod(nx))  '\n']);
    end; clear nx;
end