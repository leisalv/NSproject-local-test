function report_badchannels(badchans, params, chanlocs, outpath)
%% report_badchannels(badchannels, badreport, chanlocs, outpath)
%
% Report for bad channels
%


%% Generate maps of bad channel locations

figname=fullfile(outpath,'Rep_BadChans_Topo.png');

fax=figure('units','normalized','outerposition',[0.1 0.10 0.8 0.8]);

badsign = sort(unique([badchans.all_marked{:},badchans.bc_ica]));

map_badsign = zeros(numel(params.evaluationChannels),1);
map_badsign(badsign==1) = 1;

title('Bad electrodes'); hold on;
prc_badsign = (numel(badsign)./numel(params.evaluationChannels)).*100;
topoplot(map_badsign,chanlocs,'maplimits',[0 1],'conv','on',...
    'electrodes','on','emarker2',{badsign,'o',[.2 .2 .9],10,1},...
    'whitebk','on','style','contour','numcontour',1); hold on;
text(0, -0.7,[num2str(round(prc_badsign,1)),'%'], ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
exportfig(fax,figname,'Format','png','Color','cmyk',...
        'Resolution',300,'Renderer','opengl');
clear map_badsign;

hold off;
    

%% Save bad channel list in CSV file

longest = max([numel(sort(unique([badchans.bc_nodata{:}]))),...
    numel(sort(unique([badchans.bc_lowsnr{:}]))),...
    numel(sort(unique([badchans.bc_hfnoise{:}]))),...
    numel(sort(unique([badchans.bc_spacorr{:}]))),...
    numel(sort(unique([badchans.bc_deviat{:}]))),...
    numel(sort(unique([badchans.bc_ransac{:}]))),...
    numel(sort(unique([badchans.bc_dropout{:}]))),...
    numel(sort(unique(badchans.bc_ica)))]);
titles = {'NoData','LowSNR','HFnoise','Spacor','Deviat','Ransac','Dropout','ICAbc'};
reptable = cell(longest+1,length(titles));
reptable(1,:) = titles;
reptable(2:numel(sort(unique([badchans.bc_nodata{:}])))+1,1) = ...
    num2cell(sort(unique([badchans.bc_nodata{:}]))');
reptable(2:numel(sort(unique([badchans.bc_lowsnr{:}])))+1,2) = ...
    num2cell(sort(unique([badchans.bc_lowsnr{:}]))');
reptable(2:numel(sort(unique([badchans.bc_hfnoise{:}])))+1,3) = ...
    num2cell(sort(unique([badchans.bc_hfnoise{:}]))');
reptable(2:numel(sort(unique([badchans.bc_spacorr{:}])))+1,4) = ...
    num2cell(sort(unique([badchans.bc_spacorr{:}]))');
reptable(2:numel(sort(unique([badchans.bc_deviat{:}])))+1,5) = ...
    num2cell(sort(unique([badchans.bc_deviat{:}]))');
reptable(2:numel(sort(unique([badchans.bc_ransac{:}])))+1,6) = ...
    num2cell(sort(unique([badchans.bc_ransac{:}]))');
reptable(2:numel(sort(unique([badchans.bc_dropout{:}])))+1,7) = ...
    num2cell(sort(unique([badchans.bc_dropout{:}]))');
reptable(2:numel(sort(unique([badchans.bc_ica])))+1,8) = ...
    num2cell(sort(unique([badchans.bc_ica]))');

writecell(reptable,fullfile(outpath,'Rep_BadChans_Table.csv'))


%% Save bad channel list in text file

filename=fullfile(outpath,'Rep_BadChans_List.txt');
fid=fopen(filename,'w');
for nx = 1:length(badsign)
    fprintf(fid, [num2str(badsign(nx))  '\n']);
end
fclose(fid);