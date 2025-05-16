function report_iccleaning(EEG,reject,keep,classcomps,goodchannels,outpath)
%% report_iccleaning(EEG,reject,keep,classcomps,goodchannels,outpath)
%
% Plot components
%

[~, class]  = max(classcomps,[],2);
cprob = NaN(length(class),1);
for ni = 1:length(class)
    cprob(ni,1) = round((classcomps(ni,class(ni))).*100,1);
end; clear ni;

% Plot Rejected components
icplotrest_prova(reject,class,cprob,EEG,goodchannels,outpath)

% Plot Retained components
icplotrest_prova(keep,class,cprob,EEG,goodchannels,outpath)


%% Save bad comps list in text file

filename=fullfile(outpath,'Rep_BadComps_List.txt');
fid=fopen(filename,'w');
for nx = 1:length(reject)
    fprintf(fid, [num2str(reject(nx))  '\n']);
end
fclose(fid);


%% Save good comps list in text file

filename=fullfile(outpath,'Rep_GoodComps_List.txt');
fid=fopen(filename,'w');
for nx = 1:length(keep)
    fprintf(fid, [num2str(keep(nx))  '\n']);
end
fclose(fid);


%% Report bad component  

classes = EEG.etc.ic_classification.ICLabel.classes;
titles = [{'NICs','RejectedICs'}, classes(2:end)];
Ncomp = size(classcomps,1);
reptable = cell(Ncomp+1,length(titles));
reptable(1,:) = titles;
reptable(2,1) = num2cell(Ncomp);
[~,typeic]= max(classcomps(:,2:end),[],2);
typeic= typeic+1;

if ~isempty(reject)
    reptable(2:length(reject)+1,2) = num2cell(reject);    
    reptable(2:length(setdiff(find(typeic==2),keep))+1,3) = num2cell(setdiff(find(typeic==2),keep));
    reptable(2:length(setdiff(find(typeic==3),keep))+1,4) = num2cell(setdiff(find(typeic==3),keep));
    reptable(2:length(setdiff(find(typeic==4),keep))+1,5) = num2cell(setdiff(find(typeic==4),keep));
    reptable(2:length(setdiff(find(typeic==5),keep))+1,6) = num2cell(setdiff(find(typeic==5),keep));
    reptable(2:length(setdiff(find(typeic==6),keep))+1,7) = num2cell(setdiff(find(typeic==6),keep));
    reptable(2:length(setdiff(find(typeic==7),keep))+1,8) = num2cell(setdiff(find(typeic==7),keep));
end 
writecell(reptable,fullfile(outpath,'Rep_ICreject_Table.csv'))