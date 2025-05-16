function myicplot(iclist,class,cprob,EEG,goodchannels,outpath)
%% icplot(iclist,chanlocs,goodchannels,varname,outpath)
%
% Plot components
%

nrow = 3; 
ncol = 4;

components = EEG.icawinv;

labels = {'Brain','Muscle','Eye','Heart','LineN','ChanN','Other'};

nfig = ceil(numel(iclist)./(nrow*ncol));
nmod = mod(numel(iclist),(nrow*ncol));

for nf = 1:nfig
    ind2 = nf*(nrow*ncol);
    ind1 = ind2-((nrow*ncol)-1);
    list = ind1:ind2;
    if nf == nfig
        list = list(1:nmod);
    end
    rlis = iclist(list);
    varname = inputname(1);
    figname=fullfile(outpath,['Rep_IC',varname,'_ICs',num2str(nf),'.png']);
    fax=figure('units','normalized','outerposition',[0 0 1 1]);
    for ni = 1:numel(list)
        subplot(nrow,ncol,ni);
        title(['IC',num2str(rlis(ni))]); hold on;
        comp = components(:,rlis(ni));
        topoplot(comp,EEG.chanlocs(goodchannels),'conv','on',...
            'electrodes','on','whitebk','on','style','fill','numcontour',4,'plotrad',0.7);
        text(0, -0.65,[labels{class(rlis(ni))},' ',...
            num2str(cprob(rlis(ni))),'%'],...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end; clear ni;
    exportfig(fax,figname,'Format','png','Color','cmyk',...
        'Resolution',300,'Renderer','opengl');
    clear rlis list ind1 ind2 fax;
end; clear nf;