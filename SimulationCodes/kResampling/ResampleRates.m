clear
close all
%this script loads in k-rates along with various feature data for a
%specified region. After creating a site-matched feature array, a
%resampling procedure is performed consisting of two steps: 1. unbiased
%kmeans on the feature data (NOT on k). 2. Randomly shuffling CpG positions
%within each kmeans-generated-cluster. The new dataset, with k-rates
%reassigned to new "positions" is then saved.
RegionIDs={'CGI'};

for ii=1:numel(RegionIDs)
    setname=RegionIDs{ii};
    ratefilename=['kineticRates_' setname];
    load(ratefilename) %load the datafile
    MaxN=min(200000,numel(fittedSites));
    inds=1:MaxN;
    keepWGBS=WGBS(1:MaxN);
    keepDNase=DNase(1:MaxN);
    keepCpGd=CpGd(1:MaxN);
    keepsites=fittedSites(1:MaxN);
    Features=[keepWGBS(:),double(keepCpGd(:)),keepDNase(:)];
    
    %Not all sites have WGBS data associated. Perform smoothing to reduce
    %noise and fill in the gaps. i.e. calc avg quantities in
    %sliding window
    windowSz=200; %bp
    chunksize=round(windowSz/2);
    tileddata=zeros(size(Features));
    NumCpGSites=numel(inds);
    for loopsites=1:NumCpGSites
        starti=max(1,loopsites-chunksize);
        endi=min(NumCpGSites,loopsites+chunksize);
        chunkinds=starti:endi;
        chunk=keepsites(chunkinds);
        getinds=find(abs(keepsites(loopsites)-chunk)<=windowSz);
        WGBSchunk=Features(chunkinds(getinds),1);
        tileddata(loopsites,1)=mean(WGBSchunk(WGBSchunk>=0));
        tileddata(loopsites,2:3)=mean(Features(chunkinds(getinds),2:3),1,'omitnan');
    end
    Features=tileddata;
    
    %filter out the sites that still don't have WGBS data after smoothing,
    %and truncate (we can't use them)
    keepindsWGBS=find(~isnan(Features(:,1)));
    keepks=inferedRates(inds(keepindsWGBS));
    keepsites=keepsites(keepindsWGBS);
    
    Features=Features(keepindsWGBS,:);
    NumCpGSites=numel(keepks);
    
    Features=Features./repmat(std(Features,'omitnan'),NumCpGSites,1); %rescale the feature space by std
    %now resample the ks based on their local features. Output an
    %array of resampled ks (Resample_ks) that is NumCpGSites rows (same number of sites as
    %previous) and NSamp columns (to do it multiple times to get better
    %sampling)
    NSamp=3;
    Resample_ks=zeros(NumCpGSites,NSamp);
    for loopsamp=1:NSamp
        NClust=60;
        idx=kmeans(Features,NClust,'MaxIter',10000);
        for jj=1:NClust
            getinds=find(idx==jj);
            kchunk=keepks(getinds);
            NMembers=numel(getinds);
            Perm=randperm(NMembers);
            newkchunk=kchunk(Perm);
            Resample_ks(getinds,loopsamp)=newkchunk;
        end
    end
    
    keepks=keepks(:);
    keepsites=keepsites(:);
    savefilename=[setname '_resampled'];
    save(savefilename,'keepsites','keepks','Resample_ks')   

end



