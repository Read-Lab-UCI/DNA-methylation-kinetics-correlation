clear
close all
tic
setnames={'CGI','Enhancer','SINE'};

for loopset=1:numel(setnames)
    setname=setnames{loopset};
    
    
    %load the input data to get background methylation level and CpG site
    %positions
    
    inputfilename=['../Input_' setname];
    load(inputfilename)
    sites=double(inputsites);
    sites=mod(sites,1E9);
    sites=sites(:);
    fracmeth=inputmeth(:);
    NumCpGSites=numel(sites);
    
    %loop over replicates
    for repnum=1:3
        
        Times=[0,1,4,16]; %Requested time-points in hours
        NumReads_range=[10,5,5,5]; %number of reads requested at each of associated timepoints
        
        AllDat=zeros(NumCpGSites,4,2); %initialize simulated data array. 1st dimension: CpG sites
        %2nd dim: 4 exp timepoints
        %3rd dim: methylated reads (col 1) unmethylated reads (col 2)
        SimName=['SimRepliBSProc_' setname '_' num2str(repnum)];
        datafilename=['AllDat_' SimName];
        for j=1:numel(Times)
            TimePt=Times(j);
            NumReads=NumReads_range(j);
            
            for loopreads=1:NumReads
                TimePts=TimePt+rand; %offset each timepoint in each read w random number btw 0 and 1
                %call the simulator
                [AllDat_each]=Processive_FPTKMC_1DDiff_MultiDomain(sites,fracmeth,TimePts,SimName);
                
                if loopreads==1
                    AllDat_time=AllDat_each;
                else
                    AllDat_time=AllDat_time+AllDat_each;
                end
                SimFlag=1;
            end
            
            AllDat(:,j,:)=AllDat(:,j,:)+reshape(AllDat_time,[NumCpGSites,1,2]);
            
        end
        save(datafilename,'AllDat','sites');
    end
    
    %Simulating the WGBS experiment. Timepoints are randomly chosen
    for repnum=1
        NumReads=20;
        TrajFlag=0;
        SimFlag=1;
        AllDat=zeros(NumCpGSites,2); %initialize simulated data array,
        SimName=['SimWGBSProc_' setname '_' num2str(repnum)];
        datafilename=['AllDat_' SimName];
        for loopreads=1:NumReads
            
            TimePts=rand*30; %randomly selects a post-replication time for this read
            %call the simulator
            [AllDat_each]=Processive_FPTKMC_1DDiff_MultiDomain(sites,fracmeth,TimePts,SimName);
            
            AllDat=AllDat+AllDat_each;
        end
        save(datafilename,'AllDat','sites');
    end
end
toc
