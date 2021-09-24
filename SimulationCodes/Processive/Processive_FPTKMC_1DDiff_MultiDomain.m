function [SimDat] = Processive_FPTKMC_1DDiff_MultiDomain(sites,fracmeth,TimePt,SimName)

%Simulates post-replication methylation according to Processive Model using
%Kinetic Monte Carlo First Passage Time algorithm

%Inputs: 
%sites: position IDs of CpGs (Nx1 array, N=NumCpGSites)
%fracmeth: the steady-state fraction methylated at each site: these are used to intialize the time-0 state of DNA at start of simulation
        %(Nx1 array)
%TimePts is the value of time (in hours) at which state of system is requested 
        %will get as close as possible to requested time with Gillespie
        %algorithm which is next-event
%SimName: (string) Name of simulation, only used for saving parameters file

%Outputs:
%SimDat: (N x 2 array), 1st column is 1 if the site is methylated at the requested timepoint, 
        %0 otherwise (and the second column is opposite
    
NumCpGSites=numel(sites);
S0=round(sum(fracmeth)); %assumes that the number of substrate is only based on hemimethylated sites,
%not total CpGs
Ratio=0.022; %Ratio of enzyme molecules to hemi-CpGs in simulation. 
k1r=5; %units %hr^-1
kcat=40; %units hr^-1 

Km=0.8; %micromolar, 
k1f_det=(k1r+kcat)/Km; %in %hr^-1 muM^-1
Vnuc=1E-12; %liters (nuclear volume, used for converting micromolar to copy number)
htot=30E6; %total number of hCpG substrates in the nucleus
NA=6.022E23;
V_rxn=Vnuc*NumCpGSites/htot; %well-mixed reaction volume
k1f=k1f_det/1E-6/V_rxn/NA; %convert the deterministic k1f to stochastic for this volume

koff=1E4; %units /h. (value shouldn't matter --only ratio with D

D_k_ratio=1/(0.0267)^2; %ratio D/koff

ParamFileName=['Params_' SimName];
ParamArray=[Ratio,k1f,k1r,kcat,koff,D_k_ratio,S0];
save(ParamFileName,'ParamArray')

E=max(1,round(Ratio*S0)); %total number of Enzymes in simulation
D=D_k_ratio*koff; %Diffusion coeff, based on ratio set by user

%calc predicted mean half time
Km_stoch=(k1r+kcat)/k1f;
t50=(Km_stoch+S0)/(2*kcat*E);

%We need lookup tables for properties of the first passage time
%distributions, so we don't have to recalculate them repeatedly during the
%simulation, which could be costly. These need to be calculated once per D_k_ratio
LookupFName=['Dratio' num2str(D_k_ratio) '.mat'];
CDFName=['CDFs' LookupFName];
FPT_tName=['FPT_ts' LookupFName];
ExitName=['Exit_Probs' LookupFName];
MaxDistName=['MaxDist' LookupFName];
if exist(CDFName,'file')==2
    load(CDFName,'CDFs');
    load(FPT_tName,'FPT_ts');
    load(ExitName,'Exit_Probs');
    load(MaxDistName,'MaxDist');
else
    % if the files are not found, call the function that computes them
    [CDFs,FPT_ts,Exit_Probs,MaxDist]=MakeFPTLookup_analytical(D_k_ratio);
    disp('table done')
    save(CDFName,'CDFs','-v7.3')
    save(FPT_tName,'FPT_ts','-v7.3')
    save(ExitName,'Exit_Probs','-v7.3')
    save(MaxDistName,'MaxDist');
end

%Reactions
%E+h <-> Eh k1f, k1r %enzyme binding/unbinding reactions
%Eh -> Em kcat %catalytic step

%Diffusion/unbinding reactions
%Em + h -> m + Eh (diffusion) %The enzyme leaves current site and reaches next
%target by diffusion
%Em -> E + m koff %The enzyme unbinds to solution

NumRxns=3;  %number of reactions (only for Gillespie rxns, not diffusion)
NumSpec=5; % [u,h,Eh,Em,m] %species in the model

SimDat=zeros(NumCpGSites,2); %initialize array to store data, cols: methylated, unmeth.

%setting up reaction stoichiometries for Kinetic Monte Carlo (Gillespie
%part). (The diffusive superhops are taken care of elsewhere)
Stoich=zeros(NumRxns,NumSpec); %stoichiometry array
%Rxn 1; h + E -> Eh, k1f
Stoich(1,2)=-1;
Stoich(1,3)=1;
%Rxn 2; Eh -> h + E, k1r
Stoich(2,2)=1;
Stoich(2,3)=-1;
%Rxn 3; Eh -> Em kcat
Stoich(3,3)=-1;
Stoich(3,4)=1;


%Initialize the starting state with some hemimethylated sites,
%according to fracmeth
SitesArray0=zeros(NumCpGSites,NumSpec);
InitRandArray=rand(NumCpGSites,1);
SitesArray0(InitRandArray<=fracmeth(:),2)=1;
SitesArray0(:,1)=abs(SitesArray0(:,2)-1);
CurrState=SitesArray0;
clear SitesArray0 InitRandArray

%call the function that computes propensities for KMC steps. a is array with propensity of each reaction
%(columns) at each site (rows)
[a,a_0]=GetPropensity(CurrState);

%Initialize arrays that keep track of the diffusion domains for the FPT
%steps.
DiffusionDomains=zeros(NumCpGSites,3); %cols=Original Em position, L neighbor position, R neighbor position
DiffusionTimeTicker=zeros(NumCpGSites,3); %cols = sampled exit time, Original Em index (leaving index), Target site index (where it exited to)
DiffusionTimeTicker(:,1)=Inf; %initially set the exit times to infinity. (They will never happen.)
tau=0;
count=0;

%call the function that updates the lists based on current system state
[DiffusionDomains,DiffusionTimeTicker]=UpdateDiffusionLists(CurrState,sites,DiffusionDomains,DiffusionTimeTicker,tau,count);

Time=0;

while Time<TimePt
    [tau_KMC,I,J]=KMCStep(a,a_0);
    [tau_FPT,sInd,nInd]=FPTStep(DiffusionTimeTicker);
    
    %the next event (KMC or FPT) corresponds to min tau
    if tau_KMC<=tau_FPT %if time to KMC step is shorter, select it
        tau=tau_KMC;
        count=count+1;%count the number of steps
        if isinf(tau) %if no reaction is possible, increment the time anyway
            Time=Time+.5;
        else
            Time=Time+tau;
            CurrState(I,:)=CurrState(I,:)+Stoich(J,:); %update system state
        end
    else %else the next step is a diffusion/unbinding event
        tau=tau_FPT;
        count=count+1;
        if isinf(tau) %if no reaction is possible, increment the time anyway
            Time=Time+.5;
        else
            Time=Time+tau;
            %update the system state to reflect diffusion from an Em site to neighbor h
            %change the state of the start and stop sites for hop events
            CurrState(sInd,4)= CurrState(sInd,4)-1;
            CurrState(sInd,5)= CurrState(sInd,5)+1;
            if nInd>0
                %if nID=0, then the enzyme goes back to solution, no change to neighbor state
                CurrState(nInd,2)= CurrState(nInd,2)-1;
                CurrState(nInd,3)= CurrState(nInd,3)+1;
            end
        end
    end
    
    %Reset the propensity array and NeighborArray according to new state of system
    [a,a_0]=GetPropensity(CurrState); %call the function that computes the propensity for each site
    %call the function that updates the diffusion information: e.g. domains
    [DiffusionDomains,DiffusionTimeTicker]=UpdateDiffusionLists(CurrState,sites,DiffusionDomains,DiffusionTimeTicker,tau,count);

end
Time
SimDat(:,1)=SimDat(:,1)+CurrState(:,4)+CurrState(:,5); %store methylated at desired timepoint
SimDat(:,2)=SimDat(:,2)+CurrState(:,1)+CurrState(:,2)+CurrState(:,3); %store unmethylated at desired timepoint

    function [a,a_0]=GetPropensity(CurrState)
        %Calc the per-site propensities for each reaction. Array is
        %NumCpGs x NumRxns. This only computes propensities for Gillespie
        %reactions. Diffusion/unbinding are taken care of in the FPTStep
        
        %Inputs: CurrState array (NumCpGSitesxNumSpec)
        %Outputs: 
            % a: (NumCpGSites x NumSpec) the propensities of each reaction
                    % at each site
            % a0: (double) sum over all propensities
        Enum=E-sum(CurrState(:,3))-sum(CurrState(:,4)); %number of free enzyme copies
        %Calc the per-site propensities for each reaction.
        a=zeros(NumCpGSites,NumRxns);
        a(:,1)=k1f*Enum*CurrState(:,2); %E + h -> [Eh]
        a(:,2)=k1r*CurrState(:,3); %  [Eh] -> E + h
        a(:,3)=kcat*CurrState(:,3); %Eh -> Em
        a_0=sum(a(:));
    end

    function [tau_KMC,I,J]=KMCStep(a,a_0)
        %Inputs: a and a_0
        %Outputs: tau_KMC: time to next reaction
            %I - row (site) at which next reaction occurs
            %J - column (index of which reaction)
        %Gillespie algorithm: 2 random numbers, first to get time to next
        %reaction, 2nd to select which reaction
        Rands_KMC=rand(1,2); %carry out Gillespie step, i.e. Kinetic Monte Carlo
        tau_KMC=1/a_0*log(1/Rands_KMC(1)); %time to next reaction
        findnextrxn=cumsum(a(:))>Rands_KMC(2)*a_0;
        NextRxn=sum(1-findnextrxn)+1; %index of next reaction
        [I,J]=ind2sub(size(a),NextRxn); %I returns the row index where the next
        %reaction occurs, J is the index of next reaction from NumRxns
    end

    function [tau_FPT,sInd,nInd]=FPTStep(DiffusionTimeTicker)
        %select the next FPT event by choosing the minimum time, select the
        %miminum FPT and get the target site 
        [tau_FPT,samp_ind]=min(DiffusionTimeTicker(:,1));
        sInd=DiffusionTimeTicker(samp_ind,2); %leaving site index
        nInd=DiffusionTimeTicker(samp_ind,3); %target site index
    end

    function [DiffusionDomains,DiffusionTimeTicker]=UpdateDiffusionLists(CurrState,sites,DiffusionDomains,DiffusionTimeTicker,tau,count)
        %this function updates the DiffusionDomains and DiffusionTimeTicker
        %arrays based on current state of the system. It looks for places
        %where a state change has impacted the diffusion domains. (Diffusion
        %domains need to be kept track of for multiple timesteps because
        %diffusion is not memoryless!--i.e. can't just consider current state of system,
        %but also past history). Then it updates the ticker accordingly and
        %returns both updated lists.
        
        %subtract the latest time step from the stored times in
        %TimeTicker col 1. This is because the clock is still running for
        %domains that remain in play (an enzyme is still bound, and has
        %neither reached target or unbound to solution yet--but other
        %reactions have occured in the system so time has elapsed.)
        DiffusionTimeTicker(:,1)=DiffusionTimeTicker(:,1)-tau;
        
        %find indices of Em states currently and compare to previous
        %timestep
        Curr_Em_inds_all=find(CurrState(:,4)); %current Em
        Prev_Em_inds_all=find(DiffusionDomains(:,1)); %previous Em
        [C,ia,ib]=setxor(Curr_Em_inds_all,Prev_Em_inds_all); %get the indices of sites that have appeared or disappeared
        %ia are indices of sites that are newly Em, and thus need to be
        %added for the first time to DiffusionTimeTicker and
        %DiffusionDomains.
        %ib are indices of sites that are no longer Em and need to be
        %removed from those arrays.
        
        NewDiffusionDomains=zeros(NumCpGSites,3); %reinitialize the array for new domains
        if Curr_Em_inds_all
            [ActiveDiffDom,ActiveDiffDom_I] = BuildDiffusionDomains(CurrState,sites); %build the domains for all the current Em states
            NewDiffusionDomains(Curr_Em_inds_all,:)=ActiveDiffDom; %populate the new list based on current state
            if ia %if some of the Em sites are newly appeared (ia) then add those
                NewDiffDom=ActiveDiffDom(ia,:); %restrict the domains to the new ones only
                [sampts,exitstates]=CalcExitforDomains(NewDiffDom); %calculate the exit properties
                sIDs=ActiveDiffDom_I(ia,1); %indices of new Em sites
                for ii=1:numel(ia)
                    if exitstates(ii)>2
                        nIDs(ii)=0;
                    else
                        nIDs(ii)=ActiveDiffDom_I(ia(ii),exitstates(ii)+1);
                    end
                end
                DiffusionTimeTicker(sIDs,:)=[sampts,sIDs,nIDs]; %add the exit properties to Ticker for new domains
            end
            
        else %if there are no Em sites, then diffusion is not possible. Make sure the lists get cleared out.
            DiffusionTimeTicker=zeros(NumCpGSites,3); %times, sInd,nInd (index of neighbor whose state changes)
            DiffusionTimeTicker(:,1)=Inf;
        end
        
        if ib %remove any sites that are no longer Em
            Prev_Em_inds=Prev_Em_inds_all(ib);
            DiffusionTimeTicker(Prev_Em_inds,2:3)=zeros(numel(ib),2);
            DiffusionTimeTicker(Prev_Em_inds,1)=Inf;
        end
        
        %compare current domains to those stored previously from last
        %iteration.
        [Remaini,ia,ib]=intersect(Curr_Em_inds_all,Prev_Em_inds_all); %get the indices of rows that are common to current and previous step
        
        [rowi,coli]=find(DiffusionDomains(Remaini,2:3)-NewDiffusionDomains(Remaini,2:3));
        if rowi
            getrows=unique(rowi);
            sIDs=Remaini(getrows); %get the rows where something has changed (should be the L or R neighbor position)--should be only 1 row
            [sampts,exitstates]=CalcExitforDomains(NewDiffusionDomains(sIDs,:)); %recalculate the exit properties again for the updated domains
            
            for ii=1:numel(sIDs)
                if exitstates(ii)>2
                    nIDs(ii)=0;
                else
                    nIDs(ii)=ActiveDiffDom_I(ia(getrows(ii)),exitstates(ii)+1);
                end
            end
            DiffusionTimeTicker(sIDs,:)=[sampts(:),sIDs(:),nIDs(:)]; %add the exit properties to Ticker for new domains
        end
        DiffusionDomains=NewDiffusionDomains;
        
    end

    function [sampts,exitstates]=CalcExitforDomains(Domains)
        %takes as input an N by 3 array. Col 1: initial position x0 (initially a
        %site ID of a CpG), Col2: location of nearest neighbor to L,
        %Col 3: location of nearest neighbor to R.
        %outputs the time and exit state by diffusion
        
        DistL=abs(Domains(:,1)-Domains(:,2));
        DistR=abs(Domains(:,1)-Domains(:,3));
        %set up domain info [L,de] (L=size of domain, de=distance from left edge)
        Domains2=[DistL+DistR,DistL]; %express the domains in 2 numbers: [L,de]
        N_Dom=size(Domains,1); %number of domains to keep track of with avaiable hps
        
        %flip them if the current Em is closer to right edge (for looking up FPT
        %function)
        Domains_adj=[Domains2(:,1),min(Domains2(:,2),(Domains2(:,1)-Domains2(:,2)))];
        %Record whether L/R neighbor is closer: if CloseL==0 (right is closer), will need
        %to flip the ordering again later
        CloseL=Domains2(:,2)<=Domains2(:,1)-Domains2(:,2);
        %Now Sample from the CDFs for to get exit time. Also sample
        %from Exit_Probs to get exit location (does it hop L, R, or
        %unbind back to solution with koff)
        Rands=rand(N_Dom,2);
        NoHops=(DistL>MaxDist) & (DistR>MaxDist);
        sampts=zeros(numel(N_Dom),1);
        exitstates=zeros(numel(N_Dom),1);
        for ii=1:N_Dom
            if NoHops(ii)
                %if there are no sites that are reachable by hopping, set the
                %FPT as exponential with rate koff (weak diffusion limit)
                sampt=1/koff*log(1/rand); %sample from exponential (just like KMC)
                exitstate=3;
            else
                
                inds=Domains_adj(ii,:);
                %this is an approximation: if the domain size is too large>MaxL=MaxDist*2,
                %but we have kept it because the enzyme is <MaxDist from one edge, just assume the whole domain is MaxL
                inds(1)=min(inds(1),(MaxDist*2)-1);
                
                CDF=CDFs{inds(1),inds(2)};
                ts=FPT_ts{inds(1),inds(2)}/koff; %divide by koff to recover time in correct units
                %sampling the time
                get_ts=ts(CDF<=Rands(ii,1));
                if isempty(get_ts)
                    sampt=eps;
                    t_ind=1;
                else
                    sampt=get_ts(end); %array of FPTs
                    t_ind=numel(get_ts);
                end
                
                %get the uncorrected Exit Probs
                UnExitP=Exit_Probs{inds(1),inds(2)}(:,t_ind);
                
                %correct the ExitPs: flip the Exit Probs if R neighbor is closer
                ExitP=[UnExitP(2-CloseL(ii)),UnExitP(1+CloseL(ii)),UnExitP(3)];
                %get corresponding exit states based on Exit Probs
                findexitstate=cumsum(ExitP(:))>Rands(ii,2);
                %exitstate =1 (hopped left) =2 (hopped R) =3 (unbind)
                exitstate=sum(1-findexitstate)+1;
            end
            sampts(ii)=sampt;%%min([sampt,Inf]); %fix --if sampt is NaN because of issue with the domain, this makes it Inf %store the sampled FPT
            exitstates(ii)=exitstate;
        end
    end
    function [ActiveDomains,ActiveDomains_Indices] = BuildDiffusionDomains(CurrState,sites)
        %takes as input the current state and site IDS.
        
        %Returns the list of "active" diffusion domains, N rows (# of Em
        %positions), 3 columns: position, L neighbor position, R
        %neighbor position. Domains stores the three positions,
        %Domains_Indices stores the corresponding 3 indices
        Curr_Em_inds=find(CurrState(:,4));
        if Curr_Em_inds
            Curr_Em_sites=sites(Curr_Em_inds);
            N_Em=numel(Curr_Em_sites); %number of sites in Em state
            Curr_h_inds=find(CurrState(:,2));
            Curr_h_sites=sites(Curr_h_inds);
            repC=repmat(Curr_h_sites(:)',N_Em,1);
            ga=repC<Curr_Em_sites(:);
            NLFound=sum(ga,2)>0; %array that tracks whether a L neighor was found. 1=found
            %Neighbor_L is the location of nearest neighbor to left in h
            %state & Neighbor_R is (same) to right
            [Neighbor_L,track]=max(repC.*ga,[],2);
            Neighbor_Li=Curr_h_inds(track);
            ga=repC>Curr_Em_sites(:);
            NRFound=sum(ga,2)>0; %array that tracks whether a L neighor was found. 1=found
            za=repC.*ga;
            MaxSiteID=max(sites); %used later for vectorization
            za(za==0)=MaxSiteID+MaxDist+1; %trick used for vectorization
            [Neighbor_R,track]=min(za,[],2);
            Neighbor_Ri=Curr_h_inds(track);
            %Rarely it can happen that no neighbor is found. In this case
            %we create a "fake" neighbor which is too far away for a hop.
            %The enzyme will always unbind in this case.
            Neighbor_L(NLFound==0)=Curr_Em_sites(NLFound==0)-(MaxDist+1);
            Neighbor_R(NRFound==0)=MaxSiteID+MaxDist+1;
            %We then set the corresponding indices to 0
            Neighbor_Li(NLFound==0)=0;
            Neighbor_Ri(NRFound==0)=0;        
            ActiveDomains=[Curr_Em_sites(:),Neighbor_L(:),Neighbor_R(:)];
            ActiveDomains_Indices=[Curr_Em_inds(:),Neighbor_Li(:),Neighbor_Ri(:)];
        else
        end 
    end
end








