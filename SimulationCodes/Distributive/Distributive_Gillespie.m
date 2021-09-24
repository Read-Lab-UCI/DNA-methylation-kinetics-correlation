function [SimDat] = Distributive_Gillespie(sites,fracmeth,TimePts,SimName)

%Simulates post-replication methylation according to Distributive Model
%using Gillespie algorithm

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

ParamFileName=['Params_' SimName];
ParamArray=[Ratio,k1f,k1r,kcat,S0];
save(ParamFileName,'ParamArray')

E=max(1,round(Ratio*S0)); %total number of Enzymes in simulation

NumRxns=3; %number of reactions
NumSpec=4; %number of species, indexed as: [u,h,Eh,m]

%calc predicted mean half time
Km_stoch=(k1r+kcat)/k1f;
t50=(Km_stoch+S0)/(2*kcat*E);

SimDat=zeros(NumCpGSites,2); %initialize output file

Stoich=zeros(NumRxns,NumSpec); %initialize stoichiometry array

%Rxn 1; h + E -> Eh, k1f
Stoich(1,2)=-1;
Stoich(1,3)=1;
%Rxn 2; Eh -> h + E, k1r
Stoich(2,2)=1;
Stoich(2,3)=-1;
%Rxn 3; Eh -> E + m kcat
Stoich(3,3)=-1;
Stoich(3,4)=1;

%Initialize the starting state with some hemimethylated sites
SitesArray0=zeros(NumCpGSites,NumSpec);
%sample from the input methylation landscape to get initial states (u or h)
InitRandArray=rand(NumCpGSites,1);
AssignedFrac=InitRandArray<=fracmeth;
SitesArray0(AssignedFrac,2)=1; %identify the sites that are h in the starting state
SitesArray0(:,1)=abs(SitesArray0(:,2)-1); %the sites that are not h in the starting state are u
CurrState=SitesArray0;

%call the propensity function. a is array with propensity of each reaction
%(columns) at each site (rows)
[a,a_0]=GetPropensity(CurrState);


Time=0;
count=0;

while Time<TimePts(end)
    %Gillespie algorithm: 2 random numbers, first to get time to next
    %reaction, 2nd to select which reaction
    Rands=rand(1,2); %pick 2 random numbers
    tau=1/a_0*log(1/Rands(1)); %compute time to next reaction.
    count=count+1; %count the number of steps
    if isinf(tau) %if the time-to-next is infinite, increment anyway
        Time=Time+.5;
    else
        Time=Time+tau; %increment time
        findnextrxn=cumsum(a(:))>Rands(2)*a_0;
        NextRxn=sum(1-findnextrxn)+1; %find the position in the cumsum array of the reaction that takes place
        
        [I,J]=ind2sub(size(a),NextRxn); %find the position in a of the reaction that takes place (J) and the site in which takes places (I)
        CurrState(I,:)=CurrState(I,:)+Stoich(J,:); %update the state of the system according to reaction event
    end
    
    %Reset the propensity array according to new state of system
    [a,a_0]=GetPropensity(CurrState); %call the function that computes the propensity for each site
end
Time
SimDat(:,1)=SimDat(:,1)+CurrState(:,4); %methylated reads (m)
SimDat(:,2)=SimDat(:,2)+CurrState(:,1)+CurrState(:,2)+CurrState(:,3); %unmethylated reads (u,h, or Eh)

    function [a,a_0]=GetPropensity(CurrState)
        %Calc the per-site propensities for each reaction. Array is
        %NumCpGs x NumRxns.
        
         %Inputs: CurrState array (NumCpGSitesxNumSpec)
        %Outputs: 
            % a: (NumCpGSites x NumSpec) the propensities of each reaction
                    % at each site
            % a0: (double) sum over all propensities
        
        Enum=E-sum(CurrState(:,3)); %number of free enzyme copies
        %Calc the per-site propensities for each reaction.
        a=zeros(NumCpGSites,NumRxns);
        a(:,1)=k1f*Enum*CurrState(:,2); %E + h -> [Eh]
        a(:,2)=k1r*CurrState(:,3); %  [Eh] -> E + h
        a(:,3)=kcat*CurrState(:,3); %Eh -> E + m
        a_0=sum(a(:));
    end
end



