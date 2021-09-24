function [CDFs,FPT_ts,Exit_Probs,MaxDist] = MakeFPTLookup_analytical(D)
%this function calculates the dynamics of a protein diffusing along DNA. It diffuses
%within a 1D domain and can exit the domain either by exiting left boundary,
%exiting right boundary, or unbinding

%input args: diffusion coefficient D in units of bp^s/tau, where tau = 1/koff
%this is equiv to the ratio of D to koff (assuming consistent time units)

%output args: cell arrays containing lookup information for:
%CDFs - cumulative distribution functions for the exit times
%Exit_Probs - time-dependent relative exit probabilities to exit L,R, or to solution
%FPT_ts - time points at which CDFs and Exit_Probs are calculated. Time
%points are in time units of tau
%MaxDist - the computed max distance we need to consider for possible
%"superhop" by diffusion

%the cell arrays are computed for all possible domain sizes (up to a reasonable limit) and initial
%positions within those domains. Cell arrays are MaxL*MaxDist where rows are domain size and
%columns are initial distance from left edge

%first decide how big the arrays need to be -> decide on MaxL. How large of
%domain do we need to consider before exit by diffusion is negligible (it
%will almost always exit by unbinding). Use an analytical calculation
%comparing Mean Exit Time in pure diffusion limit vs pure unbinding limit. 
Ratio=10; %Assume the ratio of MFPT_diff/MFPT_unbind cannot exceed this number
%(i.e., if MFPT_diff is too much longer than MFPT_unbind, then diffusion to
%exit can be neglected)
MaxL=round(sqrt(8*D*Ratio)); %based on an analytical formula for MFPT_diff
MaxDist=ceil(MaxL/2); %the furthest the protein can be from either end (halfway)

%initialize the cell arrays to store all the computed CDFs, etc.
CDFs=cell(MaxL,MaxDist);
FPT_ts=cell(MaxL,MaxDist);
Exit_Probs=cell(MaxL,MaxDist);
NTerms=300; %number of terms to keep in the Fourier series solution.
ns=1:NTerms;

for L=2:MaxL %loop over the domain size defined as number of non-target states between the boundaries
    disp(['Computing FPTs for Domain Size' num2str(L)])
    MaxdeL=ceil(L/2); %max distance to either edge is halfway between
    
    for x0=1:MaxdeL %loop over the initial position within the domain
        
        %Functions below are based on analytical Fourier series solution
        %term by term
        Pterm = @(x,t,n) 2/L.*sin(n*pi*x0/L).*sin(n*pi*x/L).*exp((-D*n.^2*pi^2./L^2-1).*t);
        phi = @(n) 2./(n*pi).*sin(n*pi*x0/L).*(1-(-1).^n);
        theta = @(n) D*n.^2*pi^2./L^2 + 1;
        %survival probability S(t)=int_0_L P(x,t)dx
        Sterm = @(n,t) phi(n).*exp(-theta(n).*t);
     
        MFPTarray=phi(ns)./theta(ns);
        MFPT_ana=sum(MFPTarray); %Mean First Passage Time
        
        %set up a grid of non-uniform timepoints at which to calc FPTD (assuming
        %need more timepoints early) based on MFPT estimate
        t_end=10*MFPT_ana; %number of times the estimated MFPT to go to (max time)
        t_mid=MFPT_ana;
        ntimes=3000; %number of timepoints to compute FPTD at
        t_early=linspace(0,t_mid,ntimes/2);
        t_axis=[t_early(1:end-1),linspace(t_early(end),t_end,ntimes/2)];
        
        %compute arrays--later need to sum over terms in FOurier series
        ns=reshape(ns,1,numel(ns));
        %calculate time-evolution of Probability. Need these to compute
        %fluxes to L and R absorbing states
        Parray=Pterm(1,t_axis(:),ns);
        P_ana_L=sum(Parray,2); %the probability to be in left-most state (x=1)
        Parray=Pterm(L-1,t_axis(:),ns);
        P_ana_R=sum(Parray,2); %the probability to be in right-most state (x=L-1)
        
        %Compute the survival probability. Need this to sample First%Passage Time
        Sarray=Sterm(ns,t_axis(:));
        S_ana=sum(Sarray,2);
        
        dPsol=1*S_ana; %equals dP (prob flux) to end state (solution, due to unbinding) %here assuming koff=1, %correction of timeunits later
        
        dPL=D*P_ana_L; %relative flux to left boundary
        dPR=D*P_ana_R; %relative flux to right boundary
        
        TotFlux=dPL+dPR+dPsol;
        RelFluxes_ana=[dPL./TotFlux,dPR./TotFlux,dPsol./TotFlux];
        CDF_ana=1-S_ana;
        CDFs{L,x0}=CDF_ana;
        FPT_ts{L,x0}=t_axis;
        Exit_Probs{L,x0}=RelFluxes_ana';
    end
end
end
