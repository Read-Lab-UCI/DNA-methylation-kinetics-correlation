function [] = MLERatesInference(AllDatDir, RatesDir, StartChr, EndChr)

%function to load in Repli-BS dataset and fit site-specific remethylation rates
%(k) and steady-state methylation values (f) using maximum likelihood
%estimation, according to the model Prob(methylated read at time t at site)=
%f*(1-exp(-k*t)).
%Protocol refers to the choice of fitting parameters, including: read-depth
%requirements (only sites with sufficient depth are retained in fitting)
%and the method for dealing with edge-cases, where parameters are not
%strictly identifiable.
clc;
close all;
if ~exist(RatesDir, 'dir')
    mkdir(RatesDir)
end

%set the per-site read-depth requirement for estimating remethylation kinetics
EnoughReads0=5; %number of reads required at t=0
EnoughReadsLater=5; %number of reads required for later timepoints
EnoughMethReads=5; %total number of methylated reads required

for chromosome = StartChr : EndChr
    inputReadDataPath = strcat(AllDatDir,'/AllDat_chr',int2str(chromosome),'.mat'); % file path of read data
    inferredRateSavingPath = strcat(RatesDir, '/Rates_chr',int2str(chromosome),'.mat');
    [fittedSites, inferredRates, inferredMethyFrac, AllDat] = kineticRatesMLE(inputReadDataPath,inferredRateSavingPath, EnoughReads0, EnoughReadsLater, EnoughMethReads);
end
end

function [fittedSites, inferredRates, inferredMethyFrac, AllDat] = kineticRatesMLE(inputReadDataPath,inferredRateSavingPath, EnoughReads0, EnoughReadsLater, EnoughMethReads)
% kineticRatesMLE - Maximum likelihood estimation (MLE) for site-specific
% kinetic rates of DNA methylation maintenance.
%
% Inputs:
%    inputReadDataPath - Path of input read data in .mat format
%    inferredRateSavingPath - Path for saving the inferred kinetic rates
%
% Outputs:
%    fittedSites - The location of the fitted CpG sites on genome: N integers, N
%    is the number of CpGs
%
%    inferredRates - The inferred kinetic rates with upper and lower
%                   confidence intervals(CIs), N by 4 matrix
%
%    inferredMethyFrac - The inferred steady state methylation fraction with upper and lower
%                   confidence intervals(CIs), N by 3 matrix
%
% Example:
%    kineticRatesMLE('../data/ReadData/AllDat_Chr1WT.mat',
%                    '../data/InferedRates/kineticRateChr1.mat')

% load read data
disp(sprintf('Inferring for %s now ...', inputReadDataPath));
load(inputReadDataPath, 'AllDat', 'sites');

% AllDat is an array of size (NSites, NTimepoints, 2).
% AllDat(i, j, 1) is the number of methylated reads at site i at timepoint j.
% AllDat(i, j, 2) is the number of unmethylated reads at site i at timepoint j.
% Totalreads(i, j) = AllDat(i, j, 1) + AllDat(i, j, 2).
% sites is an CpG location array with size (NSites, 1)

% timePoints : the experimental timepoints, and they are shifted by tShift,
% which accounts for the fact that at, e.g., time=0 hrs, the captured reads
% started replicating within the previous hour. This is done also because
% the model assumes zero probability of methylation occuring by time=0.
tShift = 0.5;
timePoints = [0, 1, 4, 16] + tShift;

%determine which sites in the data set have sufficient read depth
Reads=sum(AllDat(:,:,1:2),3);
NumReads0=Reads(:,1);
NumReadsLater=sum(Reads(:,2:end),2);
NumMethReads=sum(AllDat(:,:,1),2);
Conditions=[NumReads0>=EnoughReads0,NumReadsLater>=EnoughReadsLater,NumMethReads>=EnoughMethReads];
KeepSites=find(prod(Conditions,2));
numSites=numel(KeepSites);

methyFracGrid = 0 : 0.01 : 1;    %grid of f-values at which LogLikelihood will be computed
rateGrid = 10.^(-1.5 : .01 : 1.5);   %grid of rate-values (denoted as k hereafter)

%initialize various arrays
inferredRates = zeros(numSites, 3); %The inferred kinetic rates with upper, lower confidence intervals(CIs) and Edge Flag, N by 4 matrix

%The Edge Flag(4th field) is defined as follows:
%-1: sites where k is undefined, set k=0 (generally because no methylated reads-->f=0)
%0: normal sites - k and full CI95 are identifiable
%1: sites on edge: only lower bound of k is identifiable.
%Reported k is an estimate of lower bound, corresponding to lower CI75
%2: sites where maxLL is not on boundary, but CI95 is. So k is
%identifiable, but lower CI95 bound is not. Reported upper CI95 is a lower-bound estimate

inferredMethyFrac = zeros(numSites, 3); %The inferred steady state methylation fraction with upper and lower confidence intervals(CIs), N by 3 matrix
fittedSites = sites(KeepSites)';

tic

%loop over the sites for rate and steady state methylation level MLE
parfor ii = 1 : numSites  %parfor
    %get the read data for site ii
    methylatedReadTimeCourseForSiteii = AllDat(KeepSites(ii), : , 1);
    unmethylatedReadTimeCourseForSiteii = AllDat(KeepSites(ii), : , 2);
    [inferredRates(ii, : ), inferredMethyFrac(ii, : ), LogLikelihood] = siteMLE(rateGrid, methyFracGrid, timePoints, methylatedReadTimeCourseForSiteii, unmethylatedReadTimeCourseForSiteii);
end
save(inferredRateSavingPath, 'fittedSites', 'inferredRates', 'inferredMethyFrac')
toc
end

function [inferredRateAndCIs, inferredMethyFracAndCIs, LogLikelihood] = siteMLE(rateGrid, methyFracGrid, timePoints, methylatedReadTimeCourseForSiteii, unmethylatedReadTimeCourseForSiteii)
% calcLikelihoodSurface - calculate the likelihood surface
%
% Inputs:
%    rateGrid - the rate grid to calculate likelihood, 1 x n1 double vector
%    methyFracGrid - the methylation fraction grid to calculate likelihood, 1 x n2 double vector
%    timePoints - the timepoints at which the reads are experimentally measured, 1 x 4 double vector
%    methylatedReadTimeCourse - the time course counts of methylated reads for the given timepoints, 1 x 4 double vector
%    unmethylatedReadTimeCourse - the time course counts of unmethylated reads for the given timepoints, 1 x 4 double vector
% Outputs:
%    inferredRateAndCIs - inferred kinetic rates and its confidence interval
%4th element is a flag:1=no maximum could be identified, so rate
%       value gives a bound. 2 = maximum identified, rate is identifiable
%    inferredMethyFracAndCIs - inferred methylation fraction and its confidence interval
%    LogLikelihood - the likelihood surface
% Example:
%    [inferredRateAndCIs, inferredMethyFracAndCIs] = siteMLE(10.^(-2: .01: 1), 0: 0.01: 1, [0.5, 1.5, 4.5, 16.5], [7, 0, 2, 3], [1, 0, 5, 6])

LogLikelihood = calcLikelihoodSurface(rateGrid, methyFracGrid, timePoints, methylatedReadTimeCourseForSiteii, unmethylatedReadTimeCourseForSiteii);

%CIs will be computed based on Likelihood Ratio Test, using Chi-squared values
PrcVal95=3.841; %95th of chi-squared distribution, 1 free param.
PrcValIn=1.32; %Inner percentile. Corresponds to rth %ile of chi-squared, r is user chosen
PrcValEdge=0.102; %Inner percentile. Corresponds to qth %ile of chi-squared, q is user chosen

%PrcValIn is used to determine identifiability--if the entire rate region
%lies within rth percent CI of the maximum, then the rate is considered
%unidentifiable.

%PrcValEdge is used to set the rate in edge cases--i.e. when the maximum is
%found to be on or close to the edge of the region. In this case, the qth
%percentile is chosen as a bound on the rate.


maxRate = max(rateGrid);
minRate = min(rateGrid);

logMinRate = log10(maxRate);
logMaxRate = log10(maxRate);

%Find the edge confidence intervals
[RawRate,RawFrac,CIRateEdge,CIfracEdge,MaxLL]= getCI(rateGrid, methyFracGrid, LogLikelihood, PrcValEdge);
[RawRate,RawFrac,CIRateIn,CIfracIn,MaxLL]= getCI(rateGrid, methyFracGrid, LogLikelihood, PrcValIn);
[RawRate,RawFrac,CIRate95,CIfrac95,MaxLL]= getCI(rateGrid, methyFracGrid, LogLikelihood, PrcVal95);


%CIs are computed based on Likelihood Ratio Test, using Chi-squared values

%check to see whether the inner confidence interval hits the boundary
CheckCIsL=[minRate,maxRate]-CIRateIn;
keepindlam=find(abs(CheckCIsL)>realmin);

if numel(keepindlam)==0 %this indicates the inner CI overlaps entire region-->k is considered unidentifiable.
    %Set k==0 (do not attempt to fit a rate)
    inferredRateAndCIs=[0,0,0];
    inferredMethyFracAndCIs=[0,0,0];
elseif numel(keepindlam)==1 %this indicates CI value hits one boundary-->indicated upper/lower bound on k
    %when the CI hits the boundary, use the edge CI value as the rate
    %estimate.
    inferredRateAndCIs=[CIRateEdge(keepindlam),CIRate95(1),CIRate95(2)];
    lamind=find(rateGrid==CIRateEdge(keepindlam));
    [aa,bb]=max(LogLikelihood(:,lamind));
    inferredMethyFracAndCIs=[methyFracGrid(bb),CIfrac95(1),CIfrac95(2)];
else
    inferredRateAndCIs=[RawRate,CIRate95(1),CIRate95(2)];
    inferredMethyFracAndCIs=[RawFrac,CIfrac95(1),CIfrac95(2)];
end
end


function [LogLikelihood] = calcLikelihoodSurface(rateGrid, methyFracGrid, timePoints, methylatedReadTimeCourse, unmethylatedReadTimeCourse)
% calcLikelihoodSurface - calculate the likelihood surface
%
% Inputs:
%    rateGrid - the rate grid to calculate likelihood, 1 x n1 double vector
%    methyFracGrid - the methylation fraction grid to calculate likelihood, 1 x n2 double vector
%    timePoints - the timepoints at which the reads are experimentally measured, 1 x 4 double vector
%    methylatedReadTimeCourse - the time course counts of methylated reads for the given timepoints, 1 x 4 double vector
%    unmethylatedReadTimeCourse - the time course counts of unmethylated reads for the given timepoints, 1 x 4 double vector
% Outputs:
%    LogLikelihood - The logLikelihood surface for given rate and methylation fraction grid, n1 x n2 double matrix
%
% Example:
%    calcLikelihoodSurface(10.^(-2: .01: 1), 0: 0.01: 1, [0.5, 1.5, 4.5, 16.5], [7, 0, 2, 3], [1, 0, 5, 6])

pMethyArray = zeros(numel(methyFracGrid), numel(rateGrid), numel(timePoints));
pUnmethyArray = zeros(numel(methyFracGrid), numel(rateGrid), numel(timePoints));

%loop over the timepoints to compute the LogLikelihood surface
for tind = 1 : numel(timePoints)
    timeI = timePoints(tind);
    expoentialTerms = 1 - exp(-rateGrid .* timeI);
    expoentialTerms(expoentialTerms >= 1) = 1 - eps;
    repExpoentialTerms = repmat(expoentialTerms, numel(methyFracGrid), 1);
    repMethyFrac = repmat(methyFracGrid', 1, numel(rateGrid));
    pMethyArray( : , : , tind) = methylatedReadTimeCourse(tind) .* log(repMethyFrac .* repExpoentialTerms);
    pUnmethyArray( : , : , tind) = unmethylatedReadTimeCourse(tind) .* log(1 - repMethyFrac .* repExpoentialTerms);
end
%this is the LogLikelihood surface as a function of the parameters
LogLikelihood = sum(pMethyArray, 3, 'omitnan') + sum(pUnmethyArray, 3, 'omitnan');
end


function [RawRate, RawFrac, CIRate, CIFrac, MaxLL] = getCI(rateGrid, methyFracGrid, logLikelihood, prcVal)
% getCI - get confidence intervals of rate of remethylation and logLikelihood
%
% Inputs:
%    rateGrid - the rate grid to calculate likelihood, 1 x n1 double vector
%    methyFracGrid - the methylation fraction grid to calculate likelihood, 1 x n2 double vector
%    logLikelihood -  the likelihood surface
%    prcVal - the percentile value of chi-squared distribution with 1 free param
% Outputs:
%    RawRate - the raw max likelihood parameters of rate
%    RawFrac - the raw max likelihood parameters of methylation fraction
%    CIRate - confidence intervals of rate of remethylation and logLikelihood
%    CIFrac - confidence intervals of methylation fraction
%    MaxLL - the max value of logLikelihood

%find the parameter values that maximize the logLikelihood
szL=size(logLikelihood);
[MaxLL,i]=max(logLikelihood(:));
[I,J]=ind2sub(szL,i); %I,J are the indices of ML parameters

%these are the "raw" ML parameters, but need to deal with edge cases (below)
RawRate=rateGrid(J);
RawFrac=methyFracGrid(I);

%compute the profile likelihood function.
ProfileRate=max(logLikelihood,[],1);
ProfileFrac=max(logLikelihood,[],2);

%find indices of k values with LL (logLikelihood) values within CI
InnerInds=find(-2*(ProfileRate-logLikelihood(I,J)) <= prcVal);

%Find the confidence interval
RateRegion=rateGrid(InnerInds);
FracRegion=methyFracGrid(-2*(ProfileFrac-logLikelihood(I,J)) <= prcVal);

%write an array: [CI LB, CI UB, LL LB, LL UB]
CIRate=[RateRegion(1) RateRegion(end)];
CIFrac=[FracRegion(1) FracRegion(end)];
end

