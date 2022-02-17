## Locally-correlated kinetics of post-replication DNA methylation reveals processivity and region-specificity in DNA methylation maintenance
BioRxiv: https://www.biorxiv.org/content/10.1101/2021.09.28.462223v2
### Codes

- `MLERatesInference.m` is used for inferring site-specific remethylation rates via maximum likelihood estimation in this paper.

- The `SimulationCodes/Distributive` and `SimulationCodes/Processive` folder include the stochastic models and input data for Distributive model and Processive model(based on 1D diffusion) used in this paper.

- The `SimulationCodes/kResampling` includes the code for calculating the resampling kinetic rates used as \phi in Equation 7 in the paper, and the Total correlation in Equation 7 comes from the correlation in `correlation/rate_correlations` folder.

### Data

- `correlations` : The rate correlation and wgbs correlation used in the paper (Figure 1).

- Inferred kinetic rates of DNA methylation maintenance: The kinetic rates used in the paper are stored at Google Drive because it is too large to store in Github. https://drive.google.com/drive/folders/1xJ8E6BY1ar4jIcYk62x4wjvA56LjLWyE?usp=sharing

  The files of inferred kinetic rates are stored in .bed format and can be open and edit by any text editor. Columns in each file are Tab-separated and has the following fields:
  column number: content, meaning
	1st: chromosome, the name of the chromosome
	2nd: start, the starting position of the feature in the chromosome 
	3rd: end, the ending position of the feature in the chromosome 
	4th: WGBS, the methylation fraction of the CpG in WGBS data from [GEO:GSM1112841](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1112841)
	5th: k, the estimated kinetic rate of remethylation post-replication using maximum likelihood estimation (MLE)
	6th: f, the estimated steady state methylation fraction (steady-state fraction of cells in population mehylated at that site)
	7th: CpGd, the local CpG density surrounding a given site, which is calculated by the number of neighboring CpGs within a 500bp window
	8th: chromatin_state, the H1-hESC chromatin state at the CpG site using HMM from [ENCODE/Broad](https://www.genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=regulation&hgta_track=wgEncodeBroadHmm&hgta_table=wgEncodeBroadHmmH1hescHMM&hgta_doSchema=describe+table+schema)
	9th: DNase_level, the DNase level which is used to measure the chromatin accessibility at the CpG site retrieved from [ENCODE/OpenChrom(Duke University)](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromDnase/)
	10th: NOS, Neucleo Occupancy Score at the CpG site obtained from [GSM1194220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1194220)
