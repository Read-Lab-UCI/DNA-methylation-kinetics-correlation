## Locally-correlated kinetics of post-replication DNA methylation reveals processivity and region-specificity in DNA methylation maintenance

### Codes

- `MLERatesInference.m` is used for inferring site-specific remethylation rates via maximum likelihood estimation in this paper.

- The `SimulationCodes/Distributive` and `SimulationCodes/Processive` folder include the stochastic models and input data for Distributive model and Processive model(based on 1D diffusion) used in this paper.

- The `SimulationCodes/kResampling` includes the code for calculating the resampling kinetic rates used as \phi in Equation 7 in the paper, and the Total correlation in Equation 7 comes from the correlation in `correlation/rate_correlations` folder.

### Data

- Inferred kinetic rates of DNA methylation maintenance: The kinetic rates used in the paper are stored at Google Drive because it is too large to store in Github. https://drive.google.com/drive/folders/1xJ8E6BY1ar4jIcYk62x4wjvA56LjLWyE?usp=sharing

- `correlations` : The rate correlation and wgbs correlation used in the paper (Figure 1).
