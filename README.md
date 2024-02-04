Code for implementing GameCom and reproducing results from associated paper. All scripts/notebooks that 
require solving an optimization problem have a parameter "cp_solver" at the beginning to specify which solver 
to use with cvxpy. We use Gurobi, but other options are available. See the "Choosing a solver" 
section here for more details: https://www.cvxpy.org/tutorial/advanced/index.html.

Organization of the project as follows:

- EColiModel/
  
  Models, code, and output files for reproducing results in Section 4.1 for the community of four auxotrophic *E. coli.* strains.
  - GameCom/

    Has final data for creating Figure 2 and Table 1 in the main text, and Figures 1 -- 4 in the *SI*, as well as the following Jupyter notebooks for creating this data:
    - ecoli_model_compute_ss_gne.ipynb

      Computes steady state GNE.
    - ecoli_model_stability.ipynb

      Reads in steady state GNE produced by previous notebook, and determines stability to biomass changes and stability to invasion by fifth *E. coli.* strain.
    - ecoli_model_bm_stability_analysis.ipynb

      Analyze metabolism of steady state GNE that are stable and unstable to biomass perturbations. Create Figure 2 (a) in main text and Figure 1 (a) in *SI*, as well as Table 1 in main text.
      Also do FVA analysis and create Figures 2 --4 in the *SI*.
    - ecoli_model_invasion_stability_analysis.ipynb

      Analyze metabolism of steady state GNE that are invaded and uninvaded by fifth *E. coli.* strain. Create Figure 2 (b) in main text and Figure 1 (b) in *SI*.
  - ModelFiles/

    .mat files specifying community metabolism. Taken from SteadyCom code.
    - FourSpecies/
    - FiveSpecies/
  - NECom/

    Implementation of NECom for comparison to GameCom.
    - NECom_Ecoli.ipynb
  - SteadyCom/

    .mat files for using SteadyCom code available from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005539#sec014 to run SteadyCom for comparison and to make .mat files to feed into GameCom.
    - makeEcoliModel.m

      Create .mat files specifying the community model.
    - steadyComFourSpeciesEcoli.m

      Apply SteadyCom to four *E. coli.* strain community.

- GutModel/

  Models, code, and output files for reproducing results in Section 4.2 in the main text and Section S.4.2 in the *SI* for the Gut microbiome model.
  - GameCom/

    Has final data for creating Figure 3 in the main text, as well as the following Jupyter notebooks for creating this data:
    - main.py

      Concise/general implementation of GameCom model (compared to implementation in the EColiModel folder that's specific/hardcoded for the *E. coli.* model.)
    - nine_species_compute_ss_gne.ipynb

      Application of GameCom to gut microbiome model. Computes steady state GNE and their stability to biomass perturbations.
  - ModelFiles/

    .mat files specifying gut microbiome model.
  - SteadyCom/
    - makeNineSpeciesForPython.m

      Create .mat files specifying the community model.
    - simNineSpecies

      Apply SteadyCom to the gut microbiome community.
    - SteadyComSubroutines.m

      Utility functions for the other two files.

