## md Folder
This describes the MD simulations code for the paper. This involves running simulations studying self-assembly of coarse-grained representations of NEAT1 with proteins and for transcription dynamics of NEAT1. 

## Organization
- Code has the scripts for running and analyzing MD simulations used in this paper.
- Examples has notebooks to demonstrate how to run simulations of NEAT1 with NONO, FUS, and TDP-43, and simulations for transcription dynamics of NEAT1. 


## Requirements
- List of libraries:
    - hoomd
    - gsd
    - freud
    - matplotlib
    - seaborn
    - numpy
    - pandas
    - pickle
- A complete copy of the local environment file with pinned dependencies is linked here: [paraspeckle environment](paraspeckles.yml)

## Details
1. Under Code
    - AnalyzeData: Scripts to determine average data including MIPs and RDFs for trajectories
    - SimulationCode: Scripts to initalize and run a simulation with NEAT1 and proteins

2. Under Examples:
    - NEAT1andProteins contains a notebook to generate trajectories of NEAT1 with proteins, and analyze the trajectories for the largest cluster in each simulation. The notebook will generate a folder called Data organized by experiment name, parameter type, file type, then date. 
        Example Format of Data/ExperimentName/ParameterType/
        - DataAnalysis: Folder to hold analyzed data from AnalyzingNEAT1andProteins. Organized by type of analysis then date. 
        - InputFiles: Folder for initial states made for each trajectory ran organized by date. 
        - JobInputs: Folder for the parameters used for each trajectory ran organized by date. 
        - Outputs: Trajectories organized by date.

    - TranscriptionDynamics contains a notebook to generate and analyze trajectories with NEAT1 copies added over steps of the simulaiton. 
    
    
