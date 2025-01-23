## Indirect evolutionary effects increase the maintenance of polymorphisms in mutualistic networks

This README file was generated on 23.01.2025 by Leandro Giacobelli Cosmo

GENERAL INFORMATION

1. Title of Dataset: 

Cosmo_et_al_coevo_polymorphism

2. Author information:

Leandro G. Cosmo: Department of Evolutionary Biology and Environmental Studies, University of Zurich - Zurich, Switzerland
 
Paulo R. Guimarães Jr.: Departamento de Ecologia, Instituto de Biociências, Universidade de São Paulo - USP, São Paulo, SP, Brazil.

Jordi Bascompte: Department of Evolutionary Biology and Environmental Studies, University of Zurich - Zurich, Switzerland

Corresponding author: Leandro G. Cosmo, E-Mail: Leandro.giacobellicosmo@uzh.h

DATA & FILE OVERVIEW

1. File List: 

Data files:

main_mean_results.csv

Scripts/Source functions:

main_simulations.jl\
aux_functions.jl

DATA-SPECIFIC INFORMATION:

main_mean_results.csv: full dataset containing the average result of the numerical simulations used for the analyses in the main text.\The variables in the dataset correspond to: 

(1) network: name of the mutualistic network used, as defined in the Web of Life database.\
(2) sp_id: ID of the species within the network.\
(3) type: type of the species (plant or animal).\
(4) position: position of the species within the network (core or periphery).\
(5) degree: Degree (number of direct interactions) of the species.\
(6) peq: Frequency of the morph A at the coevolutionary equilibrium.\
(7) qeq: - Frequency of the morph B at the coevolutionary equilibrium.\
(8) pvar: Morph variance at the coevolutionary equilibrium, defined as peq*qeq.\
(9) y: Parameter of the model that represents the average strength of mutualistic interactions over the two morphs.\
(10) yA: Parameter of the model that controls the strength of mutualistic interactions for morph A.\
(11) yB: Parameter of the model that controls the strength of mutualistic interactions for morph B.\
(12) T_i: Total amount of direct and indirect evolutionary effects received by the species.\
(13) F_i - Indirect evolutionary effects received by the species.\

main_functions.jl and aux_functions.jl: Julia functions used to run the model numerical simulations.\

USAGE INSTRUCTIONS:

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> Cosmo_et_al_coevo_polymorphism

It is authored by Cosmo et al.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths. 

After installing everything, run the script "main_simulations.jl" located at the "scripts" folder.
