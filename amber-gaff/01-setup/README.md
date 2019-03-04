The scripts are available from 
https://github.com/michellab/BioSimSpace/tree/devel/nodes/playground

## Step 1 BSS : parameterise a protein. The pdb has been split into 3 pdbs, one for each molecule in the protein.     
parameterise.py --input protein.pdb/protein-model-01.pdb --forcefield ff14SB --output _protein/protein-model-01
parameterise.py --input protein.pdb/protein-model-02.pdb --forcefield ff14SB --output _protein/protein-model-02
parameterise.py --input protein.pdb/protein-model-03.pdb --forcefield ff14SB --output _protein/protein-model-03
### Now combine
combine.py --system1 _protein/protein-model-01.prm7 _protein/protein-model-01.rst7 --system2 _protein/protein-model-02.prm7 _protein/protein-model-02.rst7 --output _protein/prot-1plus2
combine.py --system1 _protein/prot-1plus2.prm7 _protein/prot-1plus2.rst7 --system2 _protein/protein-model-03 --output _protein/protein-gas
## Step 2 BSS: Parameterise a list of ligands 
parameterise.py --input ligands.pdb/12.pdb --forcefield gaff2 --output _ligands/12
parameterise.py --input ligands.pdb/14d.pdb --forcefield gaff2 --output _ligands/14d
## Step 3 BSS: Solvate each ligand  # solvate.py -
solvate.py --input _ligands/12.prm7 _ligands/12.rst7 --output _solvated-ligands/12_water --water tip3p --extent 20
solvate.py --input _ligands/14d.prm7 _ligands/14d.rst7 --output _solvated-ligands/14d_water --water tip3p --extent 20
## Step 4 BSS: Assemble each protein-ligand complex and solvate 
combine.py --system1 _ligands/12.prm7 _ligands/12.rst7 --system2 _protein/proteiin-gas.prm7 _protein/protein-gas.rst7 --output prot-lig
solvate.py --input prot-lig.prm7 prot-lig.rst7 --output _solvated-complexes/12_prot_water --water tip3p --extent 20
rm prot-lig.prm7 prot-lig.rst7
combine.py --system1 _ligands/14d.prm7 _ligands/14d.rst7 --system2 _protein/protein-gas.prm7 _protein/protein-gas.rst7 --output prot-lig
solvate.py --input prot-lig.prm7 prot-lig.rst7 --output _solvated-complexes/14d_prot_water --water tip3p --extent 20
## Step 5 BSS Run an MD equilibration for each solvated system 
### This script to be written. Julien will check if existing scripts can do this. Should do minimisation, 100 ps NVT, 900 ps NPT using SOMD ideally
### equilibrate.py --input _solvated-ligands/12_water --protocol default
### equilibrate.py --input _solvated-ligands/14d_water --protocol default
### equilibrate.py --input _solvated-complexes/12_prot_water --protocol default
### equilibrate.py --input _solvated-complexes/14d_prot_water --protocol default

## Step 6 BSS: run prepareFEP on every pair of ligands
prepareFEP.py --input1 _solvated-ligands/12_water.prm7 _solvated-ligands/12_water_eq.rst7 --input2 _solvated-ligands/14d_water.prm7 _solvated-ligands/14d_water_eq.rst7 --output _solvated-perturbations-ligands/12_to_14d_free
prepareFEP.py --input1 _solvated-ligands/12_prot_water.prm7 _solvated-ligands/12_prot_water_eq.rst7 --input2 _solvated-ligands/14d_prot_water.prm7 _solvated-ligands/14d_prot_water_eq.rst7 --output _solvated-perturbations-complexes/12_to_14d_bound

This will complete the setup of the input files for a SOMD relative FEP. 
The next steps to implement should be in a separate folder that can be packed and copied around. 
To implement them, create a new subfolder in amber-gaff and create a protocol and run folder. e.g.

01-setup
02-protocol
03-run

The contents of these folders should be similar to what was done here

https://github.com/michellab/somdfreenrgbatch

You will need to edit the script in run/setup.py to collect the relevant input files from the 01-setup folder. 

