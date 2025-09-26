This repository contains scripts to analyze two databases of tripeptide conformers:

Dicks, L.; Wales, D. J. Exploiting Sequence-Dependent Rotamer Information in Global Optimization of Proteins. _J. Phys. Chem. B._ **2022**, _126_, 8381-8390.

Culka, M.; Kalvoda, T.; Gutten, O.; Rulisek, L. Mapping Conformational Space of All 8000 Tripeptides by Quantum Chemical Methods: What Strain Is Affordable within Folded Protein Chains? _J. Phys. Chem. B._ **2021**, _125_, 58-69.

Instructions to download and extract information from these databases can be found in the READMEs in the wales_tripeptide and rulisek_tripeptide folders.

No installation is necessary for this repository. Only NumPy and MDAnalysis python packages are needed to perform the analysis.

The two databases can be compared using the scatterplots_by_seq.py in this directory. Running the script as is plots the phi and psi dihedral angles of the middle residue of each sampled peptide from Dicks and Wales and the dihedral angles for each peptide from the Culka et al. MD-generated conformer set, with the conformational free energy evaluated in water.

The rulisek_tripeptide directory also has a script to build tabulated dihedral potentials for use with OpenAWSEM. I chose to use the MD-water conformational set from the Culka et al. paper because I want this term to be as physically realistic as possible, without building in any bias from the PDB, and I wanted to have as many samples as possible to construct my potential.
