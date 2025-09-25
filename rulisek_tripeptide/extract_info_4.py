import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals
from MDAnalysis.analysis import distances
import os
if not os.path.exists('extracted_info'): os.mkdir('extracted_info')
import multiprocessing

data_dir = "P-CONF_1.6M"
sequences = os.listdir(data_dir)
directories = [f"{data_dir}/{sequence}" for sequence in sequences]

# load list of unusually high-energy conformers that the authors ignored but zipped up with everything else (see https://osf.io/7nwu4/wiki/home/)
they_said_to_ignore = []
with open('they_said_to_ignore_these.txt','r') as f:
    for line in f:
        they_said_to_ignore.append(line.strip())
# load list of PDB files without energy information in their headers, which will also be ignored
structures_without_energies = []
with open('structures_without_energies.txt','r') as f:
    for line in f:
        structures_without_energies.append(line.strip())
to_ignore = they_said_to_ignore + structures_without_energies 
to_ignore.append('DIR167.pdb') # for some reason, this one structure is missing the C-terminal cap, causing errors in the analysis
# putting main logic in function so it can be parallelized
def process_all_confs_of_seq(directory):
    all_confs_info = []
    for filename in os.listdir(directory):
        # ignore unusable files
        if filename in to_ignore:
            continue
        if filename[-4:] != ".pdb": continue
        # list to hold energy and angles
        conf_info = []
        # get energies of conformers of the given sequence
        with open(f"{directory}/{filename}",'r') as f, open("failed_extraction.txt", 'a') as bad_log:
            for counter, line in enumerate(f):
                if counter < 4:
                    try:
                        conf_info.append(float(line.split("=")[1].strip()))
                        # index 0: water absolute energy
                        # index 1: water energy relative to global min for the current sequence
                        # index 2: octanol absolute energy
                        # index 3: octanol energy relative to global min for the current sequence
                    except:
                        bad_log.write(f"could not identify energies: {directory}/{filename}\n")
                        print(f"could not identify energies: {directory}/{filename}\n")
                        raise
                else:
                    break
        # get phi and psi dihedral angles of middle residue in structure
        u = mda.Universe(f"{directory}/{filename}",dt=1) # timestep value for a single frame is meaningless, just providing it to silence the warning
        rama_output = dihedrals.Ramachandran(u.atoms, c_name='C', n_name='N', ca_name='CA',check_protein=True).run()
        assert rama_output.results.angles.shape == (1,3,2) # 1 frame x 3 full AA residues x 2 rama angles per residue
        try:
            phi = rama_output.results.angles[0,1,0] # axis 0 is frame, axis 1 is residue position in selection, axis 2 is angle type
        except IndexError:
            print(f"BAD phi: {directory}, {filename}")
            raise
        try:
            psi = rama_output.results.angles[0,1,1]
        except IndexError:
            print(f"BAD psi: {directory}, {filename}")
            raise
        conf_info.append(phi)
        conf_info.append(psi)
        # get distance between CA atoms of adjacent residues
        CA_CA_dist = distances.self_distance_array(u.select_atoms("(resnum 2 or resnum 4) and name CA")) # resnum follows the topology-which is 1-indexed, so the N-terminal cap is 1 and the first real residue is 2
        conf_info.append(CA_CA_dist[0]) # CA_CA_dist is a list with 1 element
        # log combined list of energies and angles before proceeding with the next iteration of the loop
        all_confs_info.append(conf_info)
    np.save(f"extracted_info/{directory.split('/')[-1][:3]}.npy", np.array(all_confs_info)) # isolate sequence name and save it

with multiprocessing.Pool() as pool: 
    pool.map(process_all_confs_of_seq, directories)
#for directory in directories:
#    process_all_confs_of_seq(directory)

