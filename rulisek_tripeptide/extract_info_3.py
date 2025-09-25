import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals
from MDAnalysis.analysis import distances
import os
import multiprocessing

data_dir = "P-CONF_1.6M"
sequences = os.listdir(data_dir)
already_done = os.listdir('extracted_info')
directories = [f"{data_dir}/{sequence}" for sequence in sequences if f"{sequence}.npy" not in already_done]
#directories = [f"{data_dir}/{sequence}" for sequence in sequences]
# DIR peptide seems to be causing issues. I don't know why DIR was causing issues. it was some sort of issue with parsing the R residue, but nothing looks wrong with the pdb file. I'll just write a numpy array of all zeros for now as a workaround
directories = [f"{data_dir}/DIR"]

def process_all_confs_of_seq(directory):
    all_confs_info = []
    for filename in os.listdir(directory):
        conf_info = []
        if filename[-4:] != ".pdb": # i might put some of my own files in these directories
            continue
        #if "DIR" in directory:
        #    np.save(f"extracted_info/DIR.npy", np.zeros(7)) # 4 energies, phi angle, psi angle, CA_CA distance
        #    continue
        # get energies of structure
        with open(f"{directory}/{filename}",'r') as f, open(f"{directory}/also_bad.txt", 'a') as bad_log:
            for counter, line in enumerate(f):
                if counter < 4:
                    try:
                        conf_info.append(float(line.split("=")[1].strip()))
                        # index 0: water absolute energy
                        # index 1: water energy relative to global min for the current sequence
                        # index 2: octanol absolute energy
                        # index 3: octanol energy relative to global min for the current sequence
                    except:
                        bad_log.write(f"bad: {filename}\n")
                        conf_info = [0,0,0,0,] # signal that this conf should be ignored later in the analysis pipeline
                        break
                else:
                    break
        # get phi and psi dihedral angles of middle residue in structure
        u = mda.Universe(f"{directory}/{filename}")
        #selection = u.residues[2] # 5 residues total because of caps, position[2] is the middle residue
        rama_output = dihedrals.Ramachandran(u.atoms, c_name='C', n_name='N', ca_name='CA',check_protein=True).run()
        try:
            phi = rama_output.angles[0,2,0] # axis 0 is frame, axis 1 is residue position in selection, axis 2 is angle type
        except IndexError:
            print(f"BAD phi: {directory}, {filename}")
            phi = 0 # we'll correct it later
        try:
            psi = rama_output.angles[0,2,1]
        except IndexError:
            print(f"BAD psi: {directory}, {filename}")
            psi = 0 # just a placeholder; we'll fix it later
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


