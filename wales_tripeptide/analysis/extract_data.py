import numpy as np
import matplotlib.pyplot as plt
import os
if not os.path.exists('stats_by_seq'): os.mkdir('stats_by_seq')

# constants                                                  # note that we exclude HIP
aa = ["ALA",'ARG','ASN','ASP','CYS','GLN','GLU','GLY','HID','HIE','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
num_chi_angles = {"ALA":0, "ARG":5, "ASN":2, "ASP":2, "CYS":1, "GLN":3, "GLU":3, "GLY":0,
                  "HID":2, "HIE":2,          "ILE":2, "LEU":2, "LYS":4, "MET":3, "PHE":2,
                  "PRO":0, "SER":1, "THR":1, "TRP":2, "TYR":2, "VAL":1}

def make_scatterplot(phi, psi, energy, filename):
    scatterplot = plt.scatter(phi, psi, c=energy, vmin=0, s=1)
    plt.colorbar(scatterplot)
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.savefig(filename)
    plt.close()

# load data
energy = []
middle_phi = []
middle_psi = []
energy_HID = [] 
energy_HIE = [] 
HID_middle_phi = []
HID_middle_psi = [] 
HIE_middle_phi = [] 
HIE_middle_psi = [] 
for a1 in aa:
    print(f"Extracting data for peptides whose first residue is {a1}")
    middle_columns_start = 1+2+num_chi_angles[a1] 
                        # we have to go through (energy + (phi1 + psi1) + {set of all chi1}) columns before we get to the middle residue
    for a2 in aa:
        if a2 in ["ALA", "GLY"]:
            continue # sequences with ALA and GLY in position 2 were not considered in the study
        for a3 in aa:
            energy_seq = []
            middle_phi_seq = []
            middle_psi_seq = []
            #with open(f"../../../wales_tripeptide/TripeptideConformations/{a1}_{a2}_{a3}_new.csv",'r') as f:
            with open(f"../TripeptideConformations/{a1}_{a2}_{a3}_new.csv",'r') as f:
                for line in f:
                    if not line.strip():
                        continue# just an empty line
                    conf_info = [item for item in line.split(" ") if item != ""]
                    try:
                        energy_seq.append(float(conf_info[0]))
                    except IndexError:
                        print(f"conf_info: {conf_info}")
                        print(f"sequence: {a1}{a2}{a3}")
                        raise
                    middle_phi_seq.append(float(conf_info[middle_columns_start]))
                    try:
                        middle_psi_seq.append(float(conf_info[middle_columns_start+1]))
                    except IndexError: 
                        print(f"sequence: {a1}{a2}{a3}")
                        print(conf_info)
                        print(f"middle_columns_start: {middle_columns_start}")
                        print(f"num_chi_angles[a1]: {num_chi_angles[a1]}")
                        raise
            energy_seq = np.array(energy_seq)
            middle_phi_seq = np.array(middle_phi_seq)*np.pi/180
            middle_psi_seq = np.array(middle_psi_seq)*np.pi/180
            np.save(f"stats_by_seq/{a1}{a2}{a3}.npy",np.array([middle_phi_seq, middle_psi_seq, energy_seq]))
            if a2 == "HID":
                energy_HID.append(energy_seq)
                HID_middle_phi.append(middle_phi_seq)
                HID_middle_psi.append(middle_psi_seq)
                # we're going to exclude HID from the main energy, phi, and psi lists so we don't count histidine twice
            elif a2 == "HIE":
                energy_HIE.append(energy_seq)
                HIE_middle_phi.append(middle_phi_seq)
                HIE_middle_psi.append(middle_psi_seq)
                energy.append(energy_seq)
                middle_phi.append(middle_phi_seq)
                middle_psi.append(middle_psi_seq)
            else:
                energy.append(energy_seq)
                middle_phi.append(middle_phi_seq)
                middle_psi.append(middle_psi_seq)
            # we could also do separate plots for cases where the first residues is HID/HIE but not necessarily the second or third residues,
            # or the third residue is HID/HIE, but not necessarily the first or second residues

# plot energy vs ramachandran angles for HID and HIE
#     doing all sequences where HID or HIE is the middle residue
#     need to flatten arrays before inputting (this means that every conformer of every qualifying sequence is included on the same plot)
HID_scatter_phi = []
HID_scatter_psi = []
HID_scatter_energy = []
HIE_scatter_phi = []
HIE_scatter_psi = []
HIE_scatter_energy = []
assert len(HID_middle_phi) == len(HID_middle_psi) == len(energy_HID) == len(HIE_middle_phi) == len(HIE_middle_psi) == len(energy_HIE)
for counter in range(len(HID_middle_phi)):
    assert len(HID_middle_phi[counter])==len(HID_middle_psi[counter])==len(energy_HID[counter])
    assert len(HIE_middle_phi[counter])==len(HIE_middle_psi[counter])==len(energy_HIE[counter])      
    for counter2 in range(len(HID_middle_phi[counter])):
        HID_scatter_phi.append(HID_middle_phi[counter][counter2])
        HID_scatter_psi.append(HID_middle_psi[counter][counter2])
        HID_scatter_energy.append(energy_HID[counter][counter2])
    for counter2 in range(len(HIE_middle_phi[counter])):
        HIE_scatter_phi.append(HIE_middle_phi[counter][counter2])
        HIE_scatter_psi.append(HIE_middle_psi[counter][counter2])
        HIE_scatter_energy.append(energy_HIE[counter][counter2])
make_scatterplot(HID_scatter_phi, HID_scatter_psi, HID_scatter_energy, "HID_scatterplot.png")
make_scatterplot(HIE_scatter_phi, HIE_scatter_psi, HIE_scatter_energy, "HIE_scatterplot.png")
 
