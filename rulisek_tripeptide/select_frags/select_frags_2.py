import numpy as np

# user input
protein = '2ppn'
fasta = f'{protein}.fasta'
file_path = "/work/pw8/fc36/rulisek_tripeptide/extracted_info"
energy_type ='MD_and_water'#'protein_and_water'#'MD_and_water'
output_angles_name = f'{protein}_memory_angles_{energy_type}.npy'
output_gammas_name = f'{protein}_memory_weights_2_{energy_type}.npy'

# constants
#     first word: where the conformers were sampled from
#     second word: solvent used for dielectric screening and solvent interaction calculation
energy_type_to_index = {'MD_and_water': 0,
                        'MD_and_octanol': 1,
                        'protein_and_water': 2,
                        'protein_and_octanol': 3}

# load sequence
lines = []
with open(fasta,'r') as f:
    for line in f:
        lines.append(line)
if not len(lines)==2:
    raise ValueError(f"there should be 2 lines in the fasta file, the first line being a comment, beginning with '>', and the second line being sequence without whitespace")
if not lines[0][0] == ">":
    raise ValueError(f"there should be 2 lines in the fasta file, the first line being a comment, beginning with '>', and the second line being sequence without whitespace")
sequence = lines[1].strip() # we said there shouldn't be whitespace but just to be nice we'll do this

# load angles and energies for each position (second to second-to-last residues, considering self and nearest neighbors)
phi_by_pos = []
psi_by_pos = []
gammas_by_pos = []
for counter in range(1,len(sequence)-1):
    # our arrays that we load are formatted in a particular way -- see ../extract_info*.py
    data = np.load(f"{file_path}/{sequence[counter-1:counter+2]}.npy")
    phi_angles = data[:,4]*np.pi/180
    psi_angles = data[:,5]*np.pi/180
    energies = data[:,energy_type_to_index[energy_type]]    
    # some conformers don't have energies printed for some reason
    #     earlier in my data processing pipeline, I set those energies to 0
    #     now we remove them
    lessthanzero_indices = np.where(energies < 0)
    notzero_indices = np.where(energies != 0)
    assert np.all(np.array(lessthanzero_indices) == np.array(notzero_indices))
    keep_indices = notzero_indices
    phi_angles = phi_angles[keep_indices]
    psi_angles = psi_angles[keep_indices]
    gammas = energies[keep_indices]
    # there are a lot of regions on the ramachandran plot that we didn't sample
    #     these regions probably tend to be about as high in energy as the highest-energy conformer we sampled
    #     but we don't explicitly account for unsampled regions
    #     so we'll set the energy of the highest-energy sampled conformer to 0, and reference everything to that
    #     so that everything we've sampled can be given a favorable gaussian well where the depth is determined by the energy of the conformer
    gammas -= np.max(gammas)
    # there are some weird high-energy outlier conformers in the original dataset that were excluded from the analysis in that paper;
    # we'll exclude weirdly high-energy conformers, too.
    if np.min(gammas) < -1000: # this indicates that the max is a big outlier
        keep_indices = np.where(gammas<-1000)
        phi_angles = phi_angles[keep_indices]
        psi_angles = psi_angles[keep_indices]
        gammas = gammas[keep_indices]
        gammas -= np.max(gammas)
    # now we add our info for the current structure to our big list of info for all structures
    phi_by_pos.append(phi_angles)
    psi_by_pos.append(psi_angles)
    gammas_by_pos.append(gammas)
    #breakpoint()

# rearrange angles and gammas into the format expected by the rama_AM_term
#     identify position with greatest number of memories--we'll need to pad the others with 0s to reach this length
longest = 0
for phi, psi, gammas in zip(phi_by_pos, psi_by_pos, gammas_by_pos):
    assert phi.shape == psi.shape == gammas.shape
    if phi.shape[0] > longest:
        longest = phi.shape[0]
#     initialize new array with the correct dimensions and add elements wherever they belong
openawsem_format_angles = np.zeros((len(phi_by_pos), 2, longest))
openawsem_format_gammas = np.zeros((len(phi_by_pos), longest))
#breakpoint()
for counter, (phi, psi, gammas) in enumerate(zip(phi_by_pos, psi_by_pos, gammas_by_pos)):
    openawsem_format_angles[counter, 0, :phi.shape[0]] = phi
    openawsem_format_angles[counter, 1, :psi.shape[0]] = psi
    openawsem_format_gammas[counter, :gammas.shape[0]] = gammas #/ np.sum(gammas) # normalize gammas cause the magnitudes are meaningless anyway

# save arrays
np.save(output_angles_name, openawsem_format_angles)
np.save(output_gammas_name, openawsem_format_gammas)

