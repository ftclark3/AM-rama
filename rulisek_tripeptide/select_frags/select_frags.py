import numpy as np

# user input
fasta = '1r69.fasta'
file_path = "/work/pw8/fc36/rulisek_tripeptide/extracted_info"
energy_type = 'MD_and_water'#'protein_and_water'#'MD_and_water'
kT = 300*8.31446261815324/1000/4.184  # assuming energy is in units of kcal/mol (see below)
kT *= 2 # errors in the energies estimated to be about 2 kcal/mol and we don't want this noise to overwhelm our signal, so we'll keep more guys
output_angles_name = f'1r69_memory_angles_{energy_type}.npy'
output_gammas_name = f'1r69_memory_weights_{energy_type}.npy'

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
    data = np.load(f"{file_path}/{sequence[counter-1:counter+2]}.npy")
    phi_angles = data[:,4]
    psi_angles = data[:,5]
    #import pdb; pdb.set_trace()
    energies = data[:,energy_type_to_index[energy_type]]
    energies -= np.min(energies) # set the minimum to 0
    exponent = -energies/kT
    keep_indices = np.where(exponent >= -4) # anything else will be such a small weight it won't really do anything
    phi_angles = phi_angles[keep_indices]
    psi_angles = psi_angles[keep_indices]
    exponent = exponent[keep_indices]
    gammas = np.exp(exponent)# instead of machine learning the gammas, we'll just get the boltzmann weights of the physical energies
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
    openawsem_format_gammas[counter, :gammas.shape[0]] = gammas / np.sum(gammas) # divide by sum of boltzmann weights to remain comparable to SM version

# save arrays
np.save(output_angles_name, openawsem_format_angles)
np.save(output_gammas_name, openawsem_format_gammas)

