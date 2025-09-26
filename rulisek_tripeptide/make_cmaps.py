import numpy as np
import os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt; plt.rcParams.update({'font.size': 22})  
                       # causes new axes to be initiated in each subprocess so that we don't mix our plots together 
import multiprocessing

def evaluate_potential(phis, psis, gammas, sigma=np.pi/180, min_energy=10):
    grid=np.zeros((90,90)) # width of 1 degree, min seq to 10 kcal/mol
    for phi_coord in range(grid.shape[0]):
        phi_coord_adjusted = 4*phi_coord # each increment of the counter by 1 should equal 4 degrees
        phi_coord_rad = phi_coord_adjusted * np.pi/180 - np.pi # shift from (0, 2pi) to (-pi, pi)
        phi_part = (np.cos(phis - phi_coord_rad) - 1)**2
        for psi_coord in range(grid.shape[1]):
            psi_coord_adjusted = 4*psi_coord # each increment of the counter by 1 should equal 4 degrees
            psi_coord_rad = psi_coord_adjusted * np.pi/180 - np.pi  # shift from (0, 2pi) to (-pi, pi)
            psi_part = (np.cos(psis - psi_coord_rad) - 1)**2
            grid[phi_coord,psi_coord] += np.sum(gammas*np.exp(-(phi_part+psi_part)/(2*sigma**2)))
    grid = grid / np.abs(np.min(grid)) * 10 # minimum set to -10 kcal/mol
    return grid

# setup
energy_type ='MD_and_water'#'MD_and_octanol'#'MD_and_water'
if not os.path.exists('cmaps'):
    os.mkdir('cmaps')
if not os.path.exists(f'cmaps/data_{energy_type}'):
    os.mkdir(f'cmaps/data_{energy_type}')
if not os.path.exists(f'cmaps/plots_{energy_type}'):
    os.mkdir(f'cmaps/plots_{energy_type}')
file_path = "extracted_info"
seqs = []
alphabet = "ACDEFGHIKLMNPQRSTVWY"
for aa1 in alphabet:
    for aa2 in alphabet:
        for aa3 in alphabet:
            seqs.append(f"{aa1}{aa2}{aa3}")
# data array indices
#     first word: where the conformers were sampled from
#     second word: solvent used for dielectric screening and solvent interaction calculation
energy_type_to_index = {'MD_and_water': 0,
                        'MD_and_octanol': 2}
# I THINK THIS LIST WAS WRONG
#{'MD_and_water': 0,
                       # 'MD_and_octanol': 1,
                       # 'protein_and_water': 2,
                       # 'protein_and_octanol': 3}

# load angles and energies for each position (second to second-to-last residues, considering self and nearest neighbors)
phi_by_seq= []
psi_by_seq = []
gammas_by_seq = []
for seq in seqs:
    # our arrays that we load are formatted in a particular way -- see ../extract_info_4.py
    data = np.load(f"{file_path}/{seq}.npy")
    phi_angles = data[:,4]*np.pi/180
    psi_angles = data[:,5]*np.pi/180
    gammas = data[:,energy_type_to_index[energy_type]]    
    breakpoint()
    # there are a lot of regions on the ramachandran plot that we didn't sample
    #     these regions probably tend to be about as high in energy as the highest-energy conformer we sampled
    #     but we don't explicitly account for unsampled regions
    #     so we'll set the energy of the highest-energy sampled conformer to 0, and reference everything to that
    #     so that everything we've sampled can be given a favorable gaussian well where the depth is determined by the energy of the conformer 
    gammas -= np.max(gammas)
    # now we add our info for the current structure to our big list of info for all structures
    phi_by_seq.append(phi_angles)
    psi_by_seq.append(psi_angles)
    gammas_by_seq.append(gammas)

# convert angles and energies into gaussian potential
def loop_body(phis, psis, gammas, seq):
    assert phis.shape == psis.shape == gammas.shape
    assert len(phis.shape) == 1
    grid = evaluate_potential(phis, psis, gammas)
    fig, ax = plt.subplots()
    im = ax.imshow(grid.T,origin='lower',interpolation='none') # now, phi increases along the x axis (axis 1) and psi increases going up the y axis (axis 0)
    cbar = fig.colorbar(im, ax=ax, ticks=[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0])#cbar = fig.colorbar(matplotlib.cm.ScalarMappable(), ax=ax, ticks=[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0])#plt.colorbar(im,label='Energy (kcal/mol)')
    cbar.ax.set_yticklabels(['-10','-9','-8','-7','-6','-5','-4','-3','-2','-1','0'])
    ax.set_xlabel('phi')#plt.xlabel('phi')
    ax.set_ylabel('psi')#plt.ylabel('psi')
    #ax.set_xticks([0, 90, 180, 270, 360],labels=[-180,-90,0,90,180])#plt.xticks([0, 90, 180, 270, 360],labels=[-180,-90,0,90,180])
    #ax.set_yticks([0, 90, 180, 270, 360],labels=[-180,-90,0,90,180])#plt.yticks([0, 90, 180, 270, 360],labels=[-180,-90,0,90,180])
    ax.set_xticks([0, 22.5, 45, 67.5, 90,],labels=[-180,-90,0,90,180])#plt.xticks([0, 90, 180, 270, 360],labels=[-180,-90,0,90,180])
    ax.set_yticks([0, 22.5, 45, 67.5, 90,],labels=[-180,-90,0,90,180])#plt.yticks([0, 90, 180, 270, 360],labels=[-180,-90,0,90,180])
    ax.set_title(f"Energy vs. phi/psi for {seq[1]}\n in context {seq}")#plt.title(f"Energy vs. phi/psi for {seq[1]}\n in context {seq}")
    fig.savefig(f'cmaps/plots_{energy_type}/{seq}.png',bbox_inches='tight')#plt.savefig(f'cmaps/plots/{seq}.png',bbox_inches='tight')
    plt.close(fig) # this is the command to close the figure
    # 
    # openmm wants the array to go from 0 to 2pi, but currently it goes from -pi to pi, so we have to fix that
    grid = grid.T[::-1,:] # now, phi increases along the x axis (axis 1) and psi increases going up the y axis (axis 0)
    old_ul = grid[:45,:45]
    old_ur = grid[:45,45:]
    old_ll = grid[45:,:45]
    old_lr = grid[45:,45:]
    new_grid = np.zeros((90,90))
    new_grid[45:,:45] = old_ur # wrap [0,pi)x[0,pi) to [0,pi)x[0,pi), so it stays put
    new_grid[45:,45:] = old_ul # wrap [-pi,0)x[0,pi] to [pi,2pi)x[0,pi)
    new_grid[:45,:45] = old_lr # wrap [0,pi)x[-pi,0) to [0,pi)x[pi,2pi)
    new_grid[:45,45:] = old_ll # wrap [-pi,0)x[-pi,0) to [pi,2pi)x[pi,2pi)
    grid = new_grid[::-1,:].T # put grid back as we found it, with phi increasing down axis 0 and psi increasing across axis 1
    grid = grid.flatten(order='F') # This is going to be passed in to OpenMM's CMAPTorsionForce.
                                   # CMAPTorsionForces expects a 1D array where array[i+size*j] is the energy when
                                   # the "first" angle is at position i and the "second" angle is at position j.
                                   # We normally think of phi as being the first angle and psi as being the second angle,
                                   # so we should flatten our array such that phi varies within blocks of length size
                                   # and psi varies between blocks of length size. This requires flattening our array
                                   # in column-major ("F", for "fortran") order
    np.save(f'cmaps/data_{energy_type}/{seq}.npy', grid)
    
#for phis, psis, gammas, seq in zip(phi_by_seq, psi_by_seq, gammas_by_seq, seqs):
#    loop_body(phis, psis, gammas, seq)
with multiprocessing.Pool() as pool:
    pool.starmap(loop_body, zip(phi_by_seq, psi_by_seq, gammas_by_seq, seqs))
