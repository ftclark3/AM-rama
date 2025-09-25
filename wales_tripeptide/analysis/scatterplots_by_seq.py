import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    # causes new axes to be initiated in each subprocess so that we don't mix our plots together
import multiprocessing

one_to_three = {'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 'G':'GLY', 'H':'HID', 'I':'ILE',
                'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN', 'P':'PRO', 'Q':'GLN', 'R':'ARG', 'S':'SER',
                'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'}

def loop_interior(a1, a2, a3):
    rulisek_path = f'/work/pw8/fc36/rulisek_tripeptide/extracted_info/{a1}{a2}{a3}.npy' 
        # these have shape (num_conformers, 7), where the energy is column 0, phi is column 4, and psi is column 5
    wales_path = f'/work/pw8/fc36/wales_tripeptide/analysis/stats_by_seq/{one_to_three[a1]}{one_to_three[a2]}{one_to_three[a3]}.npy'
    wales_data = np.load(wales_path)
        # these have shape (3, num_conformers), where the first row is phi, second row is psi, and third is energy
    rulisek_data = np.load(rulisek_path)[:,[0,4,5]]
    rulisek_data[:,0] -= np.min(rulisek_data[:,0]) # set rulisek minimum to 0
    #print(np.max(wales_data[2,:]))
    #print(np.median(wales_data[2,:]))
    #print(np.median(rulisek_data[:,0]))
    #print((np.max(wales_data[2,:])/np.median(wales_data[2,:]))*np.median(rulisek_data[:,0]))
    rulisek_data = rulisek_data[rulisek_data[:,0]<(np.max(wales_data[2,:])/np.median(wales_data[2,:]))*np.median(rulisek_data[:,0]),:] # remove outliers
    plt.scatter(rulisek_data[:,1]*np.pi/180,rulisek_data[:,2]*np.pi/180,c=rulisek_data[:,0],cmap="Blues_r", s=4)
    plt.colorbar().set_label('Rulisek')
    plt.scatter(wales_data[0,:],wales_data[1,:],c=wales_data[2,:],cmap="Greens_r", s=4)
    plt.colorbar().set_label('Wales')
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.title(f"Rulisek vs. Wales Dihedrals and Energies, {a1}{a2}{a3}")
    plt.savefig(f'scatterplots_by_seq/{a1}{a2}{a3}.png')
    plt.close()

args = []
for a1 in one_to_three.keys():
    for a2 in one_to_three.keys():
        if a2 in ["G","A"]:
            continue # wales study didn't consider these amino acids as the middle residue
        for a3 in one_to_three.keys():
            args.append([a1,a2,a3])
#for arg in args:
#    loop_interior(arg[0],arg[1],arg[2])
#    exit()
with multiprocessing.Pool() as pool:
    pool.starmap(loop_interior, args)
