import numpy as np
import os
import multiprocessing
import copy

names = [f"extracted_info/{filename}" for filename in os.listdir('extracted_info')]
names = [names[0]]
arrays = copy.deepcopy(names)
for counter in range(len(arrays)):
    try:
        arrays[counter] = np.load(arrays[counter])
        #print("success")
    except:
        raise
        #print(f"fail: {arrays[counter]}")

def compute_gamma(array,name):
    w = 0.2 # width parameter
    well_energies = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0} # real well energies will be very negative, so we just initialize to 0 and will replace it later
    well_coordinates = {1:(180,180), 2:(180,180), 3:(180,180), 4:(180,180), 5:(180,180), 6:(180,180), 7:(180,180), 8:(180,180), 9:(180,180)} # coordinates corresponding to the minimum energies in well_energies
    indicators_at_coords_of_lowest_energy = {1:[0,0,0,0,0,0,0,0,0], 2:[0,0,0,0,0,0,0,0,0], 3:[0,0,0,0,0,0,0,0,0], 4:[0,0,0,0,0,0,0,0,0], 5:[0,0,0,0,0,0,0,0,0], 6:[0,0,0,0,0,0,0,0,0], 7:[0,0,0,0,0,0,0,0,0], 8:[0,0,0,0,0,0,0,0,0], 9:[0,0,0,0,0,0,0,0,0]} # all indicator functions evaluated at the best (phi, psi) pair for each well, "best" defined as the one that gives the lowest P-CONF energy; it is expected that all but 1 indicator in these lists will be almost 0
    for row in array:
        if np.all(row[:4] == np.array([0,0,0,0])):
            print("row skipped")
            continue # no energies were logged
        else:
            phi = row[4]*np.pi/180
            psi = row[5]*np.pi/180
            energy = row[0] # absolute free energy in water
            ind_1 = np.exp(-(w*10*np.square(np.cos(phi-(-145*np.pi/180))-1)+w*400*np.square(np.cos(psi-(160*np.pi/180))-1))) + np.exp(-(w*4500*np.square(np.cos((phi+psi)/2-(10*np.pi/180))-1)+w*2300*np.square(np.cos((phi-psi)/2-(-120*np.pi/180))-1)))
            ind_2 = np.exp(-(w*3200*np.square(np.cos((phi+psi)/2-(50*np.pi/180))-1)+w*400*np.square(np.cos((phi-psi)/2-(-90*np.pi/180))-1)))
            ind_3 = np.exp(-(w*1600*np.square(np.cos((-2*phi+psi)/2-(130*np.pi/180))-1)+w*800*np.square(np.cos((-2*phi-psi)/2-(-50*np.pi/180))-1)))
            ind_4 = np.exp(-(w*3200*np.square(np.cos((phi+psi)/2-(-50*np.pi/180))-1)+w*100*np.square(np.cos((phi-psi)/2-(-50*np.pi/180))-1)))
            ind_5 = np.exp(-(w*30*np.square(np.cos(phi-(-140*np.pi/180))-1)+w*400*np.square(np.cos(psi-(-50*np.pi/180))-1)))
            ind_6 = np.exp(-(w*3200*np.square(np.cos((phi+psi)/2-(-50*np.pi/180))-1)+w*400*np.square(np.cos((phi-psi)/2-(-0*np.pi/180))-1)))
            ind_7 = np.exp(-(w*3200*np.square(np.cos((phi+psi)/2-(-50*np.pi/180))-1)+w*3200*np.square(np.cos((phi-psi)/2-(90*np.pi/180))-1)))
            ind_8 = np.exp(-(w*1600*np.square(np.cos((-2*phi+psi)/2-(-100*np.pi/180))-1)+w*800*np.square(np.cos((-2*phi-psi)/2-(-35*np.pi/180))-1)))
            ind_9 = np.exp(-(w*8000*np.square(np.cos((phi+psi)/2-(50*np.pi/180))-1)+w*300*np.square(np.cos((phi-psi)/2-(-10*np.pi/180))-1)))
            indicators_list = [ind_1, ind_2, ind_3, ind_4, ind_5, ind_6, ind_7, ind_8, ind_9]
            if max(indicators_list) < 0.01: # ignore outlier point
                continue
            else:
                well_assignment = indicators_list.index(max(indicators_list))+1 # assign this sample to the well whose indicator is most-activated
                if energy < well_energies[well_assignment]:
                    well_energies[well_assignment] = energy # update to more favorable energy
                    well_coordinates[well_assignment] = (phi, psi) # store coordinates where the more favorable energy occurs
                    indicators_at_coords_of_lowest_energy[well_assignment] = indicators_list # going from 1-indexed to 0-indexed
    #X = np.array(well_coordinates.values()).reshape((-1,2)) # column 0 should be phi, column 1 should be psi
    #X = np.stack(np.ones(X.shape[0]), X, axis=1) # add constant coefficient to be multiplied by bias parameter
    X = np.array(list(indicators_at_coords_of_lowest_energy.values()))
    assert X.shape == (9,9), X # 9 wells, 9 indicator classes evaluated for each well
    y = np.array(list(well_energies.values())) # prepend 1 to list of values for bias variable
    y -= np.mean(y) # all the energies are very negative; we want to subtract that baseline bias energy
    assert y.shape == (9,) # 9 entries in well_energies, one for each well
    try:
        #gamma, residuals, rank, singular_values = np.linalg.lstsq(X,y)
        gamma = np.linalg.inv(X) @ y
    except: # matrix is nearly singular probably
        print(f"BAD: {name}")
        gamma = y/np.diag(X)
    print("location")
    print(name)
    print(gamma)
    print(y/np.diag(X)) # the result we would get if we ignored linear between energy of a given well and weights assigned to other wells-- should be about the same as gamma
    assert np.allclose(X@gamma, y) # just in case we had numerical issues with our matrix inversion or i made a basic coding error
    exit()
    return gamma

with multiprocessing.Pool() as pool:
    all_gammas = pool.starmap(compute_gamma, zip(arrays, names))
    all_gammas = np.array(all_gammas).reshape((20,20,20))
    np.save('all_gammas.npy',all_gammas)
