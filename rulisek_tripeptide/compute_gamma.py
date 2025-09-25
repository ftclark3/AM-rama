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
    X = []
    y = []
    for row in array:
        if np.all(row[:4] == np.array([0,0,0,0])):
            print("row skipped")
            continue # no energies were logged
        else:
            y.append(row[0]) # absolute energy in water
            phi = row[4]*np.pi/180
            psi = row[5]*np.pi/180
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
            if max(indicators_list) < 0.01:
                #pass
                y = y[:-1] # remove the last element that we just added because we're not considering this conformer
                continue
            else:
                X.append([1,ind_1, ind_2, ind_3, ind_4, ind_5, ind_6, ind_7, ind_8, ind_9]) # 1 for bias term
    X = np.array(X)
    y = np.array(y)
    try:
        gamma, residuals, rank, singular_values = np.linalg.lstsq(X,y)
    except:
        print(f"BAD: {name}")
        return (0, 0, 0) # just placeholders
    print("location")
    print(name)
    print(f"shape of y: {y.shape}")
    print(f"shape of array: {array.shape}")
    print(gamma)
    print(np.mean(np.sqrt(residuals)/np.abs(y-gamma[0])), np.max(np.sqrt(residuals)/np.abs(y-gamma[0])), np.median(np.sqrt(residuals)/np.abs(y-gamma[0]))) # i think this is the average unsigned error
    print(np.max(residuals/np.abs(y-gamma[0])))
    return gamma, np.mean(np.sqrt(residuals)/np.abs(y-gamma[0])), np.max(residuals/np.abs(y-gamma[0]))

with multiprocessing.Pool() as pool:
    all_gammas, normMSE, max_sq_rel_error = pool.starmap(compute_gamma, zip(arrays, names))
    all_gammas = all_gammas.reshape((20,20,20))
    normMSE = normMSE.reshape((20,20,20))
    max_sq_rel_error = max_sq_rel_error.reshape((20,20,20))
    np.save('all_gammas.npy',all_gammas)
    np.save('normMSE.npy',normMSE)
    np.save('max_sq_rel_error.npy',max_sq_rel_error)
