#small script to evaluate the energy differences between octanol and water surfaces

import os
import numpy as np

for filename in os.listdir('rulisek_tripeptide/cmaps/data_MD_and_water'):
    water = np.load(f'rulisek_tripeptide/cmaps/data_MD_and_water/{filename}')
    octanol = np.load(f'rulisek_tripeptide/cmaps/data_MD_and_octanol/{filename}')
    foo = octanol - water
    print(f'{filename}: {np.max(np.abs(foo))}')
