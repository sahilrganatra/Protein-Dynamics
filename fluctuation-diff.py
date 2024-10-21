import sys
from prody import *
import numpy as np
import matplotlib.pyplot as plt

def fluct_diff(pdb_id):

    prody.LOGGER._setverbosity('none')

    # download PDB files using parsePDB
    structure = parsePDB(pdb_id)
    if structure is None:
        print(f'Failed to parse or download PDB file: {pdb_id}')
    
    # find C-alphas
    calphas = structure.select('calpha')
    if len(calphas) == 0:
        print('Failed to idenfity alpha carbons in this structure.')
        return
    calpha_coords = calphas.getCoords()

    # perform GNM analysis first
    gnm = GNM(calphas)
    gnm.buildKirchhoff(calpha_coords, cutoff = 10.0, gamma = 1.0)
    gnm.calcModes()
    gnm_squared_flucts = prody.calcSqFlucts(gnm)

    # perform ANM analysis next
    anm = ANM(calphas)
    anm.buildHessian(calpha_coords, cutoff = 15.0, gamma = 1.0)
    anm.calcModes()
    anm_squared_flucts = prody.calcSqFlucts(anm)

    # plot squared fluctuations, GNM vs. ANM
    plt.figure(figsize = (10, 7))
    plt.plot(gnm_squared_flucts, label = 'GNM', lw = 2, color = 'purple')
    plt.plot(anm_squared_flucts, label = 'ANM', lw = 2, color = 'orange')
    plt.xlabel('Residue Number')
    plt.ylabel('Squared Fluctuation')
    plt.title('Cα Squared Fluctuation - GNM vs. ANM Analysis')
    plt.grid()
    plt.legend()
    plt.show()

    # compute abs. value of the difference in GNM vs. ANM for each atom, print largest difference
    fluct_differences = np.abs(gnm_squared_flucts - anm_squared_flucts)
    max_diff = round(max(fluct_differences), 3)
    print(f'Max Abs Difference: {max_diff}')

    # normalize GNM and ANM values such that max value of each is 1.0
    normalized_gnm = gnm_squared_flucts / np.max(gnm_squared_flucts)
    normalized_anm = anm_squared_flucts / np.max(anm_squared_flucts)

    # plot normalized squared flucts
    plt.figure(figsize = (10, 7))
    plt.plot(normalized_gnm, label = 'GNM', lw = 2, color = 'purple')
    plt.plot(normalized_anm, label = 'ANM', lw = 2, color = 'orange')
    plt.xlabel('Residue Number')
    plt.ylabel('Normalized Squared Fluctuation')
    plt.title('Normalized Cα Squared Fluctuation - GNM vs. ANM Analysis')
    plt.grid()
    plt.legend()
    plt.show()

    # compute abs. difference between normalized squared fluctuations
    normalized_diff = np.abs(normalized_gnm - normalized_anm)
    max_norm_diff = round(max(normalized_diff), 3)
    print(f'Max Abs Norm Difference: {max_norm_diff}')

# run the script
pdb_id = sys.argv[1]
fluct_diff(pdb_id)