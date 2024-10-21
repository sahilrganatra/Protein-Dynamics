import sys
from prody import *
from matplotlib.pylab import *

ion() # turn interactive mode on

def gnm_flux(pdb_id):

    prody.LOGGER._setverbosity('none')

    # download the PDB file using parsePDB
    structure = parsePDB(pdb_id)
    if structure is None:
        print(f'Failed to download or parse PDB file: {pdb_id}')
        return

    # select C-alpha atoms
    calphas = structure.select('calpha')
    if len(calphas) == 0:
        print('Failed to select C-alpha atoms.')
        return

    # obtain C-alpha coordinates to pass to buildKirchhoff() method
    calpha_coords = calphas.getCoords()
    
    # perform GNM analysis
    gnm = GNM(calphas)

    # cutoff = max distance at which two atoms are considered to interact (in Ångströms)
    gnm.buildKirchhoff(calpha_coords, cutoff = 10.0, gamma = 1.0)
    
    # calculate normal modes
    gnm.calcModes()

    # calculate squared fluctuations
    squared_flucts = prody.calcSqFlucts(gnm) # pass normal modes
    if squared_flucts is None:
        print('Failed to calculate squared fluctuations.')
        return

    max_gnm_flux = round(max(squared_flucts), 3) # round final value to 3 decimal points

    sys.stdout.write(f'max gnm sqfluct: {max_gnm_flux}')

# take an input argument (protein from PDB) and run the script
pdb_id = sys.argv[1]
gnm_flux(pdb_id)