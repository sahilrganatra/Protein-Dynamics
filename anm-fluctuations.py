# same format as gnm-flux, except we construct a Hessian matrix

import sys
from prody import *

def anm_flux(pdb_id):

    prody.LOGGER._setverbosity('none')

    structure = parsePDB(pdb_id)
    if structure is None:
        print(f'Failed to download or parse structure: {pdb_id}')
        return

    calphas = structure.select('calpha')
    if len(calphas) == 0:
        print('Unable to identify alpha carbons in this structure.')
        return
    calpha_coords = calphas.getCoords()
    
    # perform ANM analysis
    anm = ANM(calphas)
    
    anm.buildHessian(calpha_coords, cutoff = 15.0, gamma = 1.0)
    
    anm.calcModes()

    squared_flucts = prody.calcSqFlucts(anm)
    if squared_flucts is None:
        print('Failed to calculate squared fluctuations.')
        return

    max_anm_flux = round(max(squared_flucts), 3)

    print(f'Max ANM SqFluct: {max_anm_flux}')

# run the script
pdb_id = sys.argv[1]
anm_flux(pdb_id)