import os
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from timeit import default_timer as timer

# TO DO:
# Consider PBC

def get_triplets(solvation_shell):
    '''Returns list of all triplets of nearest neighbours water molecules in the solvation shell, corresponding to each plane of the octahedron.
    Should work also if shell is less than 6.'''
    triplets=[]
    for comb in combinations(range(len(solvation_shell)),3):    # For each possible triplet
        #print(comb)
        pairs = combinations(comb, 2)
        for i,j in pairs:
            if np.linalg.norm(solvation_shell[i].position - solvation_shell[j].position) > 4:   # Checks that distances between water molecules is < 4, meaning they are nearest neighbours
                break
        else:
            triplets.append(comb)

    return triplets

def compute_versors(triplet):
    '''Returns versor normal to the plane corresponding to the triplet'''
    a = solvation_shell[triplet[0]].position - solvation_shell[triplet[1]].position
    b = solvation_shell[triplet[1]].position - solvation_shell[triplet[2]].position

    vector = np.cross(a,b)
    vector /= np.linalg.norm(vector)    # Normalizazion

    return vector

def autocorrelation(versor, k):
    '''Returns temporal autocorrelation of versor evaluated at k. k is tau/delta_t'''
    N = versor.shape[1]     # Lenght of simulation

    acf = 0
    for i in range(N-k):
        acf += np.dot(versor[:,i],versor[:,i+k])
    return acf/(N-k)

timer_start = timer()

working_dir = os.path.dirname(__file__)
XTC = os.path.join(working_dir, 'hybr+6MN+36wat.0_200ps.xtc')
TPR = os.path.join(working_dir, 'hybr+6MN+36wat.tpr')
directory_figures = os.path.join(working_dir, 'Figures')
directory_arrays = 'racf_arrays'

#delta_t = 1 # ps

u = mda.Universe(TPR, XTC)

water_molecules = u.select_atoms('type OW')
mn_ions = u.select_atoms('resname MN')

for i,mn in enumerate(mn_ions):
    # Get solvation shell
    solvation_shell = [water for water in water_molecules if np.linalg.norm(mn.position - water.position) < 3]  # All OW that are closer than 3 to the Mn

    # Get triplets
    triplets = get_triplets(solvation_shell)
    versors = {triplet : np.array(compute_versors(triplet)).reshape((3,1)) for triplet in triplets} # Dictionary: versors[triplet] gives versor corresponding to plane

    for frame in u.trajectory[1:]:
        # Calculate versors
        for triplet in triplets:
            versors[triplet] = np.concatenate((versors[triplet], compute_versors(triplet).reshape((3,1))),axis=1)   # Each versor contains time evolution, organised as [[x1,x2,...], [y1,y2,...], [z1,z2,...]]
    
    plt.figure(figsize=(8,5))
    
    for triplet in triplets:
        racf = []
        for k in range(200):
            racf.append( autocorrelation(versors[triplet],k) )
        np.savetxt( os.path.join(directory_arrays,f'{i}{triplet}'), racf)

        plt.plot(racf, label=triplet)

    plt.title(f'Mn {i}')
    plt.ylabel(r'$\frac{1}{N-k} \sum_0^{N-k}v_n \cdot v_{n+k}$')
    plt.xlabel('k')
    plt.legend()
    plt.savefig(os.path.join(directory_figures,f'{i}.png'))
    plt.close()

print(f'Execution time: {timer()-timer_start:.1f}s' )