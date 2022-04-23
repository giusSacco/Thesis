import os
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from timeit import default_timer as timer


def get_triplets(solvation_shell):
    '''Returns list of all triplets of nearest neighbours water molecules in the solvation shell, corresponding to each plane of the octahedron.
    Should work also if shell is less than 6.'''
    triplets=[]
    for comb in combinations(range(len(solvation_shell)),3):    # For each possible triplet
        pairs = combinations(comb, 2)   # Pairs of atoms in the triplet
        for i,j in pairs:
            # Checks that distances between water molecules is < 4, meaning they are nearest neighbours
            if distances.distance_array(solvation_shell[i].position, solvation_shell[j].position, box=u.dimensions)[0,0] > 4:
                break
        else:
            triplets.append([solvation_shell[index] for index in comb])
    return triplets

def difference_pbc(r1,r2):
    '''Computes difference r1 - r2 accounting for periodic boundaric conditions.'''
    L = u.dimensions[:3]    # u.dimensions is [a, b, c, alpha, beta, gamma]
    return np.remainder(r1 - r2 + L/2., L) - L/2.

def compute_versor(mn, w1, w2, w3):
    '''Returns versor obtained as sum of versors from mn to three waters'''
    v1 = difference_pbc(mn.position, w1.position)
    v2 = difference_pbc(mn.position, w2.position)
    v3 = difference_pbc(mn.position, w3.position)
    vector = v1/np.linalg.norm(v1) + v2/np.linalg.norm(v2) + v3/np.linalg.norm(v3)
    #print(np.linalg.norm(vector / np.linalg.norm(vector)))
    return vector / np.linalg.norm(vector)    # Normalizazion

def shells_are_same(triplets_shell1, triplets_shell2):
    '''Checks if triplets_shell1 == triplets_shell2'''
    # Triplets are casted to set of sets since order is not relevant
    triplets_shell1 = set(frozenset(x) for x in triplets_shell1)
    triplets_shell2 = set(frozenset(x) for x in triplets_shell2)
    return triplets_shell1 == triplets_shell2
     
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
mn = mn_ions[0]
i = 0
j = 0

solvation_shell_before = [water for water in water_molecules if distances.distance_array(mn.position, water.position, box=u.dimensions)[0,0] < 3]  # All OW that are closer than 3 to the Mn
triplets_before = get_triplets(solvation_shell_before)
# Dictionary: versors[mn] has time evolution stored as [[x1,x2,...], [y1,y2,...], [z1,z2,...]]
versors = {mn : np.array(compute_versor(mn, *triplets_before[0])).reshape((3,1)) for mn in mn_ions}
k = 0
for frame in u.trajectory[1:]:
    k+=1
    solvation_shell_now = [water for water in water_molecules if distances.distance_array(mn.position, water.position, box=u.dimensions)[0,0] < 3]
    triplets_now = get_triplets(solvation_shell_now)

    if shells_are_same(triplets_now, triplets_before):  # Update versor
        versors[mn] = np.concatenate((versors[mn], compute_versor(mn, *triplets_before[0]).reshape((3,1))),axis=1)
    else:
        print(k, '1', solvation_shell_before == solvation_shell_now)
        if versors[mn].shape[1] > 100:   # Save only for stable shells
            print(k , 2)
            # Each array is saved as file, columns are x(t), y(t), z(t).
            # i refers to Mn, j distingueshes different files referring to the same Mn
            np.savetxt( os.path.join(directory_arrays, f'{i}_{j}'), versors[mn].reshape((versors[mn].shape[1],3)))
            j += 1
        # When shell changes we reset infos for new shell
        solvation_shell_before = solvation_shell_now
        triplets_before = triplets_now
        versors[mn] = np.array( compute_versor(mn, *triplets_before[0]) ).reshape((3,1))
else:   # at the end of the simulation save anyways
    if versors[mn].shape[1] > 100:
        np.savetxt( os.path.join(directory_arrays, f'{i}_{j}'), versors[mn].reshape((versors[mn].shape[1],3)))

print(f'Execution time: {timer()-timer_start:.1f}s' )