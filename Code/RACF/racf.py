import os, sys
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
            # Checks that distances between water molecules is < threshold_wt_wt, meaning they are nearest neighbours
            if distances.distance_array(solvation_shell[i].position, solvation_shell[j].position, box=u.dimensions)[0,0] > threshold_wt_wt:
                break
        else:
            triplets.append([solvation_shell[index] for index in comb])
    # Check properties of shell
    expected_number_triplets = {6:8, 5:4}   # expected_number_triplets[len(shell)] returns the expected number of planes for a shell of dimension len(shell)
    if len(solvation_shell) not in expected_number_triplets.keys():
        print(f'Warning: {len(solvation_shell)} water in solvation shell of Mn #{mn.id} detected at timestep {k}')
    elif len(triplets) != expected_number_triplets[len(solvation_shell)]:
        print(f'Warning: {expected_number_triplets[len(solvation_shell)]} triplets were expected but {len(triplets)} were found\
             for solvation shell of Mn #{mn.id} detected at timestep {k}')
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

def hist_water_distances(all_distances):
    plt.xlabel('Water distances')
    plt.ylabel('Count')
    plt.yscale('log')
    plt.hist(all_distances, bins=50)
    plt.savefig('hist_water_dist.png')
    plt.close()

def hist_mn_water_distances(all_distances):
    plt.xlabel('Mn-Water distances')
    plt.ylabel('Count')
    plt.yscale('log')
    plt.hist(all_distances, bins=50)
    plt.savefig('hist_mn_wat_dist.png')
    plt.close()

timer_start = timer()

# Verify Python3
if sys.version_info.major != 3:
    print('Error: Program must be run with Python3')
    sys.exit(1)

working_dir = os.path.dirname(__file__)
XTC = os.path.join(working_dir, 'hybr+6MN+36wat.0_200ps.xtc')
TPR = os.path.join(working_dir, 'hybr+6MN+36wat.tpr')
directory_figures = os.path.join(working_dir, 'Figures')
directory_arrays = 'racf_arrays'

#delta_t = 1 # ps

# Parameters
threshold_wt_wt = 3.75  # Solvation shell defined as all OW that are closer than threshold_mn_wt to the Mn
threshold_mn_wt = 2.6   # If the distance between two wt of the same shell is > 2.6 they are assumed to be on opposite sides of the octahedron
plot_histograms = True  # choose if plot histegrams

u = mda.Universe(TPR, XTC)

water_molecules = u.select_atoms('type OW')
mn_ions = u.select_atoms('resname MN')

# Initialization of lists and dictionaries, see below for description
solvation_shell_before = dict(); triplets_before = dict(); all_water_distances = list(); all_dist_mn_wat = list(); versors = dict()
savefile_index = {mn : 0 for mn in mn_ions}  # Used in file saving

for mn in mn_ions:
    # Solvation shell defined as all OW that are closer than threshold_mn_wt to the Mn
    solvation_shell_before[mn] = [water for water in water_molecules if distances.distance_array(mn.position, water.position, box=u.dimensions)[0,0] < threshold_mn_wt]
    # Triplets represent planes of the octahedron of the shell, see get_triplets() for further info
    triplets_before[mn] = get_triplets(solvation_shell_before[mn])

    if plot_histograms:
        # Distance between each pair of water molecules in the shell, will be updated at each timestep.
        all_water_distances.extend( [distances.distance_array(pair[0].position, pair[1].position, box=u.dimensions)[0,0] for pair in combinations(solvation_shell_before[mn],2)] )
        # Distance between mn and water in its shell, , will be updated at each timestep.
        all_dist_mn_wat.extend([distances.distance_array(mn.position, water.position, box=u.dimensions)[0,0] for water in solvation_shell_before[mn]])

    # Identifies rotation of shell, has time evolution stored as [[x1,x2,...], [y1,y2,...], [z1,z2,...]]
    versors[mn] = np.array(compute_versor(mn, *triplets_before[mn][0])).reshape((3,1))  # [mn][0] refers to first triplet of shell of mn

solvation_shell_now = dict(); triplets_now = dict()
k=0 # Keeps track of frame number
for frame in u.trajectory[1:]:
    k+=1

    for i,mn in enumerate(mn_ions): # i indexes mn ions
        # Update info after time evolution
        solvation_shell_now[mn] = [water for water in water_molecules if distances.distance_array(mn.position, water.position, box=u.dimensions)[0,0] < threshold_mn_wt]
        triplets_now[mn] = get_triplets(solvation_shell_now[mn])
        if plot_histograms:
            all_water_distances.extend([distances.distance_array(pair[0].position, pair[1].position, box=u.dimensions)[0,0] for pair in combinations(solvation_shell_now[mn],2)])
            all_dist_mn_wat.extend([distances.distance_array(mn.position, water.position, box=u.dimensions)[0,0] for water in solvation_shell_now[mn]])

        if k == len(u.trajectory) - 1:  # If at the end of simulation save anyways
            if versors[mn].shape[1] > 100:
                np.savetxt( os.path.join(directory_arrays, f'{i}_{savefile_index[mn]}'), versors[mn].transpose())
        elif shells_are_same(triplets_now[mn], triplets_before[mn]):  # Update versor
            versors[mn] = np.concatenate((versors[mn], compute_versor(mn, *triplets_before[mn][0]).reshape((3,1))),axis=1)
        else:   # Shell has changed
            if versors[mn].shape[1] > 100:   # Save only for stable shells
                # Each array is saved as file, columns are x(t), y(t), z(t).
                # i refers to Mn, savefile_index[mn] distingueshes different files referring to the same Mn
                np.savetxt( os.path.join(directory_arrays, f'{i}_{savefile_index[mn]}'), versors[mn].transpose())
                savefile_index[mn] += 1
            # When shell changes we reset infos for new shell
            solvation_shell_before[mn] = solvation_shell_now[mn]
            triplets_before[mn] = triplets_now[mn]
            versors[mn] = np.array( compute_versor(mn, *triplets_before[mn][0]) ).reshape((3,1))

# Plot histograms
if plot_histograms:
    hist_water_distances(all_water_distances)
    hist_mn_water_distances(all_dist_mn_wat)

print(f'Execution time: {timer()-timer_start:.1f}s' )