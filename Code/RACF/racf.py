import os, sys
from argparse import ArgumentParser
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
    if len(solvation_shell) == 4:
        pass
    elif len(solvation_shell) not in expected_number_triplets.keys():
        print(f'Warning: {len(solvation_shell)} water in solvation shell of Mn #{mn.id} detected at timestep {k}')
    elif len(triplets) != expected_number_triplets[len(solvation_shell)]  and minor_warnings: # Print this warning only if --mwarn is expressed in command line
        print(f'Warning: {expected_number_triplets[len(solvation_shell)]} triplets were expected but {len(triplets)} were found \
        for solvation shell of Mn #{mn.id} detected at timestep {k}')
    return triplets

def difference_pbc(r1,r2):
    '''Computes difference r1 - r2 accounting for periodic boundaric conditions.'''
    L = u.dimensions[:3]    # u.dimensions is [a, b, c, alpha, beta, gamma]
    return np.remainder(r2 - r1 + L/2., L) - L/2.

def compute_versors(mn, w1, w2, w3):
    '''Returns versors that point from the mn to the waters'''
    v1 = difference_pbc(mn.position, w1.position)
    v2 = difference_pbc(mn.position, w2.position)
    v3 = difference_pbc(mn.position, w3.position)
    return v1/np.linalg.norm(v1), v2/np.linalg.norm(v2), v3/np.linalg.norm(v3)    # Normalizazion

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
    plt.savefig(os.path.join(directory_output,'hist_water_dist.png'))
    plt.close()

def hist_mn_water_distances(all_distances):
    plt.xlabel('Mn-Water distances')
    plt.ylabel('Count')
    plt.yscale('log')
    plt.hist(all_distances, bins=50)
    plt.savefig(os.path.join(directory_output,'hist_mn_wat_dist.png'))
    plt.close()

program_description = '''Analysis of dynamics of Mn ions and their solvation shell.'''

timer_start = timer()

# Verify Python3
if sys.version_info.major != 3:
    print('Error: Program must be run with Python3')
    sys.exit(1)

working_dir = os.path.dirname(__file__)


# Parsing command-line inputs
PROGNAME = os.path.basename(sys.argv[0])
parser = ArgumentParser(prog = PROGNAME, description = program_description)
parser.add_argument('--xtc', dest= 'XTC', required = True, help='XTC file.')
parser.add_argument('--tpr', dest= 'TPR', required = True, help='TPR file.')
parser.add_argument('--startingindex', type = int, required = True, dest= 'starting_index', help='Analyzes trajectory of 10 ions starting from specified index.')
parser.add_argument('-n', dest= 'N', type = int, help='Number of frames to be analysed from the last. Default is all trajectory.')
parser.add_argument('--plot',action='store_true', dest= 'plot_histograms', help='Produce histograms of distances, defaults to False')
parser.add_argument('--mwarn', action='store_true', dest= 'minor_warnings', help='Print minor warning such as different number of triplets than expected. Defaults to False.')
args_parser = parser.parse_args()
XTC = args_parser.XTC
TPR = args_parser.TPR
Ion_index = args_parser.starting_index
plot_histograms = args_parser.plot_histograms
minor_warnings = args_parser.minor_warnings
N = args_parser.N

# Create output directory
directory_output = XTC[:-4]
if not os.path.exists(directory_output):
    os.mkdir(directory_output)
if not os.path.exists(os.path.join(directory_output,'output_data')):
    os.mkdir(os.path.join(directory_output,'output_data'))

# Parameters
threshold_wt_wt = 3.75  # Solvation shell defined as all OW that are closer than threshold_mn_wt to the Mn
threshold_mn_wt = 2.6   # If the distance between two wt of the same shell is > 2.6 they are assumed to be on opposite sides of the octahedron
plot_histograms = True  # choose if plot histegrams

# Initialize universe
u = mda.Universe(TPR, XTC)

if N is None:
    N = len(u.trajectory)

water_molecules = u.select_atoms('type OW')
mn_ions = u.select_atoms('resname MN')


# Initialization of lists and dictionaries, see below for description
all_water_distances = list(); all_dist_mn_wat = list()

for i,mn in enumerate(mn_ions[Ion_index:Ion_index+10], start = Ion_index): # i indexes mn ions
    # Solvation shell defined as all OW that are closer than threshold_mn_wt to the Mn
    solvation_shell_before = [water for water in water_molecules if distances.distance_array(mn.position, water.position, box=u.dimensions)[0,0] < threshold_mn_wt]
    # Triplets represent planes of the octahedron of the shell, see get_triplets() for further info
    triplets_before = get_triplets(solvation_shell_before)
    
    with open(os.path.join(directory_output,'output_data',f'{i}'),'w') as out_file:
        out_file.write(' '.join(('t','x','y','z','v1','v2','v3','flag',f'N = {N}','\n')))    # Header
        try:
            k=0 # Keeps track of frame number
            for frame in u.trajectory [-N:]:
                k+=1
                # Print execution progress
                if (k*10) % N == 0:
                    print(f'{i+1}/{len(mn_ions)}, {k/N*100:.0f}%')

                # Update info after time evolution       
                solvation_shell_now = [water for water in water_molecules if distances.distance_array(mn.position, water.position, box=u.dimensions)[0,0] < threshold_mn_wt]
                triplets_now = get_triplets(solvation_shell_now)
                if len(triplets_now) == 0:  # If no triplets found, save last versors
                    out_file.write(' '.join(str(x) for x in [u.trajectory.time,*mn.position,*versors,3,'\n']))
                    continue

                flag = 0    # Shell hasn't changed
                if not shells_are_same(triplets_now, triplets_before):     
                    if solvation_shell_before != solvation_shell_now:  # External exchange
                        flag = 1
                    else:   # Internal Exchange
                        flag = 2
                    # When shell changes we reset infos for new shell
                    solvation_shell_before = solvation_shell_now
                    triplets_before = triplets_now

                versors = compute_versors(mn, *triplets_before[0])            
                out_file.write(' '.join(str(x) for x in [u.trajectory.time,*mn.position,*versors,flag,'\n']))

                if plot_histograms:
                    all_water_distances.extend([distances.distance_array(pair[0].position, pair[1].position, box=u.dimensions)[0,0] for pair in combinations(solvation_shell_now,2)])
                    all_dist_mn_wat.extend([distances.distance_array(mn.position, water.position, box=u.dimensions)[0,0] for water in solvation_shell_now])
        except Exception as err:    # If an error happens during the trajectory analysis, pass to next ion 
            print(f'Warning!: During the analysis of {mn.id} at frame #{k} the following exception occurred: {err} \n Proceeding to the next ion.\n')
            out_file.write(f'Warning! The following exception occurred: {err}')
            continue
    print(f'Partial execution time: {timer()-timer_start:.1f}s' )

# Plot histograms
if plot_histograms:
    hist_water_distances(all_water_distances)
    hist_mn_water_distances(all_dist_mn_wat)

print(f'Total execution time: {timer()-timer_start:.1f}s' )
