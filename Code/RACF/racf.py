import os
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda

def get_triplets(solvation_shell):
    triplets=[]
    for comb in combinations(range(len(solvation_shell)),3):
        pairs = combinations(comb, 2)
        for i,j in pairs:
            if np.linalg.norm(solvation_shell[i].position - solvation_shell[j].position) > 4:
                break
        else:
            triplets.append(comb)
            #print(comb)
    #print(len(triplets))
    return triplets

def compute_versors(triplet):
    a = solvation_shell[triplet[0]].position - solvation_shell[triplet[1]].position
    b = solvation_shell[triplet[1]].position - solvation_shell[triplet[2]].position
    #c = solvation_shell[triplet[2]].position - solvation_shell[triplet[0]].position

    vector_1 = np.cross(a,b)
    #vector_2 = np.cross(a,c)
    vector_1 /= np.linalg.norm(vector_1)
    #vector_2 /= np.linalg.norm(vector_2)

    return vector_1

def autocorrelation(versor, k):
    N = versor.shape[1]

    acf = 0
    for i in range(N-k):
        acf += np.dot(versor[:,i],versor[:,i+k])
    return acf/(N-k)

working_dir = os.path.dirname(__file__)
XTC = os.path.join(working_dir, 'hybr+6MN+36wat.0_200ps.xtc')
TPR = os.path.join(working_dir, 'hybr+6MN+36wat.tpr')

u = mda.Universe(TPR, XTC)
#print(u.atoms.residues.resnames)

water_molecules = u.select_atoms('type OW')
mn_ions = u.select_atoms('resname MN')

for i,mn in enumerate(mn_ions):
    # Get solvation shell
    solvation_shell = [water for water in water_molecules if np.linalg.norm(mn.position - water.position) < 3]

    # Get triplets

    triplets = get_triplets(solvation_shell)
    versors = {triplet : np.array(compute_versors(triplet)).reshape((3,1)) for triplet in triplets}

    for frame in u.trajectory[1:]:
        # Calculate versors
        for triplet in triplets:
            versors[triplet] = np.concatenate((versors[triplet], compute_versors(triplet).reshape((3,1))),axis=1)
    plt.figure(figsize=(8,5))
    for triplet in triplets:
        racf = []
        for k in range(200):
            racf.append( autocorrelation(versors[triplet],k) )
        plt.plot(racf, label=triplet)
    plt.title(f'Mn {i}')
    plt.ylabel(r'$\frac{1}{N-k} \sum_0^{N-k}v_n \cdot v_{n+k}$')
    plt.xlabel('k')
    plt.legend()
    plt.savefig(os.path.join('Figures',f'{i}.png'))
    plt.close()
