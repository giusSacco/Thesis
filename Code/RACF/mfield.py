import os, re, sys
import numpy as np
from timeit import default_timer as timer
from argparse import ArgumentParser
import MDAnalysis as mda


class EmptyFileError(Exception):
    def __init__(self, filename : str, msg = None):
        if msg is None:
            msg = f"Warining: File '{filename}' is empty. It will be ignored."
        super().__init__(msg)

timer_start = timer()

program_description = '''Calculates magnetic field from files produced by racf.py.'''
PROGNAME = os.path.basename(sys.argv[0])
parser = ArgumentParser(prog = PROGNAME, description = program_description)
parser.add_argument('--dir', dest= 'dir', required=True, help='Results will be saved here. Input file should be in the subdirectory output_data.')
parser.add_argument('--rcut', dest = 'r_cut', required=True, type=float, help = 'cut-off radious (in Angstrom) for magnetic field calculation.')
args_parser = parser.parse_args()
dir = args_parser.dir
r_cutoff = args_parser.r_cut


mu_b = 9.274009994e-24 #J T^-1
g = -2.00231930436182
mu_0 = 1.25663706212e-06
pi = 3.141592653589793
A = mu_0/(4*pi)*g*mu_b*np.sqrt((5/2)*(5/2+1))


def magn_field(r_1,r_2):
    B=0

    if r_2[2] < 9:
        r_2[2] += 165

    n_boxes_z = int((r_cutoff-nv_pos[2])//165)  # Number of boxes to be considered along z (comprehending the #0)
    n_boxes_xy = int((r_cutoff-111.067/2)//111.067) +1  # Number of boxes to be considered along x or y (not counting the #0)
    for x in range(-n_boxes_xy, n_boxes_xy+1):
        for y in range(-n_boxes_xy, n_boxes_xy+1):
            for z in range(n_boxes_z +1):
                pbc_vector = np.array([x,y,z])*np.array([111.067,111.067,165])
                r = r_2 + pbc_vector - r_1
                r_norm = np.linalg.norm(r)
                r_versor = r/r_norm
                if r_norm < r_cutoff:
                    spin_dir = np.array([*random_dir()])
                    B += 1*(3*r_versor*(np.dot(spin_dir,r_versor)) - spin_dir)/(r_norm**3)
    return B

def random_dir():
    vec = np.random.standard_normal(size=3)
    vec /= np.linalg.norm(vec)
    return vec


XTC = [file_ for file_ in os.listdir(dir) if file_.endswith('.xtc')]
TPR = [file_ for file_ in os.listdir(dir) if file_.endswith('.tpr')]
if len(XTC)!=1 or len(TPR)!=1:
    if len (XTC) == 0:
        print(f'Error: no xtc found in {dir}')
    if len (XTC) > 1:
        print(f'Error: more than one xtc found in {dir}')
    if len (XTC) == 0:
        print(f'Error: no tpr found in {dir}')
    if len (XTC) > 1:
        print(f'Error: more than tpr xtc found in {dir}')
    sys.exit(1)

u = mda.Universe(TPR[0], XTC[0])
N = len(u.trajectory)//5
B = np.zeros((N,3))
mn_ions = u.select_atoms('resname MN')

nv_pos = np.array([111.067/2,111.067/2,-56.6])
j=0
for frame in u.trajectory[-N:]:
    for mn in mn_ions:
        B[j] += magn_field(nv_pos, mn.position)
    if (j*10) % N == 0:
        print(f' {j/N*100:.1f}%, {timer()-timer_start:.1f}s')
    j+=1

np.savetxt(os.path.join(dir,f'magnetic_field_rcut{int(r_cutoff)}.txt'), B)

print(f'Execution time: {timer()-timer_start:.1f}s' )
