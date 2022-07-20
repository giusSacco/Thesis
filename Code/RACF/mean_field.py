import os, re, sys
import numpy as np
from timeit import default_timer as timer
from argparse import ArgumentParser


class EmptyFileError(Exception):
    def __init__(self, filename : str, msg = None):
        if msg is None:
            msg = f"Warining: File '{filename}' is empty. It will be ignored."
        super().__init__(msg)

timer_start = timer()


program_description = '''Calculates mean magnetic field from files produced by racf.py. Spin is randomised each timestep.'''
PROGNAME = os.path.basename(sys.argv[0])
parser = ArgumentParser(prog = PROGNAME, description = program_description)
parser.add_argument('--dir', dest= 'dir', required=True, help='Results will be saved here. Input file should be in the subdirectory output_data.')
parser.add_argument('--rcut', dest = 'r_cut', required=True, type=float, help = 'cut-off radius (in Angstrom) for magnetic field calculation.')
args_parser = parser.parse_args()
out_dir = args_parser.dir
r_cutoff = args_parser.r_cut
# Read input file
input_dir = os.path.join(out_dir,'output_data')

mu_b = 9.274009994e-24 #J T^-1
g = -2.00231930436182
mu_0 = 1.25663706212e-06
pi = 3.141592653589793
A = mu_0/(4*pi)*g*mu_b*np.sqrt((5/2)*(5/2+1))


def magn_field(r_1,r_2,spin_dir):
    B=0
    if r_2[2] < 9:
        r_2[2] += 165

    n_boxes_z = int((r_cutoff-nv_pos[2])//160)  # Number of boxes to be considered along z (comprehending the #0)
    n_boxes_xy = int((r_cutoff-111.067/2)//111.067) +1  # Number of boxes to be considered along x or y (not counting the #0)
    for x in range(-n_boxes_xy, n_boxes_xy+1):
        for y in range(-n_boxes_xy, n_boxes_xy+1):
            for z in range(n_boxes_z +1):
                pbc_vector = np.array([x,y,z])*np.array([111.067,111.067,160])
                r = r_2 + pbc_vector - r_1
                r_norm = np.linalg.norm(r)
                r_versor = r/r_norm
                if r_norm < r_cutoff:
                    B += 1*(3*r_versor*(np.dot(spin_dir,r_versor)) - spin_dir)/(r_norm**3)
    return B

def random_dir():
    vec = np.random.standard_normal(size=3)
    vec /= np.linalg.norm(vec)
    return vec

B = np.zeros(3)

nv_pos = np.array([111.067/2,111.067/2,-56.6])
j=0

with open(os.path.join(out_dir,f'mean_field_rcut{int(r_cutoff)}.txt'),'w') as out_file:
    
    for filename in os.listdir(input_dir):
        i = 0 #t
        pattern = re.compile(r'''       # Pattern of row of file produced by racf.py
                (\d+\.\d+)\s            # t: float followed by space
                ((?:\d+\.\d+\s){3})     # position: 3 floats separated by space
                (?:\[[^\[]*\]\s){3}   # vectors: anything inside square brackets (three times)
                (\d)                    # flag: integer
                ''',re.VERBOSE)

        # Read input file
        with open(os.path.join(input_dir,filename)) as input_file:
            line = input_file.readline()
            if not line:
                raise EmptyFileError(filename)

            N = int(re.compile(r'N\s=\s(\d+)').search(line).group(1))

            line = input_file.readline() # Skip header
            while line:     # while file is not ended
                if pattern.match(line):
                    groups = pattern.match(line).groups()   # groups are t, positions, flag

                    t = float(groups[0])
                    position = np.array([float(x) for x in groups[1].split()])
                    flag = int(groups[2])
                    alpha, beta, gamma = random_dir()
                    spin_dir = np.array([alpha, beta, gamma])

                    B += magn_field(nv_pos, position, spin_dir)

                elif input_file.readline(): # If pattern is not matched and is not the last row of file
                    print('Warning, line does not match pattern!')

                line = input_file.readline()
                i+=1
            if N != i:
                print(f'Warining: File {filename} is incomplete. Missing rows. {N},{i}')

        if (j+1)%10 == 0 :
            print(f' {j+1}/{len(os.listdir(input_dir))}, {timer()-timer_start:.1f}s')
        j+=1

    out_file.write(f'B_mean = {np.average(B, axis=0)}')
    out_file.write(f'B_squared_mean = {np.average(B**2, axis=0)}')


print(f'Execution time: {timer()-timer_start:.1f}s' )
