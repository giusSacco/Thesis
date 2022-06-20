import os, re, sys
import numpy as np
from timeit import default_timer as timer
from argparse import ArgumentParser

timer_start = timer()
def read_file(filename):
    pattern = re.compile(r'''       # Pattern of row of file produced by racf.py
            (\d+\.\d+)\s            # t: float followed by space
            ((?:\d+\.\d+\s){3})     # positione: 3 floats separated by space
            ((?:\[[^\[]*\]\s){3})   # vectors: anything inside square brackets (three times)
            (\d)                    # flag: integer
            ''',re.VERBOSE)
    # Initialize lists

    v1 = list(); v2 = list(); v3 = list(); position = list(); t = list(); flag = list()
    # Read input file
    with open(filename) as input_file:
        line = input_file.readline(); line = input_file.readline() # Skip header
        while line:     # while file is not ended
            if pattern.match(line):
                groups = pattern.match(line).groups()   # groups are t, positions, vectors, flag
                # Append data of current row to lists
                t.append(float(groups[0]))
                position.append(np.array([float(x) for x in groups[1].split()]))
                vectors = groups[2].split(r'] [')
                v1.append(np.array([float(x) for x in vectors[0].replace('[','').split()]))
                v2.append(np.array([float(x) for x in vectors[1].split()]))
                v3.append(np.array([float(x) for x in vectors[2].replace(']','').split()]))
                flag.append(int(groups[3]))
            elif input_file.readline(): # If pattern is not matched and is not the last row of file
                print('Warning, line does not match pattern!')

            line = input_file.readline()

    return t, position, v1, v2, v3, flag

program_description = '''Calculates magnetic field from files produced by racf.py.'''
PROGNAME = os.path.basename(sys.argv[0])
parser = ArgumentParser(prog = PROGNAME, description = program_description)
parser.add_argument('--dir', dest= 'dir', required=True, help='Results will be saved here. Input file should be in the subdirectory output_data.')
parser.add_argument('--rcut', dest = 'r_cut', required=True, type=float, help = 'cut-off radious (in Angstrom) for magnetic field calculation.')
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

    for x in [-2,-1,0,1,2]:
        for y in [-2,-1,0,1,2]:
            pbc_vector = np.array([x,y,0])*111.067
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

N = 20001
B = np.zeros((N,3))
nv_pos = np.array([111.067/2,111.067/2,-56.6])
j=0
with open(os.path.join(out_dir,f'magnetic_field_rcut{int(r_cutoff)}.txt'),'w') as out_file:
    out_file.write('t B\n')
    for filename in os.listdir(input_dir):
        t_list, position, v1, v2, v3, flag = read_file(os.path.join(input_dir,filename))
        alpha, beta, gamma = random_dir()
        
        
        for i,t in enumerate(t_list):
            if flag[i] != 0:
                alpha, beta, gamma = random_dir()
            spin_dir = alpha*v1[i] + beta*v2[i] + gamma*v3[i]

            B[i] += magn_field(nv_pos, position[i], spin_dir)
        if (j+1)%10 == 0 :
            print(f' {j+1}/{len(os.listdir(input_dir))}, {timer()-timer_start:.1f}s')
        j+=1
    for i,t in enumerate(t_list):
        out_file.write('{} {} {} {}\n'.format(t, *B[i]))

print(f'Execution time: {timer()-timer_start:.1f}s' )
