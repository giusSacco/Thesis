import os, re
import numpy as np
from scipy.constants import mu_0, pi, hbar
import scipy.constants
from timeit import default_timer as timer

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

# Read input file
dir = 'racf_arrays'; filename = '0'
t_list, position, v1, v2, v3, flag = read_file(os.path.join(dir,filename))

mu_b = scipy.constants.physical_constants['Bohr magneton'][0]  # 9.274009994e-24 J T^-1
g = scipy.constants.physical_constants['electron g factor'][0]  # -2.00231930436182
A = mu_0/(4*pi)*g*mu_b*np.sqrt((5/2)*(5/2+1))
def magn_field(r_1,r_2,spin_dir):
    r = r_2 - r_1
    r_versor = r/np.linalg.norm(r)
    return A*(3*r_versor*(np.dot(spin_dir,r_versor)) - spin_dir)/(np.linalg.norm(r)**3)

def random_dir():
    vec = np.random.standard_normal(size=3)
    vec /= np.linalg.norm(vec)
    return vec
with open(f'magnetic_field.txt','w') as out_file:
    out_file.write('t B')
    for t in t_list:
        pass




print(f'Execution time: {timer()-timer_start:.1f}s' )
