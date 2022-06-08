import os, re
import numpy as np
from scipy.constants import mu_0, pi, hbar
import scipy.constants
from timeit import default_timer as timer
import matplotlib.pyplot as plt
from scipy.fft import fft

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
out_dir = 'filenuovo'
input_dir = os.path.join(out_dir,'output_data')

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

N = 1001
B = np.zeros((N,3))
nv_pos = np.array([111.067/2,111.067/2,1.041])

with open(os.path.join(out_dir,f'magnetic_field.txt'),'w') as out_file:
    out_file.write('t B\n')

    for filename in os.listdir(input_dir):
        t_list, position, v1, v2, v3, flag = read_file(os.path.join(input_dir,filename))
        alpha, beta, gamma = random_dir()
        
        for i,t in enumerate(t_list):
            if flag != 0:
                alpha, beta, gamma = random_dir()
            spin_dir = alpha*v1[i] + beta*v2[i] + gamma*v3[i]
            #out_file.write('{} {} {}\n'.format(*magn_field(nv_pos, position[i], spin_dir)))
            B[i] += magn_field(nv_pos, position[i], spin_dir)

    for i,t in enumerate(t_list):
        out_file.write('{} {} {} {}\n'.format(t, *B[i]))

def autocorrelation(versor, k):
    '''Returns temporal autocorrelation of versor evaluated at k. k is tau/delta_t'''
    acf = 0
    for i in range(N-k):
        acf += np.dot(versor[i],versor[i+k])
    return acf/(N-k)

'''plt.plot(B[:,0])
plt.figure()
plt.plot(B[:,1])
plt.figure()
plt.plot(B[:,2])

for x in [0,1,2]:
    racf = []
    plt.figure()
    for k in range(N//2):
        racf.append( autocorrelation(B[:,x],k) )
    plt.plot(racf)
plt.show()'''

print(f'Execution time: {timer()-timer_start:.1f}s' )
