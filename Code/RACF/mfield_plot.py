
import re, os
import matplotlib.pyplot as plt
import numpy as np
from timeit import default_timer as timer


timer_start = timer()
N = 20001
def autocorrelation(versor, k):
    #Returns temporal autocorrelation of versor evaluated at k. k is tau/delta_t
    acf = 0
    for i in range(N-k):
        acf += np.dot(versor[i],versor[i+k])
    return acf/(N-k)

# Read input file
input_dir = 'shells_0_20ns.dense'

pattern = re.compile(r'''
        magnetic_field_rcut
        (\d+)
        \.txt
        ''',re.VERBOSE)
t_list = []
Bx_list = list(); By_list= list(); Bz_list = list()
files = [filename for filename in os.listdir(input_dir) if pattern.match(filename)]

for j,file_ in enumerate(files):
    r_cutoff = int(pattern.match(file_).group(1))
    with open(os.path.join(input_dir,file_)) as f:
        line = f.readline(); line = f.readline()
        i =0
        while line:
            t, Bx, By, Bz = line.split(' ')
            t_list.append(float(t))
            Bx_list.append(float(Bx))
            By_list.append(float(By))
            Bz_list.append(float(Bz))
            line = f.readline()
            i+=1
    for x,B in zip([0,1,2], [Bx_list,By_list,Bz_list]):
        plt.figure(x, figsize=(10,10))
        plt.title(f'{B}')
        plt.plot(B,t_list, label=f'r_cut = {r_cutoff}')
        racf = []
        plt.legend()
        plt.figure(x+3)
        plt.title(f'RACF {B}')
        for k in range(500):
            racf.append( autocorrelation(B,k) )
        racf = np.array(racf)/racf[0]
        plt.plot(racf, label=f'r_cut = {r_cutoff}')
        plt.legend()
    print(f'{j+1}/{len(files)}, {timer()-timer_start:.1f}s' )
for x in range(6):
    plt.figure(x)
    plt.savefig(os.path.join(input_dir,f'{x}'))

plt.show()



