
import re, os, sys
import matplotlib.pyplot as plt
import numpy as np
from timeit import default_timer as timer
from argparse import ArgumentParser


timer_start = timer()
def autocorrelation(versor, k):
    #Returns temporal autocorrelation of versor evaluated at k. k is tau/delta_t
    acf = 0
    for i in range(N-k):
        acf += np.dot(versor[i],versor[i+k])
    return acf/(N-k)

# Read input file
PROGNAME = os.path.basename(sys.argv[0])
program_description = '''Plots magnetic field and its autocorrelations.'''
parser = ArgumentParser(prog = PROGNAME, description = program_description)
parser.add_argument('--dir', dest= 'dir', required = True, help='Directory where magnetic field files are.')
parser.add_argument('--kmax', dest= 'k_max', required = True, type = int,help='Autocorrelation is calculated from k=0 to k=k_max where k is the discretised version of tau.')
args_parser = parser.parse_args()
input_dir = args_parser.dir
k_max = args_parser.k_max

pattern = re.compile(r'''
        magnetic_field_rcut
        (\d+)
        \.txt
        ''',re.VERBOSE)

files = [filename for filename in os.listdir(input_dir) if pattern.match(filename)]
Bz2_avg = dict(); By2_avg = dict(); Bx2_avg = dict()
for j,file_ in enumerate(files):
    r_cutoff = int(pattern.match(file_).group(1))
    t_list = []
    Bx_list = list(); By_list= list(); Bz_list = list()
    with open(os.path.join(input_dir,file_)) as f:
        line = f.readline()
        if j == 0:
            N = int(re.compile(r'N\s=\s(\d+)').search(line).group(1))
        line = f.readline()
        i =0
        while line:
            t, Bx, By, Bz = line.split(' ')
            t_list.append(float(t))
            Bx_list.append(float(Bx))
            By_list.append(float(By))
            Bz_list.append(float(Bz))
            line = f.readline()
            i+=1
    Bz_array = np.array(Bz_list); Bx_array = np.array(Bx_list); By_array = np.array(By_list)
    Bz2_avg[r_cutoff] = np.average(Bz_array**2); By2_avg[r_cutoff] = np.average(By_array**2); Bx2_avg[r_cutoff] = np.average(Bx_array**2)

    for x,B,name in zip([0,1,2], [Bx_list,By_list,Bz_list], ['Bx','By','Bz']):
        plt.figure(x, figsize=(10,10))
        plt.title(f'{name}')
        plt.plot(t_list, B, label=f'r_cut = {r_cutoff}')
        racf = []
        plt.legend()
        plt.figure(x+3)
        plt.title(f'RACF {name}')
        for k in range(k_max):
            racf.append( autocorrelation(B,k) )
        racf = np.array(racf)/racf[0]
        plt.plot(racf, label=f'r_cut = {r_cutoff}')
        plt.legend()
    print(f'{j+1}/{len(files)}, {timer()-timer_start:.1f}s' )
for x, name in zip(range(6), ['Bx','By','Bz', 'Bx_corr', 'By_corr', 'Bz_corr']):
    plt.figure(x)
    plt.savefig(os.path.join(input_dir,f'{name}.png'))
    plt.close()

plt.title('<B^2>')
plt.xlabel('r_cutoff')
plt.ylabel('<B^2>')
plt.scatter(Bz2_avg.keys(), Bz2_avg.values(), label = 'z')
plt.scatter(Bx2_avg.keys(), Bx2_avg.values(), label = 'x')
plt.scatter(By2_avg.keys(), By2_avg.values(), label = 'y')
plt.legend()
plt.savefig(os.path.join(input_dir,f'B2avg_vs_cutoff.png'))
plt.close()


print(f'Execution time: {timer()-timer_start:.1f}s' )
