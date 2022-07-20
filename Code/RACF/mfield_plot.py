
import re, os, sys
import matplotlib.pyplot as plt
import numpy as np
from timeit import default_timer as timer
from argparse import ArgumentParser


timer_start = timer()

# Read input file
PROGNAME = os.path.basename(sys.argv[0])
program_description = '''Plots magnetic field and its autocorrelations.'''
parser = ArgumentParser(prog = PROGNAME, description = program_description)
parser.add_argument('--dir', dest= 'dir', required = True, help='Directory where magnetic field files are.')
args_parser = parser.parse_args()
input_dir = args_parser.dir

pattern = re.compile(r'''
        magnetic_field_rcut
        (\d+)
        \.txt
        ''',re.VERBOSE)

files = [filename for filename in os.listdir(input_dir) if pattern.match(filename)]
Bz2_avg = dict(); By2_avg = dict(); Bx2_avg = dict()
Bx_avg = dict(); By_avg = dict(); Bz_avg = dict()
B_mod_avg = dict()
for j,file_ in enumerate(files):
    r_cutoff = int(pattern.match(file_).group(1))
    t_list = []
    Bx_list = list(); By_list= list(); Bz_list = list()
    with open(os.path.join(input_dir,file_)) as f:
        line = f.readline()
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
    Bz_avg[r_cutoff] = np.average(Bz_array); By_avg[r_cutoff] = np.average(By_array); Bx_avg[r_cutoff] = np.average(Bx_array)
    Bz2_avg[r_cutoff] = np.average(Bz_array**2); By2_avg[r_cutoff] = np.average(By_array**2); Bx2_avg[r_cutoff] = np.average(Bx_array**2)
    B_mod_avg[r_cutoff] = np.average(np.sqrt(Bx_array**2 + By_array**2 + Bz_array**2))

plt.title('<B^2>')
plt.xlabel('r_cutoff')
plt.ylabel('<B^2>')
plt.scatter(Bz2_avg.keys(), Bz2_avg.values(), label = 'z')
plt.scatter(Bx2_avg.keys(), Bx2_avg.values(), label = 'x')
plt.scatter(By2_avg.keys(), By2_avg.values(), label = 'y')
plt.legend()
plt.savefig(os.path.join(input_dir,f'B2avg_vs_cutoff.png'))
plt.close()

plt.title('<B>/<|B|>')
plt.xlabel('r_cutoff')
plt.ylabel('<B>/<|B|>')
plt.scatter(Bz_avg.keys(), np.array(list(Bz_avg.values()))/np.array(list(B_mod_avg.values())), label = 'z')
plt.scatter(Bx_avg.keys(), np.array(list(Bx_avg.values()))/np.array(list(B_mod_avg.values())), label = 'x')
plt.scatter(By_avg.keys(), np.array(list(By_avg.values()))/np.array(list(B_mod_avg.values())), label = 'y')
plt.legend()
plt.savefig(os.path.join(input_dir,f'Bavg_vs_cutoff.png'))
plt.close()


print(f'Execution time: {timer()-timer_start:.1f}s' )
