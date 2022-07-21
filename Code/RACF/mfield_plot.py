
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
parser.add_argument('--dir', dest= 'dir', required = True, help='Directory where to look for "B_files" folder.')
args_parser = parser.parse_args()
input_dir = os.path.join(args_parser.dir,"B_files")

pattern = re.compile(r'''
        B\_rcut
        (\d+)
        \.txt
        ''',re.VERBOSE)

files = [filename for filename in os.listdir(input_dir) if pattern.match(filename)]
if len(files) == 0:
    print(f'Error: no magnetic field files were found in {input_dir}')
    sys.exit()

B_avg = np.zeros( (len(files), 3) )
B2_avg = np.zeros( (len(files), 3) )
Bmod_avg = np.zeros(len(files))
r_cutoffs = np.zeros(len(files))

for j,file_ in enumerate(files):
    r_cutoffs[j] = int(pattern.match(file_).group(1))
    B_ = np.loadtxt(os.path.join(input_dir,file_))
    B_avg[j] = B_[0]
    B2_avg[j] = B_[1]
    Bmod_avg[j] = B_[2,0] +1

plt.title('<B^2>')
plt.xlabel('r_cutoff')
plt.ylabel('<B^2>')
plt.scatter(r_cutoffs, B2_avg[:,0], label = 'x')
plt.scatter(r_cutoffs, B2_avg[:,1], label = 'y')
plt.scatter(r_cutoffs, B2_avg[:,2], label = 'z')
plt.legend()
plt.savefig(os.path.join(args_parser.dir,f'B2avg_vs_cutoff.png'))
plt.close()

plt.title('<B>/<|B|>')
plt.xlabel('r_cutoff')
plt.ylabel('<B>/<|B|>')
plt.scatter(r_cutoffs, B_avg[:,0]/Bmod_avg, label = 'x')
plt.scatter(r_cutoffs, B_avg[:,1]/Bmod_avg, label = 'y')
plt.scatter(r_cutoffs, B_avg[:,2]/Bmod_avg, label = 'z')
plt.legend()
plt.savefig(os.path.join(args_parser.dir,f'Bavg_vs_cutoff.png'))
plt.close()


print(f'Execution time: {timer()-timer_start:.1f}s' )
