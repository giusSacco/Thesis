import os, re, sys
from argparse import ArgumentParser
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

PROGNAME = os.path.basename(sys.argv[0])
program_description ='''Creates heatmap of <B^2>'''
parser = ArgumentParser(prog = PROGNAME, description = program_description)
parser.add_argument('--dir', dest= 'dir',default = '.',  help='Directory where B files are located.')
parser.add_argument('-n', dest = 'N', nargs=2, required=True, type=int, help = 'Range of frames to be analysed')
parser.add_argument('-ng', dest = 'ng', required=True, type=int, help = 'n of the grid')
parser.add_argument('--rcut', dest = 'r_cut', required=True, type=int, help = 'cut-off radius (in Angstrom) for magnetic field calculation.')
args_parser = parser.parse_args()
n_grid = args_parser.ng
r_cutoff = args_parser.r_cut
N_start, N_end = args_parser.N
pattern = f'''
        B\_rcut
        {r_cutoff}
        \_n{N_start}\_{N_end}
        \_nxy
        {n_grid}
        (\d)(\d)
        \.txt
        '''
file_pattern = re.compile(pattern,re.VERBOSE)

input_dir = args_parser.dir
files = [filename for filename in os.listdir(input_dir) if file_pattern.match(filename)]
if len(files) == 0:
    print(f'ERROR: No files are found in {input_dir} that match the pattern.')
    sys.exit(1)
if len(files) != n_grid**2:
    print(f'ERROR: {len(files)} files are found in {input_dir} while {n_grid**2} are expected.')
    sys.exit(1)

print(files, len(files))

Bz_2_avg = np.zeros((n_grid,n_grid))

def read_bz2(filename):
    with open(os.path.join(input_dir,filename)) as f:
        f.readline()
        line = f.readline()
        bz2 = float(line.split(' ')[-1])
        return bz2

for nx in range(n_grid):
    for ny in range(n_grid):
        Bz_2_avg[nx,ny] = read_bz2(f'''B_rcut{r_cutoff}_n{N_start}_{N_end}_nxy{n_grid}{nx}{ny}.txt''')
print(Bz_2_avg)
# 3. Plot the heatmap
plt.figure(figsize=(10,10))
heat_map = sns.heatmap( Bz_2_avg, linewidth = 1 , annot = True)
plt.title( "HeatMap <B^2>" )
plt.savefig(os.path.join(input_dir,f'heatmap_rcut{r_cutoff}_n{N_start}_{N_end}_nxy{n_grid}.png'))
plt.gca().invert_yaxis().show()
