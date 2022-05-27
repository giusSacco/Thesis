import os, re
import numpy as np
import matplotlib.pyplot as plt
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
dir = 'racf_arrays'
z_hist = []
counter_z_frames = 0
all_z = []
for filename in os.listdir(dir):
    t, position, _, _, _, flag = read_file(os.path.join(dir,filename))
    
    for i,pos in enumerate(position):
        all_z.append(pos[2])
        if flag[i] == 1:
            z_hist.append(pos[2])
    counter_z_frames += len([pos[2] for pos in position if 20<pos[2]])
print(f'Execution time: {timer()-timer_start:.1f}s' )
res_time = int(counter_z_frames)/10     # Residence time for 20<z in ns
n_transitions = len([z for z in z_hist if 20<z])
print('Interval analysed is 20<z')
print(f'Residence time is {res_time}ns, Transitions are {n_transitions}') 
print(f'Time per transition is {res_time/n_transitions:.1f}ns')
plt.xlabel('z')
plt.ylabel('External exchanges (flag=1)')
plt.hist(z_hist, log=True, bins = 50)
plt.close()
plt.hist(all_z, log=True, bins = 100)
plt.xlabel('z')
plt.ylabel('z occurrency')
plt.show()