import re

fname = 'MNshells_test.ndx'

header_pattern = re.compile(r'''    # regular expression which recognises header and captures t
.*      # Any characters
t(\d+   # "t" followed by a number, which is captured
\.
\d+)
''',re.VERBOSE)

numeric_row_pattern = re.compile(r'''((?:   # regular expression of row of integers
\d+         # integer
[\n\r\s]+   # followed by whitespace or \n
)+)         # repeated one or more times
''', re.VERBOSE)

shell = {}      # Dictionary: shell[time] is set of indices of atoms saved at a given time
times = []      # List of times in headers

with open(fname) as f:
    line = f.readline()
    while line:
        if header_pattern.match(line):
            t = header_pattern.match(line).group(1) # time
            times.append(t)     # update list of times
            shell[t]=set()      # new set corresconding to t
            line = f.readline()
            while numeric_row_pattern.match(line):
                row_numbers = {int(number) for number in numeric_row_pattern.match(line).group(1).split()}  # save integers in row
                shell[t].update(row_numbers)       # add row to shell
                line = f.readline()
            else:
                line = f.readline()
        else:
            line=f.readline()

# Plots
'''n_transition_per_frame = []
for i in range(len(shell)-1):
    n_transition_per_frame.append( len(shell[times[i]] ^ shell[times[i+1]])//3)
print(n_transition_per_frame)
import matplotlib.pyplot as plt
plt.hist(n_transition_per_frame, bins=[-.5+i for i in range(max(n_transition_per_frame) +1)])   # istogramma della frequenza di n_transition_per_frame
plt.figure()
plt.plot([len(shell[t]) for t in times])        # numero di atomi salvati per frame in funzione del tempo
plt.show()'''

# Create output file with all indices
final_set = set().union(*shell.values())    # Union of all sets

output_file = 'all_indices.ndx'
with open(output_file, 'w') as f:
    f.write('[ same_residue_as_(resname_WAT_and_within_0.25_of_name_MN)_f0_t000000.000 ]\n') # Header
    indices_sorted = sorted(final_set)      # Order indices in increasing order
    indices_sorted = [str(x) for x in indices_sorted]      # Cast to string

    list_of_rows = [indices_sorted[i:i+15] for i in range(0, len(indices_sorted), 15)]
    for row in list_of_rows:
        row.append('\n')    # Add newline
        f.write(' '.join(row))