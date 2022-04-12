import re

header_pattern = re.compile(r'''
.*
t(\d+
\.
\d+)
''',re.VERBOSE)
numeric_row_pattern = re.compile(r'''((?:
\d+
[\n\r\s]+
)+)
''', re.VERBOSE)
shell = {}
fname = 'MNshells_test.ndx'
times = []
with open(fname) as f:
    i = 0
    line = f.readline()
    while line:
        if header_pattern.match(line):
            t = header_pattern.match(line).group(1)
            times.append(t)
            shell[t]=set()
            line = f.readline()
            i +=1
            while numeric_row_pattern.match(line):
                row_numbers = {int(number) for number in numeric_row_pattern.match(line).group(1).split()}
                shell[t].update(row_numbers)
                line = f.readline()
            else:
                line = f.readline()

        else:
            line=f.readline()

print(shell[times[0]] ^ shell[times[0]])
print(len(shell[times[0]]))
                
