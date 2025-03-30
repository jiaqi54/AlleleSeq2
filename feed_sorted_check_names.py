import sys
import os
from collections import defaultdict

# Prepare expected chromosome name conventions.
cnames1 = ['chr' + str(c) for c in list(range(1, 23)) + ['X', 'Y', 'M']]
cnames2 = [str(c) for c in list(range(1, 23)) + ['X', 'Y', 'MT']]
cnamesbam = []

# Read chromosome names from a file.
with open(sys.argv[1], 'r') as cf:
    for line in cf:
        cname = line.split()[0]
        if cname in cnames1 or cname in cnames2:
            cnamesbam.append(cname)

# Check for expected number of chromosome names.
if len(cnamesbam) != 25:
    sys.exit(f"{sys.argv[0]}: unexpected chromosome names in {sys.argv[1]}")

hets_dict = defaultdict(list)

# Process each line of input.
try:
    for line in sys.stdin:
        parts = line.strip().split('\t')
        if len(parts) < 4:
            continue  # Skip lines that don't have exactly 4 parts.
        
        chrm, c2, c3, c4 = parts
        
        # Adjust chromosome names based on the conventions detected in the BAM file's chromosome names.
        if chrm not in cnamesbam:
            if chrm in cnames2 and cnamesbam[0] in cnames1:
                chrm = 'chr' + chrm
            elif chrm in cnames1 and cnamesbam[0] in cnames2:
                chrm = chrm[3:]
            elif 'M' in chrm:
                for cname in cnamesbam:
                    if 'M' in cname:
                        chrm = cname
            else:
                sys.exit(f"{sys.argv[0]}: unexpected chromosome names in stdin")
        
        # Append the processed line to the dictionary.
        hets_dict[chrm].append('\t'.join([chrm, c2, c3, c4]))

    # Output the processed lines.
    for c in cnamesbam:
        for l in hets_dict[c]:
            sys.stdout.write(l + '\n')

except BrokenPipeError:
    # Close STDERR after a BrokenPipeError to avoid traceback.
    sys.stderr.close()
    os._exit(0)  # Exit without throwing an error.
