import sys
import numpy
from collections import defaultdict

snv_dict = defaultdict(int)

for line in sys.stdin:
    inf = line.strip().split()[-1]  # Added strip() to remove newline characters
    snv_dict[inf] += 1

# Convert dict_values to a list before calculating median
rd_median = numpy.median(list(snv_dict.values()))

sys.stdout.write('#chr\tpos\trd\n')
for snv in snv_dict:
    # Ensure the result of the division is a string for concatenation
    # This will give you a float result, which seems to be what you're expecting
    rd_ratio = snv_dict[snv] / rd_median  # This division results in a float
    sys.stdout.write('\t'.join(snv.split('_')[:2] + [str(rd_ratio)]) + '\n')

sys.stderr.write('median rd:\t' + str(rd_median) + '\n')

