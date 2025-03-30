import sys
from collections import defaultdict

# Open a log file for writing output
log_file = open('/gpfs/gibbs/pi/gerstein/en325/ENTEX2/AlleleSeq2/test/script_output.log', 'w')

cnames1   = ['chr' + str(c) for c in list(range(1, 23)) + ['X', 'Y', 'M']]
cnames2   = [str(c) for c in list(range(1, 23)) + ['X', 'Y', 'MT']]
cnamesref = []

with open(sys.argv[2], 'r') as cf:
    for line in cf:
        line = line.strip()  # Strip newlines and trailing spaces
        if line.startswith('>'):
            cname = line.split()[0][1:]  # Removed redundant strip()
            if cname in cnames1 or cname in cnames2:
                cnamesref.append(cname)

if len(cnamesref) != 25:
    sys.exit(sys.argv[0] + ': unexpected? chromosome names in ' + sys.argv[2])

vcf_hets_dict = defaultdict(list)

for line in sys.stdin:
    line = line.strip()  # Strip newlines and trailing spaces before processing
    if line.startswith('#CHROM'):
        sample_col_idx = line.split().index(sys.argv[1])
    if not line.startswith('#'):
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = line.split('\t')[:9]
        if FILTER == 'PASS':
            GTl = line.split('\t')[sample_col_idx].split(':')[0].replace('/', '|').split('|')
            al = [REF] + ALT.split(',')
            print("GTl:" + str(GTl), file=log_file)
            if GTl[0].isdigit() and GTl[1].isdigit() and al[int(GTl[0])] != al[int(GTl[1])] and len(al[int(GTl[0])]) == len(al[int(GTl[1])]) == len(REF) == 1:
                print("GTl[0]:" + GTl[0] + "GTl[1]:" + GTl[1], file=log_file)
                
                if CHROM not in cnamesref:
                    if CHROM in cnames2 and cnamesref[0] in cnames1:
                        CHROM = 'chr' + CHROM
                    elif CHROM in cnames1 and cnamesref[0] in cnames2:
                        CHROM = CHROM[3:]
                    elif 'M' in CHROM:
                        for cname in cnamesref:
                            if 'M' in cname:
                                CHROM = cname
                    else:
                        continue  # don't write

                sys.stdout.write('\t'.join([
                    CHROM,
                    str(int(POS)-1),
                    POS,
                    '_'.join([CHROM, POS, REF, al[int(GTl[0])], al[int(GTl[1])]])
                ]) + '\n')

log_file.close()
