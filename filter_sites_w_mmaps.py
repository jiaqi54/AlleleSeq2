import sys
import read_pileup
import binom

#todo: needs clean-up, also, maybe not carry all columns after filters

def outwrite(l, log_mm='.'):
    sys.stdout.write('\t'.join([l.strip(), log_mm])+'\n')
                                   
def logwrite(l, mm_counts, comment):
    mm_cnts='_'.join([str(x) for x in (mm_counts['A'],mm_counts['C'],mm_counts['G'],mm_counts['T'])])
    log.write('\t'.join(l.split('\t')[:2] + [mm_cnts, comment])+'\n')

def rmsiteswrite(l, mm_counts, comment):
    mm_cnts='_'.join([str(x) for x in (mm_counts['A'],mm_counts['C'],mm_counts['G'],mm_counts['T'])])
    rm_hetSNV_f.write('\t'.join(l.split('\t')[:2] + [mm_cnts+'__'+comment])+'\n')


mode = sys.argv[1]
log = open(sys.argv[3],'w')
if mode != 'adjust': rm_hetSNV_f = open(sys.argv[2], 'w')
if mode != 'adjust': rm_hetSNV_f.write('#chr\tref_coord\tcA_cC_cG_cT__weaker_allele:change_in_ratio\n')
log.write('\t'.join(['#chr', 'pos', '(max_mm):cA_cC_cG_cT', 'mm_hap1_warn;mm_hap2_warn;mm_log\n']))

mm_pileup_dict = read_pileup.pileup_to_basecnts(sys.argv[4:])

sys.stdout.write(sys.stdin.readline().strip()+'\tmmap_log\n')

for line in sys.stdin:
    (chrm, ref_coord, hap1_coord, hap2_coord, ref_allele, hap1_allele, hap2_allele, 
            cA, cC, cG, cT, cN, ref_allele_ratio, sum_ref_n_alt_cnts, 
            p_binom, cnv) = line.split('\t')
  
    hap1_mm_basecnts = mm_pileup_dict.get(hap1_coord,{'A':0, 'C':0, 'G':0, 'T':0, 'N':0, 'warning':'.'})
    hap2_mm_basecnts = mm_pileup_dict.get(hap2_coord,{'A':0, 'C':0, 'G':0, 'T':0, 'N':0, 'warning':'.'})

    # again because of the possibility of reads flanking, say, a misphased indel, will check *both* alleles for multi-mapping reads bearing each allele
    mm_basecnts = {
            'A': max(hap1_mm_basecnts['A'], hap2_mm_basecnts['A']),
            'C': max(hap1_mm_basecnts['C'], hap2_mm_basecnts['C']),
            'G': max(hap1_mm_basecnts['G'], hap2_mm_basecnts['G']),
            'T': max(hap1_mm_basecnts['T'], hap2_mm_basecnts['T']),
            'N': max(hap1_mm_basecnts['N'], hap2_mm_basecnts['N'])
    }

    # now, only bother if the unique read counts aren't equal and the 'weaker' allele is seen in multi-mapped reads
    um_basecnts = {'A':float(cA), 'C':float(cC), 'G':float(cG), 'T':float(cT), 'N':float(cN)}
    if   um_basecnts[hap1_allele] > um_basecnts[hap2_allele]: 
        weaker = hap2_allele
    elif um_basecnts[hap1_allele] < um_basecnts[hap2_allele]: 
        weaker = hap1_allele
    else: 
        outwrite(l=line)
        logwrite(line, mm_basecnts, ';'.join([hap1_mm_basecnts['warning'], hap2_mm_basecnts['warning'], 'eq_uniq_cnts']))
        continue



    if mm_basecnts[weaker] == 0: 
        outwrite(l=line, log_mm=weaker+':'+str(mm_basecnts[weaker]))
        logwrite(line, mm_basecnts, ';'.join([hap1_mm_basecnts['warning'], hap2_mm_basecnts['warning'], weaker+':'+str(mm_basecnts[weaker])]))
    else: 
        # now, either remove the site if the allelic ratio changes by more than a threshold
        if mode != 'adjust':
            new = (um_basecnts[weaker] + float(mm_basecnts[weaker])) / (um_basecnts[hap1_allele] + um_basecnts[hap2_allele] + mm_basecnts[weaker]) 
            old = um_basecnts[weaker] / (um_basecnts[hap1_allele] + um_basecnts[hap2_allele])
            if (new - old) > float(mode):
                logwrite(line, mm_basecnts, ';'.join([hap1_mm_basecnts['warning'], hap2_mm_basecnts['warning'], weaker+'_removed']))
                rmsiteswrite(line, mm_basecnts, weaker + ':' + str(round(new-old,2)))
            else:
                logwrite(line, mm_basecnts, ';'.join([hap1_mm_basecnts['warning'], hap2_mm_basecnts['warning'], weaker+'_within_thresh']))
                outwrite(line, weaker+':'+str(mm_basecnts[weaker])+';within_thresh')
                

        # or adjust counts in the most conservative way: add the mm counts to the weaker allele:
        # all of them or until balanced with the stonger
        # to make sure the count imbalance is not caused by the multimapping reads              
        else:
            adj = min((um_basecnts[weaker] + mm_basecnts[weaker]), max(um_basecnts[hap1_allele], um_basecnts[hap2_allele])) 
            diff = adj - um_basecnts[weaker] 
            um_basecnts[weaker] = adj
            
            new_tot = int(um_basecnts[hap1_allele] + um_basecnts[hap2_allele])
            new_ratio = float(um_basecnts[ref_allele])/float(new_tot)
            new_p_binom = binom.binomtest(um_basecnts[ref_allele], new_tot, 0.5)

            sys.stdout.write('\t'.join([
                chrm, ref_coord, hap1_coord, hap2_coord, ref_allele, hap1_allele, hap2_allele,
                str(int(um_basecnts['A'])),
                str(int(um_basecnts['C'])),
                str(int(um_basecnts['G'])),
                str(int(um_basecnts['T'])),
                str(int(um_basecnts['N'])),
                str(new_ratio), str(new_tot), str(new_p_binom), 
                cnv[:-1], weaker+':+'+str(int(diff))
                ])+'\n')

            logwrite(line, mm_basecnts, ';'.join([hap1_mm_basecnts['warning'], hap2_mm_basecnts['warning'], weaker+':+'+str(int(diff))]))


if mode != 'adjust': rm_hetSNV_f.close()
log.close()
