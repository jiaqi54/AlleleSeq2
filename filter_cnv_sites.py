import sys

cnv_dict = {}
# Load CNV data into cnv_dict
with open(sys.argv[2], 'r') as cnvfile:
    for line in cnvfile:
        if not line.startswith('chrm'):
            chrm, snppos, rd = line.strip().split('\t')
            cnv_dict[chrm + '_' + snppos] = rd

# Open a file to log discarded SNP due to CNV filtering
with open(sys.argv[1], 'w') as rm_hetSNV_f:
    rm_hetSNV_f.write('#chr\tref_coord\tcA_cC_cG_cT__rd\n')

    # Process each line from standard input
    for line in sys.stdin:
        if not line.startswith('#'):
            try:
                parts = line.split()
                if len(parts) == 17:  # Ensure the line has the expected number of elements
                    chrm, ref_coord, hap1_coord, hap2_coord, ref_allele, hap1_allele, hap2_allele, cA, cC, cG, cT, cN, ref_allele_ratio, sum_ref_n_alt_cnts, p_binom, warning_hap1, warning_hap2 = parts
                    key = chrm + "_" + ref_coord
                    # Check if the key exists in cnv_dict
                    if key in cnv_dict:
                        rd = float(cnv_dict[key])
                        if 0.5 <= rd <= 1.5:
                            sys.stdout.write(line.strip()+'\t'+str(rd)+'\n')
                        else:
                            rm_hetSNV_f.write('\t'.join([chrm, ref_coord, '_'.join([cA,cC,cG,cT,'',str(rd)])])+'\n')
                    else:
                        # Handle the case where the key is not found
                        print(f"Warning: Key {key} not found in CNV dictionary.", file=sys.stderr)
                else:
                    print(f"Warning: Line skipped due to unexpected format: {line.strip()}", file=sys.stderr)
            except ValueError as e:
                print(f"Error processing line: {line.strip()}. Error: {e}", file=sys.stderr)
        else:
            sys.stdout.write(line.strip()+'\tcnv\n')

#rm_hetSNV_f.close()
