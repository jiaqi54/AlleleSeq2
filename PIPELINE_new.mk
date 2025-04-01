#######################
### PIPELINE PARAMS ###
#######################


### system / executables ##

PL                               := /vast/palmer/pi/gerstein/jl3786/AlleleSeq2-master2
NTHR                             := 8 # multithread, works for mapping, sorting, fastqc
SAMTOOLS                         := /vast/palmer/apps/avx2/software/SAMtools/1.16-GCCcore-10.2.0/bin/samtools
PICARD                           := /vast/palmer/apps/avx2/software/picard/2.25.6-Java-11/picard.jar
JAVA                             := java
JAVA_MEM                         := 50g
STAR                             := /vast/palmer/apps/avx2/software/STAR/2.7.7a-GCCcore-10.2.0/bin/STAR
FASTQC                           := /vast/palmer/apps/avx2/software/FastQC/0.11.9-Java-11/fastqc
CUTADAPT                         := /vast/palmer/apps/avx2/software/cutadapt/3.4-GCCcore-10.2.0/bin/cutadapt

### input files / paths ##

READS_R1                         :=
READS_R2                         :=

PGENOME_DIR                      := NULL
VCF_SAMPLE_ID                    := NULL


### params ##

ALIGNMENT_MODE                           := NULL # can be 'ASE', 'ASB', 'custom', 'ASCA' -- currently, for with known adapters (if present) only
RM_DUPLICATE_READS                       := on  # with 'on' duplicate reads will be removed using picard
PERFORM_FASTQC                           := on

# needed for all: ASE, ASB, custom, or ASCA:
GenomeIdx_STAR_diploid                   := $(PGENOME_DIR)/STAR_idx_diploid
STAR_outFilterMismatchNoverReadLmax      := 0.03
STAR_outFilterMatchNminOverLread         := 0.95
# STAR_readFilesCommand                    := zcat # zcat, cat, etc
STAR_readFilesCommand := $(shell \
    if echo $(READS_R1) | grep -q '\.gz$$' || echo $(READS_R2) | grep -q '\.gz$$'; then \
        echo zcat; \
    else \
        echo cat; \
    fi)
STAR_limitSjdbInsertNsj                  := 1500000 # star default is 1000000

# needed if ASE
REFGENOME_VERSION             := GRCh37   #GRCh38 or CRCh37
Annotation_diploid            := $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_diploid.gencode.v19.annotation.gtf

ifeq ($(REFGENOME_VERSION), GRCh38)
  Annotation_diploid          := $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_diploid.gencode.v24.annotation.gtf

else ifeq ($(REFGENOME_VERSION), chm13)
  Annotation_diploid          := $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_diploid.chm13v2.0_RefSeq_Liftoff_v5.1.gff3

endif

STAR_sjdbOverhang                        := 100  #STAR default will work almost as well as with the ideal (readlength -1) value according to the manual

# if custom alignment mode
STAR_parameters_file                     := $(PL)/STAR_custom_parameters_sample_file

# needed if ASCA, atac-seq:
R1_ADAPTER_SEQ                           := CTGTCTCTTATA
R2_ADAPTER_SEQ                           := CTGTCTCTTATA
# for ASCA, Nextera and transposase adapter sequences trimming, similar to what ENCODE prototype peak calling pipeline finds in ENTEx samples

# stats / counts params:

FDR_SIMS                         := 500
FDR_CUTOFF                       := 0.05
Cntthresh_tot                    := 6
Cntthresh_min                    := 0
#AMB_MODE                         := adjust # 'adjust' or allelic ratio diff threshold for filtering
# only 'adjust' mode only for now: add the 'weaker allele' base counts from multi-mapping reads, if any:
# all of them or until balanced with the stonger to make sure the imbalance is not caused by the multimapping reads
KEEP_CHR                         := # empty or e.g. 'X'


######################
### PIPELINE STEPS ###
######################

# todo: a better way than using these 'tmp's?
empty_string:=
ifeq ($(READS_R2),$(empty_string))
  tmp1 = $(notdir $(READS_R1))
  tmp2 = $(tmp1:.gz=)
  tmp3 = $(tmp2:.fastq=)
  PREFIX = $(tmp3:.fq=)
  FASTQC_out = $(PREFIX)_fastqc.html
else
  tmp11 = $(notdir $(READS_R1))
  tmp12 = $(tmp11:.gz=)
  tmp13 = $(tmp12:.fastq=)
  tmp14 = $(tmp13:.fq=)
  tmp21 = $(notdir $(READS_R2))
  tmp22 = $(tmp21:.gz=)
  tmp23 = $(tmp22:.fastq=)
  tmp24 = $(tmp23:.fq=)
  PREFIX = $(tmp14)_$(tmp24)
  FASTQC_out = $(tmp14)_fastqc.html
endif

ifeq ($(RM_DUPLICATE_READS),on)
  DEDUP_SUFFIX = rmdup.
endif

ifeq ($(PERFORM_FASTQC),off)
  FASTQC_out = $(empty_string)
endif

FINAL_ALIGNMENT_FILENAME = $(PREFIX)_$(ALIGNMENT_MODE)-params.Aligned.sortedByCoord.out.$(DEDUP_SUFFIX)bam
HetSNV_UNIQALNS_FILENAME = $(PREFIX)_$(ALIGNMENT_MODE)-params_crdsorted_uniqreads_over_hetSNVs.bam
HetSNV_MMAPALNS_FILENAME = $(PREFIX)_$(ALIGNMENT_MODE)-params_crdsorted_mmapreads_over_hetSNVs.bam

$(info PGENOME_DIR: $(PGENOME_DIR))
$(info READS_R1: $(READS_R1))
$(info READS_R2: $(READS_R2))
$(info PERFORM_FASTQC: $(PERFORM_FASTQC))
$(info PREFIX: $(PREFIX))
$(info ALIGNMENT_MODE: $(ALIGNMENT_MODE))
$(info RM_DUPLICATE_READS: $(RM_DUPLICATE_READS))
$(info FINAL_ALIGNMENT_FILENAME: $(FINAL_ALIGNMENT_FILENAME))
$(info HetSNV_UNIQALNS_FILENAME: $(HetSNV_UNIQALNS_FILENAME))
$(info HetSNV_MMAPALNS_FILENAME: $(HetSNV_MMAPALNS_FILENAME))
$(info $(empty_string))
$(info $(empty_string))
$(info $(empty_string))
$(info $(empty_string))



######################
### PIPELINE START ###
######################

all: $(FASTQC_out) $(PREFIX)_ref_allele_ratios.raw_counts.pdf $(PREFIX)_ref_allele_ratios.filtered_counts.pdf $(PREFIX)_ref_allele_ratios.filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min.pdf $(PREFIX)_interestingHets.FDR-$(FDR_CUTOFF).binom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv $(PREFIX)_interestingHets.FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv 


#currently, keeping alleleDB betabinomial scripts with as few modifications as possible

$(PREFIX)_interestingHets.FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv: $(PREFIX)_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv
	Rscript $(PL)/alleledb_calcOverdispersion.R \
		$< \
		$(PREFIX)_FDR-$(FDR_CUTOFF).betabinomial.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min
	Rscript $(PL)/alleledb_alleleseqBetabinomial.R \
		$< \
		$(PREFIX)_FDR-$(FDR_CUTOFF).betabinomial.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min \
		$(PREFIX)_counts.FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv \
		$@ \
		$(PREFIX)_FDR-$(FDR_CUTOFF).betabinom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.txt \
		$(FDR_CUTOFF)



$(PREFIX)_interestingHets.FDR-$(FDR_CUTOFF).binom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv: $(PREFIX)_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv
	python $(PL)/FalsePos.py $< $(FDR_SIMS) $(FDR_CUTOFF) > $(PREFIX)_FDR-$(FDR_CUTOFF).binom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.txt
	cat $< | python $(PL)/filter_by_pval.py $(PREFIX)_FDR-$(FDR_CUTOFF).binom.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.txt > $@


# allelic ratio distrs
$(PREFIX)_ref_allele_ratios.filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min.pdf: $(PREFIX)_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv
	Rscript $(PL)/plot_AllelicRatio_distribution.R $< $(PREFIX) filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min

# filter based on total counts and min per allele count
# and in non-autosomal chr, optionally keeping X;
$(PREFIX)_filtered_counts.chrs1-22$(KEEP_CHR).$(Cntthresh_tot)-tot_$(Cntthresh_min)-min_cnt.tsv: $(PREFIX)_filtered_counts.tsv
	cat $< | \
	python $(PL)/filter_non-autosomal_chr.py $(KEEP_CHR) | \
	python $(PL)/filter_by_counts.py $(Cntthresh_tot) $(Cntthresh_min) > $@

# allelic ratio distrs
$(PREFIX)_ref_allele_ratios.filtered_counts.pdf: $(PREFIX)_filtered_counts.tsv
	Rscript $(PL)/plot_AllelicRatio_distribution.R $< $(PREFIX) filtered_counts

# filter out sites in potential cnv regions 
# and sites with seemingly misphased/miscalled nearby variants
# filter/adjust sites imbalanced likely due to unaccounted multi-mapping reads 
# will use 'adjust' only for now
$(PREFIX)_filtered_counts.tsv: $(PREFIX)_raw_counts.tsv $(PREFIX)_hap1_mmapreads.mpileup $(PREFIX)_hap2_mmapreads.mpileup
	cat $< | \
	python $(PL)/filter_cnv_sites.py $(PREFIX)_discarded_HetSNVs_potential-CNV.log $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_rd.tab | \
	python $(PL)/filter_phase_warnings.py $(PREFIX)_discarded_HetSNVs_warn-haplotype.log | \
	#python $(PL)/filter_sites_w_mmaps.py $(AMB_MODE) $(PREFIX)_discarded_HetSNVs_warn-mmaps.log $(PREFIX)_mmap_reads_over_hetSNVs.log \
	python $(PL)/filter_sites_w_mmaps.py adjust $(PREFIX)_discarded_HetSNVs_warn-mmaps.log $(PREFIX)_mmap_reads_over_hetSNVs.log \
		$(PREFIX)_hap1_mmapreads.mpileup $(PREFIX)_hap2_mmapreads.mpileup > $@



# allelic ratio distrs
$(PREFIX)_ref_allele_ratios.raw_counts.pdf: $(PREFIX)_raw_counts.tsv
	Rscript $(PL)/plot_AllelicRatio_distribution.R $< $(PREFIX) raw_counts

# counts
$(PREFIX)_raw_counts.tsv: $(PREFIX)_hap1_uniqreads.mpileup $(PREFIX)_hap2_uniqreads.mpileup
	python $(PL)/pileup2counts.py 1 $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_ref.bed \
	$(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_hap1.bed $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_hap2.bed \
	$(PREFIX)_discarded_HetSNVs_from-pileups.log \
	$(PREFIX)_hap1_uniqreads.mpileup $(PREFIX)_hap2_uniqreads.mpileup > $@


# pileups
$(PREFIX)_hap1_mmapreads.mpileup: $(HetSNV_MMAPALNS_FILENAME)
	$(SAMTOOLS) mpileup -BQ0 --max-depth 999999 --ff UNMAP -f $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hap1.fa $< \
	--positions $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_hap1.bed --output $@

$(PREFIX)_hap2_mmapreads.mpileup: $(HetSNV_MMAPALNS_FILENAME)
	$(SAMTOOLS) mpileup -BQ0 --max-depth 999999 --ff UNMAP -f $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hap2.fa $< \
	--positions $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_hap2.bed --output $@

$(PREFIX)_hap1_uniqreads.mpileup: $(HetSNV_UNIQALNS_FILENAME)
	$(SAMTOOLS) mpileup -BQ0 --max-depth 999999 --ff UNMAP -f $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hap1.fa $< \
	--positions $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_hap1.bed --output $@

$(PREFIX)_hap2_uniqreads.mpileup: $(HetSNV_UNIQALNS_FILENAME)
	$(SAMTOOLS) mpileup -BQ0 --max-depth 999999 --ff UNMAP -f $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hap2.fa $< \
	--positions $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_hap2.bed --output $@

# by default --ff was also filtering-out some other - secondary? reads:  [UNMAP,SECONDARY,QCFAIL,DUP], leaving only UNMAP for now




# non-uniq alns over hetSNVs:
$(PREFIX)_$(ALIGNMENT_MODE)-params_crdsorted_mmapreads_over_hetSNVs.bam: $(FINAL_ALIGNMENT_FILENAME)
	cat $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_hap1.bed $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_hap2.bed | \
	$(SAMTOOLS) view -h -L - $< | awk '$$5!="255" {print $$0}' | \
	$(SAMTOOLS) view -b - > $@
	$(SAMTOOLS) index $@
	$(SAMTOOLS) flagstat $@ > $@.stat

# uniq alns over hetSNVs:
$(PREFIX)_$(ALIGNMENT_MODE)-params_crdsorted_uniqreads_over_hetSNVs.bam: $(FINAL_ALIGNMENT_FILENAME)
	cat $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_hap1.bed $(PGENOME_DIR)/$(VCF_SAMPLE_ID)_hetSNVs_hap2.bed | \
	$(SAMTOOLS) view -h -q 255 -L - $< | \
	$(SAMTOOLS) view -b - > $@
	$(SAMTOOLS) index $@
	$(SAMTOOLS) flagstat $@ > $@.stat


# if removing duplicate reads
$(PREFIX)_$(ALIGNMENT_MODE)-params.Aligned.sortedByCoord.out.rmdup.bam: $(PREFIX)_$(ALIGNMENT_MODE)-params.Aligned.sortedByCoord.out.bam
	$(JAVA) -Xmx$(JAVA_MEM) -jar $(PICARD) MarkDuplicates \
	INPUT=$< OUTPUT=$@ METRICS_FILE=$(@:.bam=.metrics) \
	REMOVE_DUPLICATES=true \
	DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES
	$(SAMTOOLS) index $@
	$(SAMTOOLS) flagstat $@ > $@.stat


## specific additional params will be read from $(STAR_parameters_file)
$(PREFIX)_custom-params.Aligned.sortedByCoord.out.bam: $(READS_R1)
	$(STAR) \
	--runThreadN $(NTHR) \
	--genomeDir $(GenomeIdx_STAR_diploid) \
	--readFilesIn $< $(READS_R2) \
	--readFilesCommand $(STAR_readFilesCommand) \
	--outFileNamePrefix $(@:Aligned.sortedByCoord.out.bam=) \
	--outSAMattributes All \
	--outFilterMultimapNmax 999999 \
	--scoreGenomicLengthLog2scale 0.0 \
	--sjdbOverhang $(STAR_sjdbOverhang) \
	--sjdbGTFfile $(Annotation_diploid) \
	--parametersFiles $(STAR_parameters_file) \
	--limitSjdbInsertNsj $(STAR_limitSjdbInsertNsj) \
	--outSAMtype BAM SortedByCoordinate
	$(SAMTOOLS) flagstat $@ > $@.stat
	$(SAMTOOLS) index $@

## params for RNA-seq; will use as default for ASE
$(PREFIX)_ASE-params.Aligned.sortedByCoord.out.bam: $(READS_R1)
	$(STAR) \
	--runThreadN $(NTHR) \
	--genomeDir $(GenomeIdx_STAR_diploid) \
	--twopassMode Basic \
	--readFilesIn $< $(READS_R2) \
	--readFilesCommand $(STAR_readFilesCommand) \
	--outFileNamePrefix $(@:Aligned.sortedByCoord.out.bam=) \
	--outSAMattributes All \
	--outFilterMismatchNoverReadLmax $(STAR_outFilterMismatchNoverReadLmax) \
	--outFilterMatchNminOverLread $(STAR_outFilterMatchNminOverLread) \
	--outFilterMultimapNmax 999999 \
	--scoreGenomicLengthLog2scale 0.0 \
	--sjdbOverhang $(STAR_sjdbOverhang) \
	--sjdbGTFfile $(Annotation_diploid) \
	--limitSjdbInsertNsj $(STAR_limitSjdbInsertNsj) \
	--outSAMtype BAM SortedByCoordinate
	$(SAMTOOLS) flagstat $@ > $@.stat
	$(SAMTOOLS) index $@	

## opts similar to AlleleSeq v1.2a bowtie1 -v 2 -m 1 mode; but with small gaps allowed, no splicing; will use as default for ASB
$(PREFIX)_ASB-params.Aligned.sortedByCoord.out.bam: $(READS_R1)
	$(STAR) \
	--runThreadN $(NTHR) \
	--genomeDir $(GenomeIdx_STAR_diploid) \
	--readFilesIn $< $(READS_R2) \
	--readFilesCommand $(STAR_readFilesCommand) \
	--outFileNamePrefix $(@:Aligned.sortedByCoord.out.bam=) \
	--outSAMattributes All \
	--outFilterMismatchNoverReadLmax $(STAR_outFilterMismatchNoverReadLmax) \
	--outFilterMatchNminOverLread $(STAR_outFilterMatchNminOverLread) \
	--outFilterMultimapNmax 999999 \
	--scoreGap -100 \
	--scoreGenomicLengthLog2scale 0.0 \
	--sjdbScore 0 \
	--limitSjdbInsertNsj $(STAR_limitSjdbInsertNsj) \
	--outSAMtype BAM SortedByCoordinate
	$(SAMTOOLS) flagstat $@ > $@.stat
	$(SAMTOOLS) index $@	


## opts for ASCA, atac-seq, similar to ASB, but will require adapter-trimmed reads
$(PREFIX)_ASCA-params.Aligned.sortedByCoord.out.bam: $(READS_R1)
	$(STAR) \
	--runThreadN $(NTHR) \
	--genomeDir $(GenomeIdx_STAR_diploid) \
	--readFilesIn $< $(READS_R2) \
	--readFilesCommand $(STAR_readFilesCommand) \
	--outFileNamePrefix $(@:Aligned.sortedByCoord.out.bam=) \
	--outSAMattributes All \
	--outFilterMismatchNoverReadLmax $(STAR_outFilterMismatchNoverReadLmax) \
	--outFilterMatchNminOverLread $(STAR_outFilterMatchNminOverLread) \
	--outFilterMultimapNmax 999999 \
	--scoreGap -100 \
	--scoreGenomicLengthLog2scale 0.0 \
	--sjdbScore 0 \
	--limitSjdbInsertNsj $(STAR_limitSjdbInsertNsj) \
	--outSAMtype BAM SortedByCoordinate
	$(SAMTOOLS) flagstat $@ > $@.stat
	$(SAMTOOLS) index $@	

$(READS_R1).trimmed.fastq.gz $(READS_R2).trimmed.fastq.gz: $(READS_R1) $(READS_R2)
	$(CUTADAPT) -m 5 \
	-a $(R1_ADAPTER_SEQ) \
	-A $(R2_ADAPTER_SEQ) \
	-o $(READS_R1).trimmed.fastq.gz \
	-p $(READS_R2).trimmed.fastq.gz \
	$(READS_R1) $(READS_R2)
	$(FASTQC) --threads $(NTHR) $(READS_R1).trimmed.fastq.gz $(READS_R2).trimmed.fastq.gz


$(FASTQC_out): $(READS_R1)
	$(FASTQC) --threads $(NTHR) $(READS_R1) $(READS_R2)


