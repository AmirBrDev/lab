__author__ = 'amirbar'


files = ["assign-type-to-all-chimeras-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_all_interactions.with-type",
         "assign-type-to-all-chimeras-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_all_interactions.with-type",
         "assign-type-to-all-chimeras-of-MG_hfq-WT101_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.with-type",
         "assign-type-to-all-chimeras-of-MG_hfq-wt202_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.with-type",
         "assign-type-to-all-chimeras-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_all_interactions.with-type",
         "assign-type-to-signif-chimeras-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_sig_interactions.with-type",
         "assign-type-to-signif-chimeras-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_sig_interactions.with-type",
         "assign-type-to-signif-chimeras-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_sig_interactions.with-type",
         "assign-type-to-single-counts-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_single_counts.with-type",
         "assign-type-to-single-counts-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_single_counts.with-type",
         "assign-type-to-single-counts-of-MG_hfq-WT101_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.with-type",
         "assign-type-to-single-counts-of-MG_hfq-wt202_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.with-type",
         "assign-type-to-single-counts-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_single_counts.with-type"]

from TableCrosser import cross_zhang, cross_raghavan, cross_bilusic, cross_tss, cross_lybecker
from TableFormatter import TableFormatter
from TableLoader import TableLoader

for name in files:
    print name
    print "*" * 100
    # TableFormatter.split_reverse_rows(TableLoader, "results/raghavan/reverse/%s.table" % name)
    # TableFormatter.split_reverse_rows(TableLoader, "results/zhang/reverse/%s.table" % name)
    TableFormatter.split_reverse_rows(TableLoader, "results/bilusic/reverse/%s.table" % name)
    # TableFormatter.split_reverse_rows(TableLoader, "results/tss/reverse/%s.table" % name)
    # TableFormatter.split_reverse_rows(TableLoader, "results/lybecker/reverse/%s.table" % name)
