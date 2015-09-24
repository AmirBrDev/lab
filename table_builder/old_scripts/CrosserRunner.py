__author__ = 'amirbar'

import subprocess
import sys

def execute(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # Poll process for new output until finished
    while True:
        nextline = process.stdout.readline()
        if nextline == '' and process.poll() != None:
            break
        sys.stdout.write(nextline)
        sys.stdout.flush()

    output = process.communicate()[0]
    exitCode = process.returncode

    if (exitCode == 0):
        return output
    else:
        print "running failed"

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

# for name in files:
#     execute("./TableCrosser.py")

from TableCrosser import cross_zhang, cross_raghavan, cross_bilusic, cross_tss, cross_lybecker

for name in files:
    # cross_raghavan(name, reverse=True)
    # cross_zhang(name, reverse=True)
    # cross_bilusic(name, reverse=True)
    cross_tss(name, reverse=True)
    # cross_lybecker(name, reverse=True)