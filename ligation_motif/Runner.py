__author__ = 'amirbar'


import subprocess
import os

# file_list= [# ("101", "MG_hfq-FLAG101_A_T1_7m_22C_cutadapt_bwa.bam", "MG_hfq-FLAG101_A_T1_7m_22C_cutadapt_bwa.bam_all_fragments_l25.txt"),
#             # ("102", "MG_hfq-FLAG102_A_T1_15m_22C_cutadapt_bwa.bam", "MG_hfq-FLAG102_A_T1_15m_22C_cutadapt_bwa.bam_all_fragments_l25.txt"),
#             # ("103", "MG_hfq-FLAG103_A_T1_5m_22C_cutadapt_bwa.bam", "MG_hfq-FLAG103_A_T1_5m_22C_cutadapt_bwa.bam_all_fragments_l25.txt")]
#             # ("104", "MG_hfq-FLAG104_A_T1_3m_22C_cutadapt_bwa.bam", "MG_hfq-FLAG104_A_T1_3m_22C_cutadapt_bwa.bam_all_fragments_l25.txt"),
#             # ("108", "MG_hfq-FLAG108_IP_50ul_Beads_cutadapt_bwa.bam", "MG_hfq-FLAG108_IP_50ul_Beads_cutadapt_bwa.bam_all_fragments_l25.txt"),
#             # ("109", "MG_hfq-FLAG109_lysed_by_MM400_cutadapt_bwa.bam", "MG_hfq-FLAG109_lysed_by_MM400_cutadapt_bwa.bam_all_fragments_l25.txt")]
#             # ("207", "MG_hfq-FLAG207_CL_Iron_limitation_cutadapt_bwa.bam", "MG_hfq-FLAG207_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt"),
#             # ("208", "MG_hfq-FLAG208_CL_Iron_limitation_cutadapt_bwa.bam", "MG_hfq-FLAG208_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt"),
#             # ("209", "MG_hfq-FLAG209_CL_Stationary_cutadapt_bwa.bam", "MG_hfq-FLAG209_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt")]
#             ("210", "MG_hfq-FLAG210_CL_Stationary_cutadapt_bwa.bam", "MG_hfq-FLAG210_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt"),
#             ("305", "MG_hfq-FLAG305-CL-Iron_limitation_cutadapt_bwa.bam", "MG_hfq-FLAG305-CL-Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt"),
#             ("312", "MG_hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam", "MG_hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt")]
#             # ("combined", "combined.bam", "combined_all_fragments_l25.txt")]

file_list = [
			 ("101", "Total_MG-hfq-FLAG101-CL-Log_cutadapt_bwa.bam", "Total_MG-hfq-FLAG101-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt"),
			 ("102", "Total_MG-hfq-FLAG102-CL-Log_cutadapt_bwa.bam", "Total_MG-hfq-FLAG102-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt"),
			 ("103", "Total_MG-hfq-FLAG103-CL-Log_cutadapt_bwa.bam", "Total_MG-hfq-FLAG103-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt"),
			 ("104", "Total_MG-hfq-FLAG104-CL-Log_cutadapt_bwa.bam", "Total_MG-hfq-FLAG104-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt"),
			 ("108", "Total_MG-hfq-FLAG108-CL-Log_cutadapt_bwa.bam", "Total_MG-hfq-FLAG108-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt"),
			 ("109", "Total_MG-hfq-FLAG109-CL-Log_cutadapt_bwa.bam", "Total_MG-hfq-FLAG109-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt"),
			 ("207", "Total_MG-hfq-FLAG207-CL-Log_cutadapt_bwa.bam", "Total_MG-hfq-FLAG207-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt"),
			 ("208", "Total_MG-hfq-FLAG208-CL-Log_cutadapt_bwa.bam", "Total_MG-hfq-FLAG208-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt"),
			 ("305", "Total_MG-hfq-FLAG305-CL-Log_cutadapt_bwa.bam", "Total_MG-hfq-FLAG305-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt"),
			 ("312", "Total_MG-hfq-FLAG312-CL-Log_cutadapt_bwa.bam", "Total_MG-hfq-FLAG312-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt"),
			]

genome_file = "./genome.fa"

proc_list = []
fl_list = []

for workdir, bam_file, bam_fragments_file in file_list:
       os.system("mkdir -p %s" % workdir)

       cmd = ["python",
              "../detect_fusion_point.py",
              "-l",
              bam_fragments_file,
              "-b",
              bam_file,
              "-g",
              genome_file,
              "-p",
              "/".join([workdir, "reads_fusion.txt"]),
              "-o",
              "10",
              "-w",
              "1000"]

       print "*" * 100
       print " ".join(cmd)

       fl = open("/".join([workdir, "stdout.txt"]), "wb")
       p = subprocess.Popen(cmd, stdout=fl, stderr=subprocess.PIPE)
       proc_list.append(p)
       fl_list.append(fl)
       # p.wait()

       # print "*" * 100
       # print "err"
       # print "*" * 100
       #
       # for line in p.stderr.readlines():
       #     print line

for p in proc_list:
       p.wait()

for fl in fl_list:
       fl.close()