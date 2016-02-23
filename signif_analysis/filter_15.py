import csv

file_list = ["/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-no-unified-lib-counts/with-iron-no-unified-count-libs-assign-type-to-signif-chimeras-of-MG_hfq-FLAG207_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_sig_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-no-unified-lib-counts/with-iron-no-unified-count-libs-assign-type-to-signif-chimeras-of-MG_hfq-FLAG208_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_sig_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-no-unified-lib-counts/with-iron-no-unified-count-libs-assign-type-to-signif-chimeras-of-MG_hfq-FLAG305-CL-Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_sig_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-no-unified-lib-counts/with-stat-no-unified-count-libs-assign-type-to-signif-chimeras-of-MG_hfq-FLAG209_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_sig_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-no-unified-lib-counts/with-stat-no-unified-count-libs-assign-type-to-signif-chimeras-of-MG_hfq-FLAG210_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_sig_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-no-unified-lib-counts/with-stat-no-unified-count-libs-assign-type-to-signif-chimeras-of-MG_hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_sig_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-no-unified-lib-counts/with-log-no-unified-count-libs-assign-type-to-signif-chimeras-of-MG_hfq-FLAG101_A_T1_7m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_sig_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-no-unified-lib-counts/with-log-no-unified-count-libs-assign-type-to-signif-chimeras-of-MG_hfq-FLAG102_A_T1_15m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_sig_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-no-unified-lib-counts/with-log-no-unified-count-libs-assign-type-to-signif-chimeras-of-MG_hfq-FLAG103_A_T1_5m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_sig_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-no-unified-lib-counts/with-log-no-unified-count-libs-assign-type-to-signif-chimeras-of-MG_hfq-FLAG104_A_T1_3m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_sig_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-no-unified-lib-counts/with-log-no-unified-count-libs-assign-type-to-signif-chimeras-of-MG_hfq-FLAG108_IP_50ul_Beads_cutadapt_bwa.bam_all_fragments_l25.txt_sig_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-no-unified-lib-counts/with-log-no-unified-count-libs-assign-type-to-signif-chimeras-of-MG_hfq-FLAG109_lysed_by_MM400_cutadapt_bwa.bam_all_fragments_l25.txt_sig_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-all-chim-lib-counts/with-iron-count-all-chim-libs-assign-type-to-all-interactions-of-MG_hfq-FLAG207_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-all-chim-lib-counts/with-iron-count-all-chim-libs-assign-type-to-all-interactions-of-MG_hfq-FLAG208_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-all-chim-lib-counts/with-iron-count-all-chim-libs-assign-type-to-all-interactions-of-MG_hfq-FLAG305-CL-Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-all-chim-lib-counts/with-stat-count-all-chim-libs-assign-type-to-all-interactions-of-MG_hfq-FLAG209_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-all-chim-lib-counts/with-stat-count-all-chim-libs-assign-type-to-all-interactions-of-MG_hfq-FLAG210_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-all-chim-lib-counts/with-stat-count-all-chim-libs-assign-type-to-all-interactions-of-MG_hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-all-chim-lib-counts/with-log-count-all-chim-libs-assign-type-to-all-interactions-of-MG_hfq-FLAG101_A_T1_7m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-all-chim-lib-counts/with-log-count-all-chim-libs-assign-type-to-all-interactions-of-MG_hfq-FLAG102_A_T1_15m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-all-chim-lib-counts/with-log-count-all-chim-libs-assign-type-to-all-interactions-of-MG_hfq-FLAG103_A_T1_5m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-all-chim-lib-counts/with-log-count-all-chim-libs-assign-type-to-all-interactions-of-MG_hfq-FLAG104_A_T1_3m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-all-chim-lib-counts/with-log-count-all-chim-libs-assign-type-to-all-interactions-of-MG_hfq-FLAG108_IP_50ul_Beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-all-chim-lib-counts/with-log-count-all-chim-libs-assign-type-to-all-interactions-of-MG_hfq-FLAG109_lysed_by_MM400_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-single-lib-counts/with-iron-count-single-libs-assign-type-to-only-single-of-MG_hfq-FLAG207_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_only_singles.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-single-lib-counts/with-iron-count-single-libs-assign-type-to-only-single-of-MG_hfq-FLAG208_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_only_singles.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-single-lib-counts/with-iron-count-single-libs-assign-type-to-only-single-of-MG_hfq-FLAG305-CL-Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_only_singles.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-single-lib-counts/with-stat-count-single-libs-assign-type-to-only-single-of-MG_hfq-FLAG209_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_singles.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-single-lib-counts/with-stat-count-single-libs-assign-type-to-only-single-of-MG_hfq-FLAG210_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_singles.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-single-lib-counts/with-stat-count-single-libs-assign-type-to-only-single-of-MG_hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_singles.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-single-lib-counts/with-log-count-single-libs-assign-type-to-only-single-of-MG_hfq-FLAG101_A_T1_7m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_only_singles.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-single-lib-counts/with-log-count-single-libs-assign-type-to-only-single-of-MG_hfq-FLAG102_A_T1_15m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_only_singles.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-single-lib-counts/with-log-count-single-libs-assign-type-to-only-single-of-MG_hfq-FLAG103_A_T1_5m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_only_singles.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-single-lib-counts/with-log-count-single-libs-assign-type-to-only-single-of-MG_hfq-FLAG104_A_T1_3m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_only_singles.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-single-lib-counts/with-log-count-single-libs-assign-type-to-only-single-of-MG_hfq-FLAG108_IP_50ul_Beads_cutadapt_bwa.bam_all_fragments_l25.txt_only_singles.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/tables-with-single-lib-counts/with-log-count-single-libs-assign-type-to-only-single-of-MG_hfq-FLAG109_lysed_by_MM400_cutadapt_bwa.bam_all_fragments_l25.txt_only_singles.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/with-iron-no-unified-count-libs-assign-type-to-all-interactions-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_all_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/with-log-no-unified-count-libs-assign-type-to-all-interactions-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_all_interactions.txt.with-known",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/with-stat-no-unified-count-libs-assign-type-to-all-interactions-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_all_interactions.txt.with-known"]

for file_name in file_list:
	
	with open(file_name, "rb") as fl:
		line = fl.readline()

		reader = csv.reader(fl, delimiter='\t', quotechar='|')
		
		with open(file_name + ".threshold_15", "wb") as flw: 
			writer = csv.writer(flw, delimiter="\t", quotechar='|')
			writer.writerow(line[:-1].split("\t"))

			for line in reader:
				if int(line[16]) >= 15:
					writer.writerow(line)
