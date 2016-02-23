import csv

file_list = ["/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG101_A_T1_7m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG101_A_T1_7m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG102_A_T1_15m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG102_A_T1_15m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG103_A_T1_5m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG103_A_T1_5m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG104_A_T1_3m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG104_A_T1_3m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG107_TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG107_TAP_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG108_IP_50ul_Beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG108_IP_50ul_Beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG109_lysed_by_MM400_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG109_lysed_by_MM400_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG207_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG207_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG208_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG208_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG209_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG209_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG210_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG210_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG305-CL-Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG305-CL-Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/users/amirbar/lab/table_builder/our_files/refactor/set_1/MG_hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt]"]

for file_name in file_list:
	
	with open(file_name, "rb") as fl:
		line = fl.readline()

		reader = csv.reader(fl, delimiter='\t', quotechar='|')
		
		with open(file_name + ".threshold_15", "wb") as flw: 
			writer = csv.writer(flw, delimiter="\t", quotechar='|')
			writer.write(line)

			for line in reader:
				if int(line[14]) >= 15:
					writer.write(line)