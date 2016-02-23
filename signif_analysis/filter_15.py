import csv

file_list = ["/home/amirbar/lab/table_builder/our_files/refactor/set_1/Total_MG-hfq-FLAG101-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/set_1/Total_MG-hfq-FLAG101-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG102-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG102-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG103-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG103-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG104-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG104-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG107-CL-Log-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG107-CL-Log-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG108-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG108-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG109-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG109-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG207-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG207-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG208-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG208-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG211-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG211-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG212-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG212-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG305-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG305-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
"/home/amirbar/lab/table_builder/our_files/refactor/Total_MG-hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt"]

for file_name in file_list:
	
	with open(file_name, "rb") as fl:
		reader = csv.reader(fl, delimiter='\t', quotechar='|')
		
		line = reader.readline()
		with open(file_name + ".threshold_15", "wb") as flw: 
			writer = csv.writer(flw, delimiter="\t", quotechar='|')
			writer.write(line)

			for line in reader:
				if int(line[14]) >= 15:
					writer.write(line)
