files = ["/home/users/raya/deepSeq/RILseq/Sahar/final_for_paper/tables_with_RNAup/assign-type-to-RNAup-unified-of-flipped-rnaup-MG_hfq-FLAG101-104_108_109_unified_Log_bwa.bam_all_fragments_l25.txt_all_interactions_with-total.txt.with-known",
		 "/home/users/raya/deepSeq/RILseq/Sahar/final_for_paper/tables_with_RNAup/assign-type-to-RNAup-unified-of-flipped-rnaup-MG_hfq-FLAG207_208_305_unified_Iron_bwa.bam_all_fragments_l25.txt_all_interactions_with-total.txt.with-known",
		 "/home/users/raya/deepSeq/RILseq/Sahar/final_for_paper/tables_with_RNAup/assign-type-to-RNAup-unified-of-flipped-rnaup-MG_hfq-FLAG209_210_312_bwa.bam_unified_Stationary_all_fragments_l25.txt_all_interactions_with-total.txt.with-known"]

def foo():

	new_names = []

	for path in files:

		with open(path, "rb") as fl:
			lines = fl.readlines()

		for line in lines[1:]:
			sline = line.split("\t")
			
			if sline[4]>sline[5]:
				new_names.append("%s::%s" % (sline[5], sline[4]))
			else:
				new_names.append("%s::%s" % (sline[4], sline[5]))


	new_names = set(new_names)


	missed_names = []

	with open("/home/users/amirbar/lab/signif_analysis/table.csv", "rb") as fl:
		lines = fl.readlines()

	for line in lines[1:]:
		sline = line.split("\t")

		print sline
		if sline[1] == '0':
			missed_names.append(sline[0])

	missed_names = set(missed_names)

	print "missed:", len(missed_names)
	print "new_names:", len(new_names)
	print "intersection:", len(new_names.intersection(missed_names))

import csv
import os


# First get the known names
known_names = []

with open("/home/users/amirbar/lab/signif_analysis/Table S3-Known-srna-target-chimeras.csv", "rb") as fl:

	line = fl.readline()

	for line in fl:
		sline = line.strip().split("\t")

		known_names.append(sline[0])
		known_names.append(sline[1])

known_names = list(set(known_names))

for path in files:

	name = "%s_with_code.csv" % os.path.basename(path)

	with open(path, "rb") as fl:
		with open("/home/hosts/disk19/amirbar/our_tables/24_3_2016/all_interactions/%s" % name, "wb") as fout:
			fout.write("\t".join(fl.readline().strip().split("\t") + ["code"]) + "\n")

			reader = csv.reader(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			writer = csv.writer(fout, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			
			for row in reader:
				first = row[4]
				second = row[5]

				# tu should be named according to known genes
				if first.lower().endswith(".tu"):
					for name_part in first.split(".")[:2]:
						if name_part in known_names:
							first = name_part
							break

				if second.lower().endswith(".tu"):
					for name_part in first.split(".")[:2]:
						if name_part in known_names:
							second = name_part
							break

				if row[4] > row[5]:
					name = "%s::%s" % (second, first)
				else:
					name = "%s::%s" % (first, second)

				name = name.replace(".5UTR", "").replace(".EST5UTR", "")
				row.append(name)
				writer.writerow(row)
