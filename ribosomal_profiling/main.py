import csv
import os

SUMMARY_FILE_TARGET_NAME_INDEX = 1
SUMMARY_FILE_TARGET_VALUE_INDEX = 8

SUMMARY_FILE_PATH = "./summary_of_RyhB_riobosome_profiling.csv"

summary_target_to_val = {}

with open(SUMMARY_FILE_PATH, "rb") as fl:

	fl.readline()

	reader = csv.reader(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

	for row in reader:
		if row[SUMMARY_FILE_TARGET_NAME_INDEX] == "":
			continue

		# print row[SUMMARY_FILE_TARGET_NAME_INDEX], row[SUMMARY_FILE_TARGET_VALUE_INDEX]
		summary_target_to_val[row[SUMMARY_FILE_TARGET_NAME_INDEX]] = \
			row[SUMMARY_FILE_TARGET_VALUE_INDEX]


RILSEQ_FILE_PATH = "RyhB_RILseq_targets_Iron_lim.csv"
OUTPUT = "./result.csv"

with open(RILSEQ_FILE_PATH, "rb") as fl:
	with open(OUTPUT, "wb") as fout:
		fl.readline()

		reader = csv.reader(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer = csv.writer(fout, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

		writer.writerow(["target", "NOR", "I"])

		for row in reader:
			try:
				writer.writerow([row[0], row[2], summary_target_to_val[row[0]]])
			except KeyError:
				writer.writerow([row[0], row[2], "-"])

	# with open("/tmp/%s" % name, "wb") as fout:
	# 	fout.write("\t".join(fl.readline().replace("\n", "").split("\t") + ["code"]) + "\n")

	# 	reader = csv.reader(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
	# 	writer = csv.writer(fout, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		
	# 	for row in reader:
	# 		if row[4] > row[5]:
	# 			name = "%s::%s" % (row[5], row[4])
	# 		else:
	# 			name = "%s::%s" % (row[4], row[5])
	# 		name = name.replace(".5UTR", "").replace(".EST5UTR", "")
	# 		row.append(name)
	# 		writer.writerow(row)