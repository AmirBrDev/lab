from Table import Table
from sequence_extractor import get_seq_list
from alignment_handler import run_local_alignment, \
							  prepare_for_multiple_aignment, \
							  run_multiple_sequence_alignment, \
							  build_project_id_to_tax_id_dictionary, \
							  build_secondary_project_id_to_tax_id_dictionary, \
							  generate_score_histogram
from Rate4SiteRunner import Rate4Site
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def run(base_genome_file, 
		candidate_list_file, 
		workdir, 
		genome_list_file,
		project_to_taxamonies_file,
		tree_file):

	seq_list = get_seq_list(base_genome_file, candidate_list_file)
	conv_table = build_project_id_to_tax_id_dictionary(project_to_taxamonies_file)
	secondary_conv_table = \
		build_secondary_project_id_to_tax_id_dictionary(project_to_taxamonies_file)

	# seq_list holds tuples of (name, seq_obj)
	for seq in seq_list:
		print seq

		# Generate sequence file
		seq_file_name = "%s/%s.fasta" % (workdir, seq[0])
		base_record = SeqRecord(seq[1], id=seq[0], name=seq[0])
		SeqIO.write(base_record, seq_file_name, "fasta")

		# Create local alignments for the sequence
		local_alignment_list = \
			run_local_alignment(seq_file_name, genome_list_file, WORK_DIR)

		print "*" * 100
		print "local alignments results"
		print "*" * 100

		for row in local_alignment_list:
			print row

		# Generate histogram
		generate_score_histogram("%s/stats/%s.png" % (workdir, seq[0]), local_alignment_list)

		# Exctract only the significant alignments
		records_for_emma = \
			prepare_for_multiple_aignment(local_alignment_list,
										  conv_table,
										  secondary_conv_table)

		print "*" * 100
		print "records for emma"
		print "*" * 100

		for record in records_for_emma:
			print record

		print "*" * 100
		print "Running Emma"
		print "*" * 100

		# Run emma
		run_multiple_sequence_alignment(records_for_emma,
										"%s/emma.aln" % WORK_DIR,
										"%s/emma.dnd" % WORK_DIR,
										WORK_DIR)

		# Run Rate4Site
		rate_runner = Rate4Site("%s/emma.aln" % workdir, 
								tree_file,
								"./rate4site64")

		rate_runner.runRate(outname="results/%s.rate" % seq[0])



if __name__ == "__main__":
	
	BASE_GENOME_FILE = "./NC_000913.fna.gz"
	CANDIDATE_LIST_FILE = "./candidate_list.table"
	WORK_DIR = "./tmp"
#	GENOME_LIST_FILE = "./genome_list.table"
	GENOME_LIST_FILE = "./genome_list_exist.table"
#	GENOME_LIST_FILE = "./genome_list_partial.table"
	PROJ_TO_TAXAMONIES_FILE = "./microbial/lproks_1.txt.gz"
	TREE_FILE = "speciesTree.cgi"

	run(BASE_GENOME_FILE,
		CANDIDATE_LIST_FILE,
		WORK_DIR,
		GENOME_LIST_FILE,
		PROJ_TO_TAXAMONIES_FILE,
		TREE_FILE)


#	rate_runner = Rate4Site("%s/emma.aln" % WORK_DIR, 
#							TREE_FILE,
#							"./rate4site64")
#
#	print rate_runner.runRate()









#from alignment_handler import get_id_by_name
#import gzip
#

#
#	conversion_table = build_project_id_to_tax_id_dictionary(PROJ_TO_TAXAMONIES_FILE)
#	sec_table = build_secondary_project_id_to_tax_id_dictionary(PROJ_TO_TAXAMONIES_FILE)
#
#	with open(GENOME_LIST_FILE, "rb") as fl:
#		header = fl.readline()
#
#		for name in fl.readlines():
#
#			try:
#				project_id = get_id_by_name(name.strip(), conversion_table)
##				print "#" * 50
##				print "matched project id %s to name: %s" % (project_id, name)
#
#			except BaseException as ex:
#
#				try:
#					project_id = get_id_by_name(name.strip(), sec_table)
#					print "#" * 50
#					print "matched project id  %s to %s" % (project_id, name)
#
#				except BaseException as ex:
#					print "*" * 50
#					print "Unmatched project id to %s" % name
