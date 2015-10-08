from Table import Table
from sequence_extractor import get_seq_list
from alignment_handler import run_local_alignment, \
							  prepare_for_multiple_aignment, \
							  run_multiple_sequence_alignment
from Rate4SiteRunner import Rate4Site
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def run():
	seq_list = get_seq_list("./NC_000913.fna.gz", "./candidate_list.table")

	# seq_list holds tuples of (name, seq_obj)
	for seq in seq_list:
		print seq

		# Generate sequence file
		seq_file_name = "./tmp/%s.fasta" % seq[0]
		base_record = SeqRecord(seq[1], id=seq[0], name=seq[0])
		SeqIO.write(base_record, seq_file_name, "fasta")

		# Create local alignments for the sequence
		local_alignment_list = \
			run_local_alignment(seq_file_name, "./genome_list.table")

		print "*" * 100
		print "local alignments results"
		print "*" * 100

		for row in local_alignment_list:
			print row

		# Exctract only the significant alignments
		records_for_emma = \
			prepare_for_multiple_aignment(base_record, local_alignment_list)

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
										"tmp/emma.aln",
										"tmp/emma.dnd")

		# Run Rate4Site
		rate_runner = Rate4Site("tmp/emma.aln", "tmp/emma.dnd")

if __name__ == "__main__":
	run()