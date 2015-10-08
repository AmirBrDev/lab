from Table import Table
from sequence_extractor import get_seq_list
from Bio.Emboss.Applications import NeedleCommandline, WaterCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from EmbossCommands import EmmaCommandline, SuperMatcherCommandline
import gzip
import os
import numpy

def run_local_alignment(seq_file, genome_list):

	genome_table = Table(genome_list)
	result = []

	for row in genome_table:
		genome_name = row["genome_file"]

		# Find the genome file according to the genome directory
		for file_name in os.listdir("./microbial/%s" % genome_name):
			if file_name.endswith(".fna.gz"):
				genome_file = "./microbial/%s/%s" % (genome_name, file_name)

		with gzip.open(genome_file, "rb") as fl:
			records = SeqIO.read(fl, "fasta")
			SeqIO.write(records, "./tmp/curr_genome.fasta", "fasta")

		alignment_file = "tmp/aln/%s.aln" % genome_name

		# Build the matching command
		cmd = SuperMatcherCommandline(asequence=seq_file, 
									  bsequence="./tmp/curr_genome.fasta", 
									  gapopen=10, 
									  gapextend=0.5,
									  outfile=alignment_file)

		# Excecute the command
		stdout, stderr = cmd()
		
		#print(stdout + stderr)

		# Get the score for the genome and the matching alignment
		align_seq = list(AlignIO.read(alignment_file, "emboss"))[1]
		align_seq._set_seq(align_seq.seq.ungap("-"))  # remove the gaps
		score = extract_score(alignment_file)

		result.append((genome_name, score, align_seq))


	os.remove("./tmp/curr_genome.fasta")

	return result


def extract_score(file_name):

	with open(file_name, "rb") as fl:
		
		line = fl.readline()

		while (not line.startswith("# Score: ")):
			line = fl.readline()

	return float(line.replace("# Score: ", ""))


def prepare_for_multiple_aignment(base_record, local_alignment_list):
	NAME_INDEX = 0
	SCORE_INDEX = 1
	RECORD_INDEX = 2

	score_list = [genome_tuple[SCORE_INDEX] for genome_tuple in local_alignment_list]

	median = numpy.median(score_list)
	std = numpy.std(score_list)

	segnificant_alignments_list = []

	# Find all the segnificant alignments
	for genome_tuple in local_alignment_list:
		name, score, record = genome_tuple

		if is_score_segnificant(score, median, std):
			segnificant_alignments_list.append(genome_tuple)

	# Each of the segnificant alignments should be added to a fasta file
	# to run within emma, here we create the records for writing
	records = []

	for genome_tuple in segnificant_alignments_list:
		record = genome_tuple[RECORD_INDEX]
		record.id = genome_tuple[NAME_INDEX]
		record.name = genome_tuple[NAME_INDEX]
		records.append(record)
	

	records.append(base_record)

	return records


def is_score_segnificant(score, median, std):

	return score > median + 2 * std


def run_multiple_sequence_alignment(records, aln_path, dnd_path):

	SeqIO.write(records, "tmp/emma.fasta", "fasta")

	cmd = EmmaCommandline(sequence="tmp/emma.fasta",
						  outseq=aln_path,
						  dendoutfile=dnd_path)

	print cmd
	stdout, stderr = cmd()
		
	print(stdout + stderr)

def test_run_local_alignment():

	seq_list = get_seq_list("./NC_000913.fna.gz", "./candidate_list.table")

	# seq_list holds tuples of (name, seq_obj)
	seq_file_name = "./tmp/%s.fasta" % seq_list[0][0]
	base_record = SeqRecord(seq_list[0][1], id=seq_file_name, name=seq_file_name)
	SeqIO.write(base_record, seq_file_name, "fasta")

	local_alignment_list = \
		run_local_alignment(seq_file_name, "./genome_list.table")

	print "*" * 100
	print "local alignments results"
	print "*" * 100

	for row in local_alignment_list:
		print row

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
	run_multiple_sequence_alignment(records_for_emma)

if __name__ == "__main__":

	test_run_local_alignment()
