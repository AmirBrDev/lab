from Table import Table
from sequence_extractor import get_seq_list
from Bio.Emboss.Applications import NeedleCommandline, WaterCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
#from Bio import AlignIO
import AlignIO
from EmbossCommands import EmmaCommandline, SuperMatcherCommandline
import gzip
import os
import numpy
import pylab

def run_local_alignment(seq_file, genome_list, workdir):

	genome_table = Table(genome_list)
	result = []

	for row in genome_table:
		genome_name = row["genome_file"]
		print genome_name

		# Find the genome file according to the genome directory
		for file_name in os.listdir("./microbial/%s" % genome_name):
			if file_name.endswith(".fna.gz"):
				genome_file = "./microbial/%s/%s" % (genome_name, file_name)

		# Run on both strands
		result.append(directed_local_alignment(genome_file,
											   seq_file,
											   genome_name,
											   True,
											   workdir))

		result.append(directed_local_alignment(genome_file,
											   seq_file,
											   genome_name,
											   False,
											   workdir))

	os.remove("./tmp/curr_genome.fasta")

	return result

def directed_local_alignment(genome_file, seq_file, genome_name, is_positive, workdir):

	if is_positive:
		name_extension = "positive"

	else:
		name_extension = "negative"

	with gzip.open(genome_file, "rb") as fl:
		records = SeqIO.read(fl, "fasta")

		if not is_positive:
			records = records.reverse_complement(id=True, 
												 name=True, 
												 description=True, 
												 features=True, 
												 annotations=True, 
												 letter_annotations=True, 
												 dbxrefs=True)

		SeqIO.write(records, "%s/curr_genome.fasta" % workdir, "fasta")

	alignment_file = "%s/aln/%s_%s.aln" % (workdir, genome_name, name_extension)

	# Build the matching command
	cmd = SuperMatcherCommandline(asequence=seq_file, 
								  bsequence="%s/curr_genome.fasta" % workdir, 
								  gapopen=10, 
								  gapextend=0.5,
								  outfile=alignment_file)

	# Excecute the command
	stdout, stderr = cmd()
	
	#print(stdout + stderr)

	# Get the score for the genome and the matching alignment
	align_seq = list(AlignIO.read(alignment_file, "amir_emboss"))[1]
	align_seq._set_seq(align_seq.seq.ungap("-"))  # remove the gaps
	score = extract_score(alignment_file)

	return ("_".join([genome_name, name_extension]), score, align_seq)


def extract_score(file_name):

	with open(file_name, "rb") as fl:
		
		line = fl.readline()

		while (not line.startswith("# Score: ")):
			line = fl.readline()

	return float(line.replace("# Score: ", ""))


def generate_score_histogram(path, local_alignment_list):
	SCORE_INDEX = 1

	score_list = \
		[genome_tuple[SCORE_INDEX] for genome_tuple in local_alignment_list]

	pylab.hist(numpy.array(score_list))
	pylab.savefig(path)
	pylab.close()


def prepare_for_multiple_aignment(local_alignment_list, 
								  proj_id_to_tax_id_dict,
								  secondary_proj_id_to_tax_id_dict):
	NAME_INDEX = 0
	SCORE_INDEX = 1
	RECORD_INDEX = 2

	score_list = \
		[genome_tuple[SCORE_INDEX] for genome_tuple in local_alignment_list]

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
		record.id = \
			get_id_by_name(genome_tuple[NAME_INDEX], proj_id_to_tax_id_dict)

		if record.id is None:
			record.id = \
				get_id_by_name(genome_tuple[NAME_INDEX], secondary_proj_id_to_tax_id_dict)

		if record.id is None:
			raise BaseException("Unmatched taxonomy id to given name %s" % \
				genome_tuple[NAME_INDEX])


		record.name = genome_tuple[NAME_INDEX]
		records.append(record)

	return records


def build_project_id_to_tax_id_dictionary(file_path):

	with gzip.open(file_path, "rb") as fl:

		# Remove two header lines
		fl.readline()
		fl.readline()

		result = {}

		# Build the table
		for line in fl.readlines():

			args = line.split("\t")
			result[args[0]] = args[2]

	return result

def build_secondary_project_id_to_tax_id_dictionary(file_path):

	with gzip.open(file_path, "rb") as fl:

		# Remove two header lines
		fl.readline()
		fl.readline()

		result = {}

		# Build the table
		for line in fl.readlines():

			args = line.split("\t")
			result[args[1]] = args[2]

	return result

def get_id_by_name(name, proj_id_to_tax_id_dict):
	appendix = name[name.rfind("_") + 1:]
	name = name.replace("_positive", "").replace("_negative", "")
	project_id = name[name.find("uid") + 3:]

	try:
		if appendix == "positive":
			return proj_id_to_tax_id_dict[project_id]
		else:
			return "%s_neg" % proj_id_to_tax_id_dict[project_id]

	except KeyError:
		return None


def is_score_segnificant(score, median, std):

	return score > median + 2 * std


def run_multiple_sequence_alignment(records, aln_path, dnd_path, workdir):

	SeqIO.write(records, "%s/emma.fasta" % workdir, "fasta")

	cmd = EmmaCommandline(sequence="%s/emma.fasta" % workdir,
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
