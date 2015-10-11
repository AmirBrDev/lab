import gzip
import Bio.SeqIO
from Table import Table

def get_seq(genome_file, start, end, strand, file_format="fasta"):

	# Retrieve the sequence
	with gzip.open(genome_file, "rb") as fl:

		genome_seq = Bio.SeqIO.read(fl, file_format).seq

	# Cut it to the specified coordinate
	if strand == "+":
		return genome_seq[start - 1 : end]

	else:
		return genome_seq[start - 1 : end].reverse_complement()


def test_get_seq():

	expected = "AGCTT"
	seq = get_seq("./NC_000913.fna.gz", 1, 5, "+")

	print expected == str(seq)

	expected = "GCTT"
	seq = get_seq("./NC_000913.fna.gz", 2, 5, "+")

	print expected == str(seq)

	expected = "AAGC"
	seq = get_seq("./NC_000913.fna.gz", 2, 5, "-")

	print expected == str(seq)


def get_seq_list(genome_file, seq_list_file):

	table = Table(seq_list_file)

	return [(entry["name"], get_seq(genome_file, 
									int(entry["start"]), 
									int(entry["end"]),
									entry["strand"])) for entry in table]


def test_get_seq_list():

	seq_list = get_seq_list("./NC_000913.fna.gz", "./candidate_list.table")
	print seq_list, len(seq_list[0])

	print len(seq_list) == 1

if __name__ == "__main__":

	test_get_seq_list()