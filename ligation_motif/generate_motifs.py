__author__ = 'amirbar'

import os
from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet import Gapped, IUPAC
import matplotlib.pyplot as plt
import optparse
import sys
import os

GAP = "-"
ALPHABET = Gapped(IUPAC.unambiguous_dna)

def process_command_line(argv):
	"""
	Return a 2-tuple: (settings object, args list).
	`argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
	"""
	if argv is None:
		argv = sys.argv[1:]

	# initialize the parser object:
	parser = optparse.OptionParser(
		formatter=optparse.TitledHelpFormatter(width=100),
		add_help_option=None)

	parser.add_option(
		"-r", "--reads_fusion",
		help="the output file from Asaf's script 'detect_fusion_point.py'")

	parser.add_option(
		"-l", "--logo_file",
		help="file name for the output logo file")

	parser.add_option(
		"-s", "--score_diagram",
		help="file name for the output score diagram")

	parser.add_option(
		"-c", "--count_diagram",
		help="file name for the output count diagram")

	parser.add_option(
		"-d", "--distance",
		help="maximal distance from ligation point", default=10, type="int")

	parser.add_option(	  # customized description; put --help last
		"-h", "--help", action="help",
		help="Show this help message and exit.")

	settings, args = parser.parse_args(argv)

	# check number of arguments, verify values, etc.:
	if args:
		parser.error('program takes no command-line arguments; '
					 '"%s" ignored.' % (args,))

	# further process settings & args if necessary

	return settings, args

def format_input(chimera_frags_list):

	result = []

	for chimera in chimera_frags_list:
		first_frag, second_frag = chimera
		result.append((Seq(first_frag, ALPHABET), Seq(second_frag, ALPHABET)))

	return  result


def format_file_input(path_to_file):

	FRAG_1_SEQ_INDEX = 6
	FRAG_2_SEQ_INDEX = 8
	FRAG_1_SCORE_START = 9

	fragment_as_seq_list = []
	positive_pos_score_list = []
	negative_pos_score_list = []

	with open(path_to_file, "rb") as fl:

		for line in fl:
			row = line.replace("\n", "").split("\t")
			FRAG_2_SCORE_START = FRAG_1_SCORE_START + len(row[FRAG_1_SEQ_INDEX])

			if len(row[FRAG_1_SEQ_INDEX]) + len(row[FRAG_2_SEQ_INDEX]) != len(row[FRAG_1_SCORE_START:]):
				print "-" * 100
				print "skipped row"
				print row[FRAG_1_SEQ_INDEX]
				print row[FRAG_2_SEQ_INDEX]
				print "-" * 100
				continue

			# Get the fragments
			fragment_as_seq_list.append((Seq(row[FRAG_1_SEQ_INDEX], ALPHABET), Seq(row[FRAG_2_SEQ_INDEX], ALPHABET)))

			# Get negative scores
			# print "*" * 50
			# print "row length:", len(row)
			# print "start-end:", FRAG_1_SCORE_START, FRAG_2_SCORE_START - 1
			score = [float(row[index]) for index in range(FRAG_2_SCORE_START - 1, FRAG_1_SCORE_START - 1, -1)]
			negative_pos_score_list.append(score)

			# Get positive scores
			# print FRAG_2_SCORE_START, len(row[FRAG_2_SEQ_INDEX]), len(row)
			score = [float(row[index]) for index in range(FRAG_2_SCORE_START, FRAG_2_SCORE_START + len(row[FRAG_2_SEQ_INDEX]))]

			positive_pos_score_list.append(score)


	return fragment_as_seq_list, negative_pos_score_list, positive_pos_score_list


def get_max_len_in_chimera(fragment_as_seq_list):

	first_frag_list, second_frag_list = zip(*fragment_as_seq_list)

	return len(max(first_frag_list, key=len)), len(max(second_frag_list, key=len))


def pad_sequences(chimera_frags_list):

	max_first, max_second = get_max_len_in_chimera(chimera_frags_list)

	result = []

	for chimera in chimera_frags_list:
		first_frag, second_frag = chimera

		first_frag = (max_first - len(first_frag)) * GAP + first_frag
		second_frag += (max_second - len(second_frag)) * GAP

		result.append((first_frag, second_frag))

	return result


def merge_frgaments(padded_fragments_list):

	result = []

	for chimera in padded_fragments_list:
		first_frag, second_frag = chimera
		result.append(Seq(str(first_frag + second_frag), ALPHABET))

	return result


def create_sequnce_file(filename, sequence_list, lower_bound, upper_bound, ligation_index):

	with open(filename, "wb") as fl:

		for sequence in sequence_list:
			fl.write(">\n")
			fl.write(str(sequence)[ligation_index + lower_bound:ligation_index] + \
					 "-" + \
					 str(sequence)[ligation_index:ligation_index+upper_bound + 1] + \
					 "\n")


def get_position_statistics(fragment_as_seq_list, negative_pos_score_list, positive_pos_score_list):

	max_first, max_second = get_max_len_in_chimera(fragment_as_seq_list)

	pos_sum_dct = {index: 0 for index in range(-max_first, 0)}
	pos_sum_dct.update({index: 0 for index in range(1, max_second + 1)})
	pos_count_dct = {index: 0 for index in range(-max_first, 0)}
	pos_count_dct.update({index: 0 for index in range(1, max_second + 1)})

	# Sum over the negatives
	for score_list in negative_pos_score_list:
		for index, score in enumerate(score_list):

			pos_sum_dct[-(index + 1)] += score
			pos_count_dct[-(index + 1)] += 1

	# Sum score over the positives
	for score_list in positive_pos_score_list:
		for index, score in enumerate(score_list):

			pos_sum_dct[index + 1] += score
			pos_count_dct[index + 1] += 1

	return pos_sum_dct, pos_count_dct


def print_settings(settings):

	print "-" * 100
	print "| settings:"
	print "-" * 100

	for key, value in settings.items():
		print "| %s = %s" % (key, value)

	print "-" * 100

def run(argv=None):

	settings, args = process_command_line(None)

	fragment_as_seq_list, negative_pos_score_list, positive_pos_score_list = \
		format_file_input(settings.reads_fusion)

	# Score diagram
	pos_sum_dct, pos_count_dct = get_position_statistics(fragment_as_seq_list,
														 negative_pos_score_list,
														 positive_pos_score_list)

	stats_dict = {position: float(pos_sum_dct[position]) / float(pos_count_dct[position])
				  for position in pos_sum_dct.keys()}


	ligation_index = abs(min(stats_dict.keys()))
	lower_bound = min(stats_dict.keys())
	upper_bound = max(stats_dict.keys())

	if settings.distance:
		lower_bound = max(lower_bound, -settings.distance)
		upper_bound = min(upper_bound, settings.distance)

	# Generate average base score plot
	points = [[i, stats_dict[i]] for i in range (lower_bound, 0)] + \
			 [[i, stats_dict[i]] for i in range (1, upper_bound + 1)]

	plt.xlim(lower_bound - 1, upper_bound + 1)
	plt.ylim(0, 1)
	plt.hold(True)

	for pt in points:
		plt.plot([pt[0], pt[0]], [0,pt[1]], "b")

	if settings.score_diagram != None:
		plt.savefig(settings.score_diagram)

	# Generate count per position plot
	points = [[i, pos_count_dct[i]] for i in range (lower_bound, 0)] + \
			 [[i, pos_count_dct[i]] for i in range (1, upper_bound + 1)]

	plt.xlim(lower_bound - 1, upper_bound + 1)
	plt.ylim(0, max(pos_count_dct.values()))
	plt.hold(True)

	for pt in points:
		plt.plot([pt[0], pt[0]], [0,pt[1]], "r")

	if settings.count_diagram != None:
		plt.savefig(settings.count_diagram)

	# Generate logo
	if settings.logo_file != None:
		padded_fragments_list = pad_sequences(fragment_as_seq_list)

		merged_fragments_list = merge_frgaments(padded_fragments_list)

		chimera_motifs = motifs.create(merged_fragments_list, ALPHABET)

		print chimera_motifs
		print chimera_motifs.counts

		create_sequnce_file("logo_data", merged_fragments_list, lower_bound, upper_bound, ligation_index)

		os.system("/home/users/amirbar/.local/bin/weblogo -A dna -c classic --resolution 600 --errorbars NO -i %(start_index)s --format PNG < logo_data  > %(logo_file)s" % \
				  {"logo_file": settings.logo_file,
				   "start_index": str(lower_bound)})
run()