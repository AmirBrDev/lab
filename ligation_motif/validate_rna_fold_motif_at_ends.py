import optparse
import sys
from RNA.structure import computeACC
import numpy as np

def process_command_line(argv):
	"""
	Return a 2-tuple: (settings object, args list).
	`argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
	"""
	if argv is None:
		argv = sys.argv[1:]

	# initialize the parser object:
	parser = optparse.OptionParser(
		formatter=optparse.TitledHelpFormatter(width=78),
		add_help_option=None)

	parser.add_option(
		'-r', '--reads_file',
		help='The reads_fusion.txt file - the output of detect_fusion_point.py script')

	parser.add_option(	  # customized description; put --help last
		'-h', '--help', action='help',
		help='Show this help message and exit.')

	settings, args = parser.parse_args(argv)

	# check number of arguments, verify values, etc.:
	if args:
		parser.error('program takes no command-line arguments; '
					 '"%s" ignored.' % (args,))

	# further process settings & args if necessary

	return settings, args

def generate_permutation(perm_length):

	letters = ['A', 'C', 'G', 'T']

	curr_perms = [""]

	for j in range(perm_length):
		curr_perms = [perm + i for perm in curr_perms for i in letters]

	return curr_perms

def run(argv=None):

	try:
		settings, args = process_command_line(argv)

	except BaseException:
		return 1

	FRAG_1_SEQ_INDEX = 6
	FRAG_2_SEQ_INDEX = 8
	FRAG_1_SCORE_START = 9

	positive_pos_score_list = []
	negative_pos_score_list = []

	# Get the reads before and after ligation point
	with open(settings.reads_file, "rb") as fl:

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

			negative_pos_score_list.append(row[FRAG_1_SEQ_INDEX])
			positive_pos_score_list.append(row[FRAG_2_SEQ_INDEX])

	NUM_OF_BASES_TO_EXCHANGE = 4

	perms = generate_permutation(NUM_OF_BASES_TO_EXCHANGE)

	scores_per_fragment = []

	i = 1

	# Calculate the score for each possible end for the fragments before
	# the ligation point
	for frag in negative_pos_score_list:

		base_seq = frag[:-NUM_OF_BASES_TO_EXCHANGE]

		for extension in perms:
			scores_per_fragment.append(computeACC(base_seq + extension, winlen=1))
			row = [str(i), base_seq + extension] + [str(val) for val in scores_per_fragment[-1]]
			print "\t".join(row)
			i += 1

	# calculate statistics for last postions
	postions_scores = [[] for i in range(NUM_OF_BASES_TO_EXCHANGE)] 

	for score_list in scores_per_fragment:
		for i in range(NUM_OF_BASES_TO_EXCHANGE):
			postions_scores[i].append(score_list[-NUM_OF_BASES_TO_EXCHANGE + i])

	for i in range(NUM_OF_BASES_TO_EXCHANGE):
		print "position:"
		print "-" * 50
		print "mean:", np.mean(postions_scores[i])
		print "median:", np.median(postions_scores[i])
		print "max:", max(postions_scores[i])
		print "min:", min(postions_scores[i])
		print "\n"

	print "*" * 100
	print "*" * 100
	print "finished negatives, starting positives"
	print "*" * 100
	print "*" * 100

	# Calculate the score for each possible end for the fragments before
	# the ligation point
	for frag in positive_pos_score_list:

		base_seq = frag[:-NUM_OF_BASES_TO_EXCHANGE]

		for extension in perms:
			scores_per_fragment.append(computeACC(base_seq + extension, winlen=1))
			row = [str(i), base_seq + extension] + [str(val) for val in scores_per_fragment[-1]]
			print "\t".join(row)
			i += 1

	# calculate statistics for first postions
	postions_scores = [[] for i in range(NUM_OF_BASES_TO_EXCHANGE)] 

	for score_list in scores_per_fragment:
		for i in range(NUM_OF_BASES_TO_EXCHANGE):
			postions_scores[i].append(score_list[i])


	for i in range(NUM_OF_BASES_TO_EXCHANGE):
		print "position:"
		print "-" * 50
		print "mean:", np.mean(postions_scores[i])
		print "median:", np.median(postions_scores[i])
		print "max:", max(postions_scores[i])
		print "min:", min(postions_scores[i])
		print "\n"

	return 0

if __name__ == "__main__":
	exit (run())