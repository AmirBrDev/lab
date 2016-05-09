__author__ = 'amirbar'

import optparse
import sys
import csv
from math import log, floor
import matplotlib.pyplot as plt

LENGTH_FIELD = 1
COUNT_FIELD = 2

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
		'-i' , '--input_file',
		help='Path to info file (output of Yaels perl script) containing the ')

	parser.add_option(
		'-o' , '--output_file',
		help='Path of output file')

	parser.add_option(
		'-e' , '--log_exponent',
		default=10,
		type=int,
		help='the exponent to use for creating buckets')

	parser.add_option(	  # customized description; put --help last
		'-h', '--help', action='help',
		help='Show this help message and exit.')

	settings, args = parser.parse_args(argv)

	# check number of arguments, verify values, etc.:
	if args:
		parser.error('program takes no command-line arguments; '
					 '"%s" ignored.' % (args,))

	# further process settings & args if necessary
	if not settings.input_file:
		parser.error("Missing missing name/path for input file.")

	if not settings.output_file:
		parser.error("Missing missing name/path for output_file file.")

	return settings, args


def run(argv=None):

	try:
		settings, args = process_command_line(argv)

	except BaseException:
		return 1

	points = {}

	# Get dictionary of sizes to counts of size
	with open(settings.input_file, "rb") as fl:

		# ignore first 5 rows
		for i in range(5):
			fl.readline()

		reader = csv.reader(fl, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)

		for row in reader:
			size = float(row[1])
			count = int(row[2])

			# calculate the bucket
			bucket_id = int(floor(log(size) / log(settings.log_exponent)))


			# Create new bucket if not exists
			if bucket_id not in points.keys():
				points[bucket_id] = 0

			points[bucket_id] += count

	# Generate the plot
	buckets_list = []
	counts_list = []
	
	# Convert to graph points - sizes as usual, counts are in logaritmic scale
	for bucket, count in points.items():
		buckets_list.append(bucket)
		counts_list.append(log(count) / log(settings.log_exponent))

	plt.scatter(buckets_list,counts_list)
	plt.plot(buckets_list, counts_list)
	plt.savefig(settings.output_file)



if __name__ == "__main__":
	exit(run())