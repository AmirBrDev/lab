import optparse
import sys
import os

SEQ_FIELD = 4
DELIMITER = " "

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
		'-f' , '--input_file',
		help='Path to file with sequences for logo')

	parser.add_option(
		'-o', '--output_file',
		help='Path to the logo to create.')

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


def run(argv=None):

	try:
		settings, args = process_command_line(argv)

	except BaseException:
		return 1

	if not settings.output_file:
		settings.output_file = "%s_logo.png" % os.path.basename(settings.input_file)

	with open(settings.input_file, "rb") as in_file:

		logo_data_file = "/tmp/%s.fasta" % os.path.basename(settings.input_file)

		with open (logo_data_file, "wb") as out_file:

			for line in in_file:

				if line == "\n":
					continue

				sline = line.replace("\n", "").split(DELIMITER)				

				while "" in sline:
					sline.remove("")

				out_file.write(">\n")
				out_file.write("%s\n" % sline[SEQ_FIELD].replace("T", "U"))

	cmd = "/home/users/amirbar/.local/bin/weblogo " \
		  "-A rna " \
		  "-c classic " \
		  "--resolution 600 " \
		  "--errorbars NO " \
		  "--format PNG " \
		  "--show-xaxis NO" \
		  "< %(logo_data)s  > %(logo_file)s" % {"logo_file": settings.output_file,
		  										"logo_data": logo_data_file}

	os.system(cmd)

	return 0

if __name__ == "__main__":
	exit (run())