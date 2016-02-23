from db_utils import upload_table
from config_handler import read_config_file_as_table

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
        "-c", "--config_file",
        help="The configuration file for the program")

    parser.add_option(      # customized description; put --help last
        "-h", "--help", action="help",
        help="Show this help message and exit.")

    settings, args = parser.parse_args(argv)

    # check number of arguments, verify values, etc.:
    if args:
        parser.error('program takes no command-line arguments; '
                     '"%s" ignored.' % (args,))

    # further process settings & args if necessary

    return settings, args

def upload_tables(table_list):

	for our_file, table_name in table_list:
		print "*" * 100
		print "Uploading %s with name %s" % (our_file, table_name)
		#upload_table(our_file, table_name)

	print "*" * 100



if __name__ == "__main__":
	settings, args = process_command_line(None)

	try: 
		config_rows = read_config_file_as_table(settings.config_file)

	except 

	table_list = zip([row["file_name"] for row in config_rows],
			 [row["table_name"] for row in config_rows])

	upload_tables(table_list)

	exit(0)
