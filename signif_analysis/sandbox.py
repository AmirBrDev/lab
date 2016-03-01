import csv
import MySQLdb
import sys
import optparse
from config_handler import ConfigException, read_config_file_as_table

def get_known_table(cursor):

	cursor.execute("SELECT sRNA, target FROM known_interactions")

	row = cursor.fetchone()

	table = []

	while row is not None:

		if row["sRNA"] > row["target"]:
			name = "%s::%s" % (row["target"], row["sRNA"])
		else:
			name = "%s::%s" % (row["sRNA"],row["target"])

		table.append(name)

		row = cursor.fetchone()

	return table


def get_total_interactions(code_name, threshold, table_list, cursor):

	union_statement = ""

	for table_name in table_list[1:]:
		union_statement += """ UNION
	SELECT interactions
	FROM %(table_name)s
	WHERE interactions >= %(threshold)s AND
		  code_name = '%(code_name)s'""" % {"table_name": table_name,
								   			"threshold": threshold,
								   			"code_name": code_name}

	query = """SELECT interactions
	FROM %(table_name)s
	WHERE interactions >= %(threshold)s AND code_name = '%(code_name)s'
	%(union_statement)s""" % {"table_name": table_list[0],
							  "union_statement": union_statement,
							  "threshold": threshold,
							  "code_name": code_name}

	cursor.execute(query)

	row = cursor.fetchone()

	count = 0

	while row is not None:

		count += int(row["interactions"])

		row = cursor.fetchone()

	return count


def get_total_libs(code_name, threshold, table_list, cursor):

	count = 0

	for table_name in table_list:
		query = """SELECT code_name 
		FROM %(table_name)s 
		WHERE interactions >= '%(threshold)s' AND 
		code_name = '%(code_name)s'""" % {"table_name": table_name,
										  "threshold": threshold,
										  "code_name": code_name}

		cursor.execute(query)

		row = cursor.fetchone()

		if row is not None:
			count += 1

	return count



def get_table_code_names(table_name, threshold, cursor):

	cursor.execute("SELECT code_name FROM %s WHERE interactions >= %s" % (table_name, threshold))

	row = cursor.fetchone()

	names = []

	while row is not None:

		names.append(row["code_name"])

		row = cursor.fetchone()

	return list(set(names))


def get_interactions_list(table_name_list, threshold, cursor):
	interactions_list = []

	for table_name in table_name_list:
		interactions_list += get_table_code_names(table_name, threshold, cursor)

	interactions_list = list(set(interactions_list))

	return interactions_list


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

	parser.add_option("-t",
					  "--threshold",
					  help="threshold of minimal reads count allowed")

	parser.add_option(  # customized description; put --help last
						"-h", "--help", action="help",
						help="Show this help message and exit.")

	settings, args = parser.parse_args(argv)

	# check number of arguments, verify values, etc.:
	if args:
		parser.error('program takes no command-line arguments; '
					 '"%s" ignored.' % (args,))

	# further process settings & args if necessary

	return settings, args


def run():
	
	settings, args = process_command_line(None)

	try:
		if settings.threshold:
			threshold = int(settings.threshold)

		else:
			threshold = 15

	except ConfigException as ex:
		print ex
		return 1

	log_threshold = 3
	iron_threshold = 2
	stationary_threshold = 2

	db = MySQLdb.connect(host="localhost", user="amirbar", db="RILseq_seperated_annotations")
	cursor = db.cursor(MySQLdb.cursors.DictCursor)

	known_interactions_list = get_known_table(cursor)

	print "threshold: %d" % threshold

	print ""
	print "*" * 50
	print "Unified counts"
	print "*" * 50

	# unified_log_names = get_table_code_names("unified_all_interactions_log", threshold, cursor)
	# unified_iron_names = get_table_code_names("unified_all_interactions_iron", threshold, cursor)
	# unified_stationary_names = get_table_code_names("unified_all_interactions_stationary", threshold, cursor)

	unified_log_names = get_table_code_names("unified_signif_log", threshold, cursor)
	unified_iron_names = get_table_code_names("unified_signif_iron", threshold, cursor)
	unified_stationary_names = get_table_code_names("unified_signif_stationary", threshold, cursor)

	total_interactions_names = list(set(unified_log_names + 
										unified_iron_names + 
										unified_stationary_names))

	total_interaction_count = len(total_interactions_names)
	known_count = sum([1 for interaction in total_interactions_names if interaction in known_interactions_list])

	print "total interactions=", total_interaction_count
	print "known without lib requirement=", known_count
	print "unknown=", total_interaction_count - known_count

	print "\n"
	print "*" * 50
	print "Generating Table"
	print "*" * 50

	# with open("output/names_to_counts_by_libs_%d" % threshold, "wb") as fl:
	with open("output/signif_names_to_counts_by_libs_%d" % threshold, "wb") as fl:


		writer = csv.writer(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

		log_repr_tables_names = ['log_101',
								 'log_102',
								 'log_103',
								 'log_104',
								 'log_108',
								 'log_109']

		stationary_repr_tables_names = ['stationary_209',
										'stationary_210',
										'stationary_312']

		iron_repr_tables_names = ['iron_207',
								  'iron_208',
								  'iron_305']

		table_repr_names = log_repr_tables_names + \
						   iron_repr_tables_names + \
						   stationary_repr_tables_names

		# table_repr_name_to_db_name = {"log_101": "all_interactions_log_101",
		# 							  "log_102": "all_interactions_log_102",
		# 							  "log_103": "all_interactions_log_103",
		# 							  "log_104": "all_interactions_log_104",
		# 							  "log_108": "all_interactions_log_108",
		# 							  "log_109": "all_interactions_log_109",
		# 							  "iron_207": "all_interactions_iron_207",
		# 							  "iron_208": "all_interactions_iron_208",
		# 							  "iron_305": "all_interactions_iron_305",
		# 							  "stationary_209": "all_interactions_stationary_209",
		# 							  "stationary_210": "all_interactions_stationary_210",
		# 							  "stationary_312": "all_interactions_stationary_312"}

		table_repr_name_to_db_name = {"log_101": "signif_log_101",
									  "log_102": "signif_log_102",
									  "log_103": "signif_log_103",
									  "log_104": "signif_log_104",
									  "log_108": "signif_log_108",
									  "log_109": "signif_log_109",
									  "iron_207": "signif_iron_207",
									  "iron_208": "signif_iron_208",
									  "iron_305": "signif_iron_305",
									  "stationary_209": "signif_stationary_209",
									  "stationary_210": "signif_stationary_210",
									  "stationary_312": "signif_stationary_312"}

		writer.writerow(["interaction"] + table_repr_names + ["log_libs", "iron_libs", "stationary_libs"])

		row_count = 0
		unique_with_lib_requirement = 0

		for interaction_name in total_interactions_names:

			row = [interaction_name]

			log_libs = get_total_libs(interaction_name, threshold, [table_repr_name_to_db_name[repr_name] for repr_name in log_repr_tables_names], cursor)
			iron_libs = get_total_libs(interaction_name, threshold, [table_repr_name_to_db_name[repr_name] for repr_name in iron_repr_tables_names], cursor) 
			stationary_libs = get_total_libs(interaction_name, threshold, [table_repr_name_to_db_name[repr_name] for repr_name in stationary_repr_tables_names], cursor)

			if log_libs < log_threshold and iron_libs < iron_threshold and stationary_libs < stationary_threshold:
				continue

			unique_with_lib_requirement += 1


			if interaction_name not in known_interactions_list:
				continue

			for repr_name in table_repr_names:
				row.append(get_total_interactions(interaction_name, threshold, [table_repr_name_to_db_name[repr_name]], cursor))

			row.append(log_libs)
			row.append(iron_libs)
			row.append(stationary_libs)

			writer.writerow(row)
			row_count += 1

		print "unique names above threshold: %d" % len(total_interactions_names)
		print "unique names above threshold and libs requirement: %d" % unique_with_lib_requirement
		print "known above threshold and libs rquirement: %d" % row_count

	print "Done!"

	return 0

if __name__ == "__main__":
	exit(run())
