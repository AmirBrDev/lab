from samba.dcerpc.nbt import db_change_info

__author__ = 'amirbar'

import optparse
import sys
import MySQLdb
import csv

DELIMITER = "\t"


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
		'-k', '--known_table',
		help='Path to file with sequences for logo')

	parser.add_option(
		'-d', '--database',
		default='article_refactor_24_3_2016',
		help='database to work with')

	parser.add_option(
		'-o', '--output_file',
		help='Path to the program output file')

	parser.add_option(	  # customized description; put --help last
		'-h', '--help', action='help',
		help='Show this help message and exit.')

	settings, args = parser.parse_args(argv)

	# check number of arguments, verify values, etc.:
	if args:
		parser.error('program takes no command-line arguments; '
					 '"%s" ignored.' % (args,))

	# further process settings & args if necessary
	if not settings.known_table:
		parser.error("Missing known interactions table.")

	if not settings.output_file:
		parser.error("Missing missing name/path for output file.")

	return settings, args


def get_interactions(file_path):

	code_names = []
	is_found = []

	# Go over the known file - table s3
	with open(file_path, "rb") as fl:

		line = fl.readline()

		for line in fl:
			sline = line.strip().split(DELIMITER)

			first = sline[0]
			second = sline[1]

			# choose correct order and append to interaction names
			if first > second:
				code_names.append("%s::%s" % (second, first))

			else:
				code_names.append("%s::%s" % (first, second))

			# check whether the interaction were found
			is_found.append(sline[5] != "0")

	return code_names, is_found


def get_srna_names(cursor, db_tables):

	union_statement = ""

	for table_name in db_tables[1:]:
		union_statement += """ UNION
	SELECT rna1_name, first_type
	FROM %(table_name)s
	WHERE first_type='srna'
	UNION
	SELECT rna2_name, second_type
	FROM %(table_name)s
	WHERE second_type='srna'""" % {"table_name": table_name}

	query = """SELECT rna1_name, first_type
	FROM %(table_name)s
	WHERE first_type='srna'
	UNION
	SELECT rna2_name, second_type
	FROM %(table_name)s
	WHERE second_type='srna'
	%(union_statement)s""" % {"table_name": db_tables[0],
							  "union_statement": union_statement}

	cursor.execute(query)

	row = cursor.fetchone()

	names = []

	while row is not None:
		if row["rna1_name"] not in names:
			names.append(row["rna1_name"])

		row = cursor.fetchone()

	return names


def get_interaction_count(cursor, db_tables, code_name):

	union_statement = ""

	for table_name in db_tables[1:]:
		union_statement += """ UNION
	SELECT interactions
	FROM %(table_name)s
	WHERE code='%(code_name)s'""" % {"table_name": table_name,
									 "code_name": code_name}

	query = """SELECT interactions
	FROM %(table_name)s
	WHERE code='%(code_name)s'
	%(union_statement)s""" % {"table_name": db_tables[0],
							  "code_name": code_name,
							  "union_statement": union_statement}

	cursor.execute(query)

	row = cursor.fetchone()

	total_interactions = 0

	while row is not None:
		total_interactions += int(row["interactions"])

		row = cursor.fetchone()

	return total_interactions


def get_interaction_fold_energy(cursor, db_tables, code_name):

	union_statement = ""

	for table_name in db_tables[1:]:
		union_statement += """ UNION
	SELECT free_energy_of_hybridization
	FROM %(table_name)s
	WHERE code='%(code_name)s'""" % {"table_name": table_name,
									 "code_name": code_name}

	query = """SELECT free_energy_of_hybridization
	FROM %(table_name)s
	WHERE code='%(code_name)s'
	%(union_statement)s""" % {"table_name": db_tables[0],
							  "code_name": code_name,
							  "union_statement": union_statement}

	cursor.execute(query)

	row = cursor.fetchone()

	fold_energy = []

	while row is not None:

		fold_energy.append(float(row["free_energy_of_hybridization"]))

		row = cursor.fetchone()

	if len(fold_energy) > 0:
		return min(fold_energy)

	else:
		return "-"


def get_target_total_reads(cursor, condition_table, code_name, srna_name):

	query = """SELECT rna1_name, rna2_name, total_rna_reads1, total_rna_reads2
	FROM %(table_name)s
	WHERE code='%(code_name)s'""" % {"table_name": condition_table,
							  		 "code_name": code_name}

	cursor.execute(query)

	row = cursor.fetchone()

	while row is not None:

		if row["rna1_name"] == srna_name:
			return row["total_rna_reads2"]

		elif row["rna2_name"] == srna_name:
			return row["total_rna_reads1"]

		row = cursor.fetchone()

	return "-"


def get_target_normalized_ip_interactions(cursor, condition_table, code_name, srna_name):

	query = """SELECT rna1_name, rna2_name, lib_norm_ip_rna1, lib_norm_ip_rna2
	FROM %(table_name)s
	WHERE code='%(code_name)s'""" % {"table_name": condition_table,
							  		 "code_name": code_name}

	cursor.execute(query)

	row = cursor.fetchone()

	while row is not None:

		if row["rna1_name"] == srna_name:
			return row["lib_norm_ip_rna2"]

		elif row["rna2_name"] == srna_name:
			return row["lib_norm_ip_rna1"]

		row = cursor.fetchone()

	return "-"


def run(argv=None):

	try:
		settings, args = process_command_line(argv)

	except BaseException:
		return 1


	log_table = "signif_chimeras_of_log_phase_cl"
	stat_table = "signif_chimeras_of_stationary_cl"
	iron_table = "signif_chimeras_of_iron_limitation_cl"

	db_tables = [log_table,
				 iron_table,
				 stat_table]

	all_interactions_log_table = "all_interactions_with_fold_log"
	all_interactions_stat_table = "all_interactions_with_fold_stat"
	all_interactions_iron_table = "all_interactions_with_fold_iron"

	all_interactions_db_tables = [all_interactions_log_table,
				 				  all_interactions_iron_table,
				 				  all_interactions_stat_table]

	db = MySQLdb.connect(host="localhost", user="amirbar", db=settings.database)
	cursor = db.cursor(MySQLdb.cursors.DictCursor)

	srna_names = get_srna_names(cursor, db_tables)
	interaction_names, status_list = get_interactions(settings.known_table)

	interaction_to_found_status = {interaction: status for interaction, status in zip(interaction_names, status_list)}

	with open(settings.output_file, "wb") as fl:

		writer = csv.writer(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer.writerow(["code",
						 "is_found",
						 "interaction_count",
						 "free_energy",
						 "log_target_total_reads",
						 "stationary_target_total_reads",
						 "iron_target_total_reads",
						 "log_norm_ip",
						 "stationary_norm_ip",
						 "iron_norm_ip"])

		# Go over the interactions and generate the table
		for interaction in interaction_names:

			first_name, second_name = interaction.lower().split("::")

			if first_name in srna_names:
				srna = first_name

			elif second_name in srna_names:
				srna = second_name

			else:
				print "skipped interaction between two sRNAs - %s" % interaction
				continue

			# row = [interaction,
			# 	   get_interaction_count(cursor, db_tables, interaction),
			# 	   get_interaction_fold_energy(cursor, db_tables, interaction),
			# 	   get_target_total_reads(cursor, log_table, interaction, srna),
			# 	   get_target_total_reads(cursor, stat_table, interaction, srna),
			# 	   get_target_total_reads(cursor, iron_table, interaction, srna),
			# 	   get_target_normalized_ip_interactions(cursor, log_table, interaction, srna),
			# 	   get_target_normalized_ip_interactions(cursor, stat_table, interaction, srna),
			# 	   get_target_normalized_ip_interactions(cursor, iron_table, interaction, srna)]
			row = [interaction,
				   interaction_to_found_status[interaction],
				   get_interaction_count(cursor, all_interactions_db_tables, interaction),
				   get_interaction_fold_energy(cursor, all_interactions_db_tables, interaction),
				   get_target_total_reads(cursor, all_interactions_log_table, interaction, srna),
				   get_target_total_reads(cursor, all_interactions_stat_table, interaction, srna),
				   get_target_total_reads(cursor, all_interactions_iron_table, interaction, srna),
				   get_target_normalized_ip_interactions(cursor, all_interactions_log_table, interaction, srna),
				   get_target_normalized_ip_interactions(cursor, all_interactions_stat_table, interaction, srna),
				   get_target_normalized_ip_interactions(cursor, all_interactions_iron_table, interaction, srna)]
			

			if row[2] > 0 and row[4] == "-" and row[5] == "-" and row[6] == "-":
				print "[Warning] %s" % str(row)

			if row[2] == 0:
				print "missed %s" % row[0]
				# print "retrieving data from other tables..."

				# row = [interaction,
				#	get_interaction_count(cursor, all_interactions_db_tables, interaction),
				#	get_interaction_fold_energy(cursor, all_interactions_db_tables, interaction),
				#	get_target_total_reads(cursor, all_interactions_log_table, interaction, srna),
				#	get_target_total_reads(cursor, all_interactions_stat_table, interaction, srna),
				#	get_target_total_reads(cursor, all_interactions_iron_table, interaction, srna),
				#	get_target_normalized_ip_interactions(cursor, all_interactions_log_table, interaction, srna),
				#	get_target_normalized_ip_interactions(cursor, all_interactions_stat_table, interaction, srna),
				#	get_target_normalized_ip_interactions(cursor, all_interactions_iron_table, interaction, srna)]
				# print row

			writer.writerow(row)

	return 0

if __name__ == "__main__":
	exit(run())
