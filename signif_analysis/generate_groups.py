from config_handler import ConfigException
from config_handler import read_config_file_as_table
import sys
import optparse
import csv
import itertools
from file_handler import FileProcessor, RowHandler

SIGNIF_TYPE = "signif"
INTERACTION_TYPE_KEY = "interaction_type"
CONDITION_KEY = "condition_name"
FILE_PATH_KEY = "file_path"
CODE_NAME_INDEX = -1
LIBS_INDEX = -2
INTERACTIONS_INDEX = 17 - 1

DIV_OPTION_SIGNIF = "signif"
DIV_OPTION_KNOWN_TARGETS = "known"
DIV_OPTION_BINDING = "no_binding"
DIV_OPTION_QUESTIONABLE = "questionable"

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

	parser.add_option("-c",
					  "--config_file",
					  help="path for the configuration file for the program")

	parser.add_option("-u",
					  "--unified_config_file",
					  help="path for the unified configuration file for the program")

	parser.add_option("-k",
					  "--known_targets_file",
					  help="path for the the file containing the known targets")

	parser.add_option("-b",
					  "--no_binding_file",
					  help="path for the the file containing the binding pairs")

	parser.add_option("-s",
					  "--signif",
					  action="store_true",
					  help="whether to add a group category for signif")

	parser.add_option("-q",
					  "--questionable_file",
					  help="path for the the file containing the questionable targets")

	parser.add_option(  # customized description; put --help last
						"-h", "--help", action="help",
						help="Show this help message and exit.")

	settings, args = parser.parse_args(argv)

	# check number of arguments, verify values, etc.:
	if args:
		parser.error('program takes no command-line arguments; '
					 '"%s" ignored.' % (args,))

	if not settings.config_file:
		parser.error('Missing required argument config_file.')

	if settings.unified_config_file and not settings.known_targets_file:
		parser.error('unified config file requires known target file to be given.')

	# further process settings & args if necessary

	return settings, args


def get_condition_files(config_rows):
	dct = {}

	# Go over the conditions and load
	for row in config_rows:

		# Assign if not signif interactions
		if row[INTERACTION_TYPE_KEY] != SIGNIF_TYPE:

			# Add new interaction type
			if not dct.has_key(row[INTERACTION_TYPE_KEY]):
				dct[row[INTERACTION_TYPE_KEY]] = {}

			# Add new condition if not exists
			if not dct[row[INTERACTION_TYPE_KEY]].has_key(row[CONDITION_KEY]):
				dct[row[INTERACTION_TYPE_KEY]][row[CONDITION_KEY]] = []

			dct[row[INTERACTION_TYPE_KEY]][row[CONDITION_KEY]].append(row[FILE_PATH_KEY])

	return dct


def generate_signif_per_conditions_dictionary(config_rows):
	dct = {}

	# Go over the conditions and load
	for row in config_rows:

		# Assign only signif interactions
		if row[INTERACTION_TYPE_KEY] == SIGNIF_TYPE:

			# Add new condition if not exists
			if not dct.has_key(row[CONDITION_KEY]):
				dct[row[CONDITION_KEY]] = []

			# Read the file and add the signif name according to the code
			with open(row[FILE_PATH_KEY], "rb") as fl:

				reader = csv.reader(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

				dct[row[CONDITION_KEY]] += [file_row[CODE_NAME_INDEX] for file_row in reader]

	# Make sure that there are no duplicate names
	for key, val in dct.items():
		dct[key] = list(set(val))

	return dct


def generate_max_lib_per_condition_dictionary(config_rows):
	dct = {}

	# Go over the conditions and load
	for row in config_rows:

		# Assign only signif interactions
		if row[INTERACTION_TYPE_KEY] == SIGNIF_TYPE:

			# Add new condition if not exists
			if not dct.has_key(row[CONDITION_KEY]):
				dct[row[CONDITION_KEY]] = {}

			# Read the file and add the signif name according to the code
			with open(row[FILE_PATH_KEY], "rb") as fl:

				reader = csv.reader(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

				# Go over the rows and for each name add the maximal amount of reads given
				for file_row in reader:
					name = file_row[CODE_NAME_INDEX]

					if name not in dct[row[CONDITION_KEY]]:
						dct[row[CONDITION_KEY]][name] = file_row[LIBS_INDEX]

					else:
						dct[row[CONDITION_KEY]][name] = max(file_row[LIBS_INDEX],
															dct[row[CONDITION_KEY]][name])

	# Make sure that there are no duplicate names
	for key, val in dct.items():
		dct[key] = list(set(val))

	return dct


def generate_known_pairs_list(known_path_file):
	FIRST_NAME_INDEX = 0
	SECOND_NAME_INDEX = 1
	result = []

	with open(known_path_file, "rb") as fl:

		fl.readline()
		reader = csv.reader(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

		for row in reader:

			if row[FIRST_NAME_INDEX] > row[SECOND_NAME_INDEX]:
				result.append("%s::%s" % (row[SECOND_NAME_INDEX], row[FIRST_NAME_INDEX]))

			else:
				result.append("%s::%s" % (row[FIRST_NAME_INDEX], row[SECOND_NAME_INDEX]))

	return result


def generate_binding_pairs_list(binding_path_file):
	FIRST_NAME_INDEX = 0
	SECOND_NAME_INDEX = 3
	result = []

	with open(binding_path_file, "rb") as fl:

		fl.readline()
		reader = csv.reader(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

		for row in reader:

			if row[FIRST_NAME_INDEX] > row[SECOND_NAME_INDEX]:
				result.append("%s::%s" % (row[SECOND_NAME_INDEX], row[FIRST_NAME_INDEX]))

			else:
				result.append("%s::%s" % (row[FIRST_NAME_INDEX], row[SECOND_NAME_INDEX]))

	return result


def generate_questionable_pairs_list(questionable_path_file):
	FIRST_NAME_INDEX = 0
	SECOND_NAME_INDEX = 3
	result = []

	with open(questionable_path_file, "rb") as fl:

		fl.readline()
		reader = csv.reader(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

		for row in reader:

			if row[FIRST_NAME_INDEX] > row[SECOND_NAME_INDEX]:
				result.append("%s::%s" % (row[SECOND_NAME_INDEX], row[FIRST_NAME_INDEX]))

			else:
				result.append("%s::%s" % (row[FIRST_NAME_INDEX], row[SECOND_NAME_INDEX]))

	return result


def count_interactions_for_file(file_path):
	dct = {}

	with open(file_path, "rb") as fl:

		fl.readline()

		reader = csv.reader(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

		# Go over the rows and fill the dictionary code name to counts
		for row in reader:
			# Add the code name if we didn't have
			if row[CODE_NAME_INDEX] not in dct:
				dct[row[CODE_NAME_INDEX]] = 0

			dct[row[CODE_NAME_INDEX]] += int(row[INTERACTIONS_INDEX])

	return dct


def max_libs_for_file(file_path):
	dct = {}

	with open(file_path, "rb") as fl:

		fl.readline()

		reader = csv.reader(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

		# Go over the rows and fill the dictionary code name to counts
		for row in reader:
			# Add the code name if we didn't have
			if row[CODE_NAME_INDEX] not in dct:
				dct[row[CODE_NAME_INDEX]] = 0

			dct[row[CODE_NAME_INDEX]] = max(dct[row[CODE_NAME_INDEX]], row[LIBS_INDEX])

	return dct


def is_signif(name, condition, signif_dictionary):
	return name in signif_dictionary[condition]


def is_known(name, known_list):
	return name in known_list


def is_questionable(name, questionable_list):
	return name in questionable_list


def is_binding(name, binding_list):
	return name in binding_list


def generate_table(interaction_type,
				   condition_totals,
				   lib_counts,
				   signif_lib_counts,
				   categories_results,
				   categories_lists):

	func_table = {DIV_OPTION_SIGNIF: is_signif,
				  DIV_OPTION_KNOWN_TARGETS: is_known,
				  DIV_OPTION_BINDING: is_binding,
				  DIV_OPTION_QUESTIONABLE: is_questionable}

	for condition_name, total_dict, in condition_totals.items():

		# Build the file name
		file_name = "output/%s_%s" % (interaction_type, condition_name)

		for category, result in categories_results.items():
			file_name += "_%s_%s" % (category, result)

		file_name += ".table"

		with open(file_name, "wb") as fl:

			writer = csv.writer(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			writer.writerow(["code_name", "interactions", "max_libs"])

			#
			for name, count in total_dict.items():

				should_write = True

				# write row only if match all conditions
				for category, category_list in categories_lists.items():

					# Check special option of signif
					if category == DIV_OPTION_SIGNIF:
						if func_table[category](name, condition_name, category_list) != categories_results[category]:
							should_write = False
							break

					# Handle the rest
					else:
						if func_table[category](name, category_list) != categories_results[category]:
							should_write = False
							break

				if should_write:
					writer.writerow([name,
									 count,
									 lib_counts[condition_name][name],
									 signif_lib_counts[condition_name][name]])


class CountHandler(RowHandler):

	def __init__(self):
		RowHandler.__init__(self)
		self.dct = {}

	def handle(self, row):
		# Add the code name if we didn't have
		if row[CODE_NAME_INDEX] not in self.dct:
			self.dct[row[CODE_NAME_INDEX]] = 0

		self.dct[row[CODE_NAME_INDEX]] += int(row[INTERACTIONS_INDEX])


class MaxLibHandler(RowHandler):
	def __init__(self):
		RowHandler.__init__(self)
		self.dct = {}

	def handle(self, row):
		# Add the code name if we didn't have
		if row[CODE_NAME_INDEX] not in self.dct:
			self.dct[row[CODE_NAME_INDEX]] = 0

		self.dct[row[CODE_NAME_INDEX]] = max(self.dct[row[CODE_NAME_INDEX]], row[LIBS_INDEX])


class UnifiedLibHandler(RowHandler):
	def __init__(self, name_to_libs_dictionary):
		RowHandler.__init__(self)
		self.dct = {}
		self.name_to_libs_dictionary = name_to_libs_dictionary
		
	def handle(self, row):
		name = row[CODE_NAME_INDEX]

		if name in self.dct:
			print "-->", name

		# Add the code name if we didn't have
		if  name not in self.name_to_libs_dictionary:
			self.dct[name] = 0

		else:
			self.dct[name] = self.name_to_libs_dictionary[name]


def run():

	settings, args = process_command_line(None)

	try:
		config_rows = read_config_file_as_table(settings.config_file)

		if settings.unified_config_file:
			unified_config_rows = read_config_file_as_table(settings.unified_config_file)

	except ConfigException as ex:
		print ex
		return 1

	interaction_type_to_condition_files_dictionary = get_condition_files(config_rows)

	signif_dictionary = generate_signif_per_conditions_dictionary(config_rows)
	# max_lib_dictionary = generate_max_lib_per_condition_dictionary(config_rows)

	# Mark whether option was set
	categories_lists = {}

	# categories_lists = {DIV_OPTION_SIGNIF: signif_dictionary}
	if settings.signif:
		categories_lists[DIV_OPTION_SIGNIF] = signif_dictionary

	if settings.known_targets_file:
		known_pair_list = generate_known_pairs_list(settings.known_targets_file)
		categories_lists[DIV_OPTION_KNOWN_TARGETS] = known_pair_list

	if settings.no_binding_file:
		binding_pair_list = generate_binding_pairs_list(settings.no_binding_file)
		categories_lists[DIV_OPTION_BINDING] = binding_pair_list

	if settings.questionable_file:
		questionable_pair_list = generate_questionable_pairs_list(settings.questionable_file)
		categories_lists[DIV_OPTION_QUESTIONABLE] = questionable_pair_list

	interaction_to_condition_totals = {}
	interaction_to_condition_libs = {}
	interaction_to_condition_signif_libs = {}

	# Go over the files and make the total count of interactions per interaction type
	for interaction_type, condition_files_dictionary in interaction_type_to_condition_files_dictionary.items():

		interaction_to_condition_totals[interaction_type] = {}
		interaction_to_condition_libs[interaction_type] = {}
		interaction_to_condition_signif_libs[interaction_type] = {}

		for condition_name, file_list in condition_files_dictionary.items():

			interaction_to_condition_totals[interaction_type][condition_name] = {}
			interaction_to_condition_libs[interaction_type][condition_name] = {}
			interaction_to_condition_signif_libs[interaction_type][condition_name] = {}

			for file_path in file_list:

				fp = FileProcessor(file_path, {"read_count": CountHandler(),
											   "signif_lib_count": MaxLibHandler()})

				fp.process()

				interactions_dct = fp.row_handlers["read_count"].dct
				max_signif_lib_dct = fp.row_handlers["signif_lib_count"].dct

				# Handle total count
				for name, count in interactions_dct.items():

					# Add the name for the first time
					if name not in interaction_to_condition_totals[interaction_type][condition_name]:
						interaction_to_condition_totals[interaction_type][condition_name][name] = 0

					# Increase the count
					interaction_to_condition_totals[interaction_type][condition_name][name] += count

				# Handle signif libs
				for name, lib_count in max_signif_lib_dct.items():

					# Add the name for the first time
					if name not in interaction_to_condition_signif_libs[interaction_type][condition_name]:
						interaction_to_condition_signif_libs[interaction_type][condition_name][name] = 0

					# Set the maximal lib
					interaction_to_condition_signif_libs[interaction_type][condition_name][name] = \
						max(interaction_to_condition_signif_libs[interaction_type][condition_name][name], lib_count)

				# Handle total libs
				for name, lib_count in interactions_dct.items():

					# Add the name for the first time
					if name not in interaction_to_condition_libs[interaction_type][condition_name]:
						interaction_to_condition_libs[interaction_type][condition_name][name] = 0

					# Set the maximal lib
					interaction_to_condition_libs[interaction_type][condition_name][name] += 1


	# Only if unified config file was inserted
	if settings.unified_config_file:

		unified_interaction_type_to_condition_files_dictionary = get_condition_files(unified_config_rows)

		unified_interaction_to_condition_totals = {}
		unified_interaction_to_condition_libs = {}

		for interaction_type, condition_files_dictionary in unified_interaction_type_to_condition_files_dictionary.items():

			unified_interaction_to_condition_totals[interaction_type] = {}
			unified_interaction_to_condition_libs[interaction_type] = {}

			for condition_name, file_list in condition_files_dictionary.items():

				unified_interaction_to_condition_totals[interaction_type][condition_name] = {}

				for file_path in file_list:

					fp = FileProcessor(file_path, {"lib_count": UnifiedLibHandler(interaction_to_condition_libs[interaction_type][condition_name])})

					fp.process()

					unified_interaction_to_condition_libs[interaction_type][condition_name] = \
						fp.row_handlers["lib_count"].dct

	# Generate group tables per condition
	for interaction_type, condition_totals in interaction_to_condition_totals.items():

		perms = []

		# Generate all permutations for the flags
		for i in range(len(categories_lists) + 1):
			values = [True] * i + [False] * (len(categories_lists) - i)
			perms += list(set(itertools.permutations(values)))

		# Generate the different groups according to t
		for permutation in perms:

			categories_results = {category: permutation[i] for i, category in enumerate(categories_lists.keys())}

			generate_table(interaction_type,
						   condition_totals,
						   interaction_to_condition_libs[interaction_type],
						   interaction_to_condition_signif_libs[interaction_type],
						   categories_results,
						   categories_lists)

	# Generate unified tables
	if settings.unified_config_file:

		for interaction_type, condition_to_libs in unified_interaction_to_condition_libs.items():

			for condition_name, names_to_libs in condition_to_libs.items():

				with open("output/unified_%s_%s.table" % (interaction_type, condition_name), "wb") as fl:

					writer = csv.writer(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
					writer.writerow(["libs", "known", "unknown"])

					libs_max = max(set([int(val) for key, val in names_to_libs.items()]))

					for i in range(libs_max + 1):
						writer.writerow([i,
										 sum([1 for name, val in names_to_libs.items() if val == i and is_known(name, categories_lists[DIV_OPTION_KNOWN_TARGETS])]),
										 sum([1 for name, val in names_to_libs.items() if val == i and not is_known(name, categories_lists[DIV_OPTION_KNOWN_TARGETS])])])

	return 0


if __name__ == "__main__":
	exit(run())
