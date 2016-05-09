import csv
import os

class ConfigException(BaseException):
	def __init__(self, msg=""):
		BaseException.__init__(self, "[Configuration Exception] %s" % msg)

def __row_valid(row):
	return True

def read_config_file_as_table(config_file, validate_row_func=__row_valid):

	if not os.path.exists(config_file):
		raise ConfigException("File %s does not exists" % config_file)

	rows = []
	
	with open(config_file, "rb") as fl:
		
		reader = csv.DictReader(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

		for i, row in enumerate(reader):
			
			if not validate_row_func(row):
				raise ConfigException("Invalid row in input file. Row number %d" % (i + 1))

			rows.append(row)


	return rows


def __false(row):
	return False

if __name__ == "__main__":
	read_config_file_as_table("bla")
