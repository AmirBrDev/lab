__author__ = 'amirbar'

import csv

class FileProcessor(object):

	def __init__(self, file_path, row_handlers_dct):

		self.file_path = file_path
		self.row_handlers = row_handlers_dct

		for handler in row_handlers_dct:
			if type(handler) is not type(RowHandler):
				assert "Invalid class type received"

	def process(self):

		with open(self.file_path, "rb") as fl:
			fl.readline()

			reader = csv.reader(fl, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

			# Go over the rows and run the handlers
			for row in reader:

				for handler in self.row_handlers.values():
					handler.handle(row)


class RowHandler(object):

	def __init__(self):
		pass

	def handle(self, row):
		raise NotImplementedError

class PrintRowHandler(RowHandler):
	def __init__(self):
		RowHandler.__init__(self)
		self.prints = []

	def handle(self, row):
		print row
		self.prints.append(str(row))

class BlaHandler(RowHandler):
	def handle(self, row):
		print "bla"


if __name__ == "__main__":

	path = "/home/users/amirbar/lab/table_builder/our_files/MEME_results.csv"

	fp = FileProcessor(path, {"printer": PrintRowHandler(),
							  "bla": BlaHandler()})
	fp.process()

	print "-" * 50
	print fp.row_handlers["printer"].prints