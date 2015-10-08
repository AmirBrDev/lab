class Table(object):
	"""
	This class represents a table where the first line is the header
	and the each row represents an entry where the data for the column 
	is seperated by 'seperator'
	"""

	class TableIterator(object):
		"""
		This class is used to iterate over the table rows
		"""

		def __init__(self, table):
			self._rows = table._rows
			self._iterator = iter(self._rows)

		def next(self):
			return self._iterator.next()

		def __iter__(self):
			return self._iterator

	def __init__(self, path, seperator="\t"):

		with open(path, "rb") as fl:
			header = fl.readline().strip().split(seperator)

			self._rows = [{key : val for key, val in 
				zip(header, line.strip().split(seperator))} \
				for line in fl.readlines()]

	def __iter__(self):
		return Table.TableIterator(self)


def test_table():

	for row in Table("./candidate_list.table"):
		print row
		
if __name__ == "__main__":
	test_table()