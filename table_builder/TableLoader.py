from Bio.TogoWS import entry

__author__ = 'amirbar'

from Table import Table, GeneTable
from Globals import TableGlobals
from BuilderExceptions.StructureException import StructureException

class TableLoader(object):
    """
    Loads tables from a file. The loaded tables has the following format:
    first row is the column names delimited by tabs.
    following rows are the values for the column delimited by tabs.

    Example:
    COL1    COL2
    VAL1    VAL2
    VAL3    VAL4

    where COL1, COL2 are the column names and the VAL#X is the values
    """

    def __init__(self):
        """
        Sets the default table type
        :return:
        """
        self._tableType = Table

    def load(self, path):
        """
        Loads a table of the format described in the class description
        :param path: the path to the file
        :return: a list containing the data as dictionary
        """
        result = []

        fl = open(path, "rb")

        name_list = fl.readline().lower().replace("\n", "").split(Table.TABLE_DELIMITER)

        for line in fl.readlines():

            line = line.lower()

            dct = {}

            value_list = line.split(Table.TABLE_DELIMITER)

            if (len(name_list) != len(value_list)):
                raise StructureException("Unmatched arguments count in table.")

            for key, val in zip(name_list, value_list):
                dct[key.replace("\n", "")] = val.replace("\n", "")

            result.append(dct)

        fl.close()

        return  result

    def loadUnprocessed(self, path):
        """
        Loads unprocessed tables
        This function should be implemented for each type of inheriting loader
        :param path: the file path
        :return:
        """
        raise NotImplementedError

    def createTable(self, name, row_list, append_dup = False):
        """
        Creates a table from output of the load
        :param name: the name for the table
        :param row_list: the output of the load method
        :return: a table object
        """
        result = self._tableType(name)

        for row in row_list:

            first = (row.pop(TableGlobals.FIRST_START_BASE_KEY),
                     row.pop(TableGlobals.FIRST_END_BASE_KEY),
                     row.pop(TableGlobals.FIRST_STRAND_KEY))

            second = (row.pop(TableGlobals.SECOND_START_BASE_KEY),
                     row.pop(TableGlobals.SECOND_END_BASE_KEY),
                     row.pop(TableGlobals.SECOND_STRAND_KEY))


            result.insert(first, second, append_dup, **row)

        return result

class TssMasterTableLoader(TableLoader):
    def loadUnprocessed(self, path):

        result = self.load(path)

        for entry in result:
            entry["name"] = entry.pop('locus_tag')
            entry[TableGlobals.FIRST_START_BASE_KEY] = entry.pop('pos')
            entry[TableGlobals.FIRST_END_BASE_KEY] = entry[TableGlobals.FIRST_START_BASE_KEY]
            entry[TableGlobals.FIRST_STRAND_KEY] = entry.pop('strand')
            entry[TableGlobals.SECOND_START_BASE_KEY] = entry[TableGlobals.FIRST_START_BASE_KEY]
            entry[TableGlobals.SECOND_END_BASE_KEY] = entry[TableGlobals.FIRST_END_BASE_KEY]
            entry[TableGlobals.SECOND_STRAND_KEY] = entry[TableGlobals.FIRST_STRAND_KEY]

        return result

class LybeckerTableLoader(TableLoader):
    def loadUnprocessed(self, path):

        result = self.load(path)

        for entry in result:
            entry["name"] = "None"

            start, end = entry.pop("coordinates").split("-")

            print start, end, (int(start) > int(end))

            if int(start) > int(end):
                entry[TableGlobals.FIRST_STRAND_KEY] = Table.NEG_STRAND
                entry[TableGlobals.FIRST_START_BASE_KEY] = end
                entry[TableGlobals.FIRST_END_BASE_KEY] = start
            else:
                entry[TableGlobals.FIRST_STRAND_KEY] = Table.POS_STRAND
                entry[TableGlobals.FIRST_START_BASE_KEY] = start
                entry[TableGlobals.FIRST_END_BASE_KEY] = end

            entry[TableGlobals.SECOND_START_BASE_KEY] = entry[TableGlobals.FIRST_START_BASE_KEY]
            entry[TableGlobals.SECOND_END_BASE_KEY] = entry[TableGlobals.FIRST_END_BASE_KEY]
            entry[TableGlobals.SECOND_STRAND_KEY] = entry[TableGlobals.FIRST_STRAND_KEY]

        return result

class LybeckerS2TableLoader(TableLoader):
    def loadUnprocessed(self, path):

        result = self.load(path)

        for entry in result:
            entry["name"] = "None"

            start, end = entry.pop("coordinates ip-dsrnas").replace(",", "").split("-")

            if int(start) > int(end):
                entry[TableGlobals.FIRST_STRAND_KEY] = Table.NEG_STRAND
                entry[TableGlobals.FIRST_START_BASE_KEY] = end
                entry[TableGlobals.FIRST_END_BASE_KEY] = start
            else:
                entry[TableGlobals.FIRST_STRAND_KEY] = Table.POS_STRAND
                entry[TableGlobals.FIRST_START_BASE_KEY] = start
                entry[TableGlobals.FIRST_END_BASE_KEY] = end

            entry[TableGlobals.SECOND_START_BASE_KEY] = entry[TableGlobals.FIRST_START_BASE_KEY]
            entry[TableGlobals.SECOND_END_BASE_KEY] = entry[TableGlobals.FIRST_END_BASE_KEY]
            entry[TableGlobals.SECOND_STRAND_KEY] = entry[TableGlobals.FIRST_STRAND_KEY]

        return result

class BilusicTableLoader(TableLoader):

    def loadUnprocessed(self, path):

        result = []

        fl = open(path, "rb")

        name_list = fl.readline().lower().split(Table.TABLE_DELIMITER)

        for line in fl.readlines():

            line = line.lower()

            value_list = line.split(Table.TABLE_DELIMITER)

            dct = {}

            for key, val in zip(name_list, value_list):
                dct[key.replace("\n", "")] = val.replace("\n", "")


            # append only non empty rows
            if dct["rna"] != "":
                result.append(dct)

        fl.close()

        for entry in result:
            entry["name"] = entry.pop("rna")

            entry[TableGlobals.FIRST_START_BASE_KEY] = entry.pop('start')
            entry[TableGlobals.FIRST_END_BASE_KEY] = entry.pop('end')
            entry[TableGlobals.FIRST_STRAND_KEY] = entry.pop('strand')
            entry[TableGlobals.SECOND_START_BASE_KEY] = entry[TableGlobals.FIRST_START_BASE_KEY]
            entry[TableGlobals.SECOND_END_BASE_KEY] = entry[TableGlobals.FIRST_END_BASE_KEY]
            entry[TableGlobals.SECOND_STRAND_KEY] = entry[TableGlobals.FIRST_STRAND_KEY]

        return result

class ZhangTableLoader(TableLoader):

    def loadUnprocessed(self, path):

        result = self.load(path)

        for entry in result:

            entry["name"] = entry.pop("gene")

            entry[TableGlobals.FIRST_START_BASE_KEY] = entry.pop('gene start')
            entry[TableGlobals.FIRST_END_BASE_KEY] = entry.pop('gene stop')
            entry[TableGlobals.FIRST_STRAND_KEY] = entry.pop('strand')
            entry[TableGlobals.SECOND_START_BASE_KEY] = entry[TableGlobals.FIRST_START_BASE_KEY]
            entry[TableGlobals.SECOND_END_BASE_KEY] = entry[TableGlobals.FIRST_END_BASE_KEY]
            entry[TableGlobals.SECOND_STRAND_KEY] = entry[TableGlobals.FIRST_STRAND_KEY]


        return result

class PDFTableLoader(TableLoader):

    S5_NAME_FIELD = "srna-name"

    def loadUnprocessed(self, path):
        result = self.load()

        for entry in result:

            if (entry.has_key(PDFTableLoader.S5_NAME_FIELD)):
                entry["name"] = entry.pop(PDFTableLoader.S5_NAME_FIELD)

            entry[TableGlobals.FIRST_STRAND_KEY] = None
            entry[TableGlobals.SECOND_STRAND_KEY] = None
            entry[TableGlobals.FIRST_START_BASE_KEY] = entry.pop('right_end')
            entry[TableGlobals.FIRST_END_BASE_KEY] = entry.pop('left_end')
            entry[TableGlobals.SECOND_START_BASE_KEY] = entry[TableGlobals.FIRST_START_BASE_KEY]
            entry[TableGlobals.SECOND_END_BASE_KEY] = entry[TableGlobals.FIRST_END_BASE_KEY]

            if (not entry.has_key("name")):
                entry["name"] = "igr_%s-%s" % \
                                (entry[TableGlobals.FIRST_START_BASE_KEY], entry[TableGlobals.FIRST_END_BASE_KEY])

        return result

class GeneTableLoader(TableLoader):

    def __init__(self):
        self._tableType = GeneTable

    def loadUnprocessed(self, path):

        ignore = ["REPLICON".lower(),
                  "SWISS-PROT-ID".lower(),
                  "GENE-CLASS".lower(),
                  "PRODUCT-NAME".lower(),
                  "GENE-CLASS".lower()]

        result = []

        fl = open(path, "rb")

        name_list = fl.readline().lower().split(Table.TABLE_DELIMITER)

        for line in fl.readlines():

            line = line.lower()

            value_list = line.split(Table.TABLE_DELIMITER)

            dct = {"other_names" : []}

            for key, val in zip(name_list, value_list):

                if (key == "SYNONYMS".lower()):
                    if (val != ""):
                        dct["other_names"].append(val)
                elif (key in ignore):
                    continue
                else:
                    dct[key.replace("\n", "")] = val.replace("\n", "")

            result.append(dct)

        fl.close()

        for entry in result:

            if (int(entry['START-BASE'.lower()]) < int(entry['END-BASE'.lower()])):
                entry[TableGlobals.FIRST_STRAND_KEY] = TableGlobals.STRAND_POSITIVE
                entry[TableGlobals.SECOND_STRAND_KEY] = TableGlobals.STRAND_POSITIVE
            else:
                entry[TableGlobals.FIRST_STRAND_KEY] = TableGlobals.STRAND_NEGATIVE
                entry[TableGlobals.SECOND_STRAND_KEY] = TableGlobals.STRAND_NEGATIVE

            entry[TableGlobals.FIRST_START_BASE_KEY] = entry.pop('START-BASE'.lower())
            entry[TableGlobals.FIRST_END_BASE_KEY] = entry.pop('END-BASE'.lower())
            entry[TableGlobals.SECOND_START_BASE_KEY] = entry[TableGlobals.FIRST_START_BASE_KEY]
            entry[TableGlobals.SECOND_END_BASE_KEY] = entry[TableGlobals.FIRST_END_BASE_KEY]


        return result

class OurTableLoader(TableLoader):
    def loadUnprocessed(self, path):
        ignore = []

        result = []

        fl = open(path, "rb")

        name_list = fl.readline().lower().split(Table.TABLE_DELIMITER)

        for line in fl.readlines():

            line = line.lower()

            dct = {}

            value_list = line.split(Table.TABLE_DELIMITER)

            for key, val in zip(name_list, value_list):
                if (key in ignore):
                    continue
                else:
                    dct[key.replace("\n", "")] = val.replace("\n", "")

            result.append(dct)

        fl.close()

        for entry in result:
            entry[TableGlobals.FIRST_START_BASE_KEY] = entry.pop('RNA1 from'.lower())
            entry[TableGlobals.SECOND_START_BASE_KEY] = entry.pop('RNA2 from'.lower())
            entry[TableGlobals.FIRST_END_BASE_KEY] = entry.pop('RNA1 to'.lower())
            entry[TableGlobals.SECOND_END_BASE_KEY] = entry.pop('RNA2 to'.lower())
            entry[TableGlobals.FIRST_STRAND_KEY] = entry.pop('RNA1 strand'.lower())
            entry[TableGlobals.SECOND_STRAND_KEY] = entry.pop('RNA2 strand'.lower())


        return result

if (__name__ == "__main__"):

    loader = LybeckerTableLoader()

    result = loader.loadUnprocessed("./lybecker/csv/sd01.csv")

    for entry in result:

        print entry


    # loader = OurTableLoader()
    # result = loader.load("./our_files/assign-type-to-all-chimeras-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_all_interactions.with-type")
    #
    # j = 0
    # for i in result:
    #    if j == 10:
    #        break
    #    print i
    #    j += 1


    # table = loader.createTable("chimera", result)
    #
    #
    # matches = table.is_overlaps(2348425, 2348445, "+")
    #
    # for match in matches:
    #     res, id = match
    #     print res, table.findById(id)





    # result = PDFTableLoader().loadUnprocessed("./parsed_files/s5.table")
    #
    # for entry in result:
    #     print entry