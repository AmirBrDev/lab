from Bio.TogoWS import entry

__author__ = 'amirbar'

from Table import Table

class TableLoader(object):

    FIRST_START_BASE_KEY = "start_1"
    FIRST_END_BASE_KEY = "end_1"
    SECOND_START_BASE_KEY = "start_2"
    SECOND_END_BASE_KEY = "end_2"
    FIRST_STRAND_KEY = "strand_1"
    SECOND_STRAND_KEY = "strand_2"
    STRAND_POSITIVE = "+"
    STRAND_NEGATIVE = "-"

    def load(self, path):
        raise NotImplementedError

    def createTable(self, name, row_list):
        result = Table(name)

        for row in row_list:

            first = (row.pop(TableLoader.FIRST_START_BASE_KEY),
                     row.pop(TableLoader.FIRST_END_BASE_KEY),
                     row.pop(TableLoader.FIRST_STRAND_KEY))

            second = (row.pop(TableLoader.SECOND_START_BASE_KEY),
                     row.pop(TableLoader.SECOND_END_BASE_KEY),
                     row.pop(TableLoader.SECOND_STRAND_KEY))


            result.insert(first, second, **row)

        return result

class PDFTableLoader(TableLoader):

    def load(self, path):
        result = []

        fl = open(path, "rb")

        name_list = fl.readline().split("\t")

        for line in fl.readlines():

            value_list = line.split("\t")

            dct = {}

            for key, val in zip(name_list, value_list):
                dct[key.replace("\n", "")] = val.replace("\n", "")

            result.append(dct)

        fl.close()

        for entry in result:

            entry[TableLoader.FIRST_STRAND_KEY] = None
            entry[TableLoader.SECOND_STRAND_KEY] = None
            entry[TableLoader.FIRST_START_BASE_KEY] = entry.pop('right_end')
            entry[TableLoader.FIRST_END_BASE_KEY] = entry.pop('left_end')
            entry[TableLoader.SECOND_START_BASE_KEY] = entry[TableLoader.FIRST_START_BASE_KEY]
            entry[TableLoader.SECOND_END_BASE_KEY] = entry[TableLoader.FIRST_END_BASE_KEY]

        return result

class GeneTableLoader(TableLoader):
    def load(self, path):

        ignore = ["REPLICON",
                  "SWISS-PROT-ID",
                  "GENE-CLASS",
                  "PRODUCT-NAME",
                  "GENE-CLASS"]

        result = []

        fl = open(path, "rb")

        name_list = fl.readline().split("\t")

        for line in fl.readlines():

            value_list = line.split("\t")

            dct = {"other_names" : []}

            for key, val in zip(name_list, value_list):

                if (key == "SYNONYMS"):
                    if (val != ""):
                        dct["other_names"].append(val)
                elif (key in ignore):
                    continue
                else:
                    dct[key.replace("\n", "")] = val.replace("\n", "")

            result.append(dct)

        fl.close()

        for entry in result:

            if (int(entry['START-BASE']) < int(entry['END-BASE'])):
                entry[TableLoader.FIRST_STRAND_KEY] = TableLoader.STRAND_POSITIVE
                entry[TableLoader.SECOND_STRAND_KEY] = TableLoader.STRAND_POSITIVE
            else:
                entry[TableLoader.FIRST_STRAND_KEY] = TableLoader.STRAND_NEGATIVE
                entry[TableLoader.SECOND_STRAND_KEY] = TableLoader.STRAND_NEGATIVE

            entry[TableLoader.FIRST_START_BASE_KEY] = entry.pop('START-BASE')
            entry[TableLoader.FIRST_END_BASE_KEY] = entry.pop('END-BASE')
            entry[TableLoader.SECOND_START_BASE_KEY] = entry[TableLoader.FIRST_START_BASE_KEY]
            entry[TableLoader.SECOND_END_BASE_KEY] = entry[TableLoader.FIRST_END_BASE_KEY]


        return result

class OurTableLoader(TableLoader):
    def load(self, path):
        ignore = []

        result = []

        fl = open(path, "rb")

        name_list = fl.readline().split("\t")

        for line in fl.readlines():

            dct = {}

            value_list = line.split("\t")

            for key, val in zip(name_list, value_list):
                if (key in ignore):
                    continue
                else:
                    dct[key.replace("\n", "")] = val.replace("\n", "")

            result.append(dct)

        fl.close()

        for entry in result:
            entry[TableLoader.FIRST_START_BASE_KEY] = entry.pop('RNA1 from')
            entry[TableLoader.SECOND_START_BASE_KEY] = entry.pop('RNA2 from')
            entry[TableLoader.FIRST_END_BASE_KEY] = entry.pop('RNA1 to')
            entry[TableLoader.SECOND_END_BASE_KEY] = entry.pop('RNA2 to')
            entry[TableLoader.FIRST_STRAND_KEY] = entry.pop('RNA1 strand')
            entry[TableLoader.SECOND_STRAND_KEY] = entry.pop('RNA2 strand')


        return result

if (__name__ == "__main__"):



    # loader = OurTableLoader()
    # result = loader.load("./our_files/assign-type-to-all-chimeras-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_all_interactions.with-type")
    #
    # j = 0
    # for i in result:
    #    if j == 10:
    #        break
    #    print i
    #    j += 1
    #
    #
    # table = loader.createTable("chimera", result)
    #
    #
    # matches = table.is_overlaps(2348425, 2348445, "+")
    #
    # for match in matches:
    #     res, id = match
    #     print res, table.findById(id)



    loader = GeneTableLoader()
    result = loader.load("genes.col")

    table = loader.createTable("geneDB", result)

    count = {}

    for entry in result:

        name = entry["NAME"].lower()

        if not count.has_key(name):
            count[name] = 1
        else:
            count[name] += 1

        if count[name] > 1:
            print "warning"

    old_count = {}

    warnings = ""
    old_warnings = ""

    for entry in result:

        for name in entry["other_names"]:

            lower_name = name.lower()

            if count.has_key(lower_name):
                count[lower_name] += 1
                warnings += ("warning: %s is an old name of %s = %d\n" % (lower_name, entry["NAME"], count[lower_name]))
            elif (old_count.has_key(lower_name)):
                old_count[lower_name] += 1
                old_warnings += ("warning: old multiplicity for %s = %d\n" % (lower_name, old_count[lower_name]))
            else:
                old_count[lower_name] = 1

    print "-" * 100
    print warnings
    print "-" * 100
    print old_warnings
    print "-" * 100

    # matches = table.is_overlaps(510860, 510863, "+")
    #
    # for match in matches:
    #     res, id = match
    #     print res, table.findById(id)
    #
    # print table.findByName("ybaS")

    pdfLoader = PDFTableLoader()

    result = pdfLoader.load("./parsed_files/s7.table")
    pdfTable = pdfLoader.createTable("s5", result)

    print "-" * 100
    print "-" * 100
    print "-" * 100


    other_name_matches = ""

    total_count = {}
    total_count.update(count)
    total_count.update(old_count)


    # table_name = "sRNA-name"
    table_name = "name"

    # find matching candidates to copy their directionality
    for entry in result:


        match = table.findByField("NAME", entry[table_name])


        if (match != (None, None)):

            first_entry_start, first_entry_end, dummyA,\
                second_entry_start, second_entry_end, dummyB = \
                pdfTable.findByField(table_name, entry[table_name])[0].split(";")

            startDiff = int(first_entry_start)
            endDiff = int(first_entry_end)

            first_entry_start, first_entry_end, dummyA,\
                second_entry_start, second_entry_end, dummyB = match[0].split(";")

            startDiff -= int (int(first_entry_start))
            endDiff -= int(first_entry_end)

            print "%s from PDF match to: %s diff is: %d;%d" % (entry[table_name], match[0], startDiff, endDiff)

            if(total_count[entry[table_name].lower()] > 1):
                print "%s has multiple instances" % entry[table_name]
        else:

            match = table.findByOtherNames(entry[table_name])

            other_name_matches += "%s match to other names: %s\n" % (entry[table_name], match[0])

            if (match[0] != None):
                first_entry_start, first_entry_end, dummyA,\
                    second_entry_start, second_entry_end, dummyB = \
                    pdfTable.findByField(table_name, entry[table_name])[0].split(";")

                startDiff = int(first_entry_start)
                endDiff = int(first_entry_end)

                first_entry_start, first_entry_end, dummyA,\
                    second_entry_start, second_entry_end, dummyB = match[0].split(";")

                startDiff -= int (int(first_entry_start))
                endDiff -= int(first_entry_end)

                other_name_matches += "the diff is %d;%d\n" % (startDiff, endDiff)

            if (total_count.has_key(entry[table_name])):

                other_name_matches += "%s appears %d in countings\n" % (entry[table_name], total_count[entry[table_name]])

    print "-" * 100
    print other_name_matches

    # result = PDFTableLoader().load("./parsed_files/s5.table")
    #
    # for entry in result:
    #     print entry