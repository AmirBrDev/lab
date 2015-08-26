__author__ = 'amirbar'

from TableLoader import GeneTableLoader, PDFTableLoader
from Globals import TableGlobals
from Table import Table

def addDirectionalityToPdfTables():

    loader = GeneTableLoader()
    geneTable = loader.createTable("geneDB", loader.loadUnprocessed("genes.col"))


    pdfLoader = PDFTableLoader()

    # S5 handling
    result = pdfLoader.loadUnprocessed("./parsed_files/s5.table")

    # match for each entry in the in s5 table a strand
    for entry in result:

        id, data = geneTable.findByField(Table.NAME_FIELD, entry["name"])

        strand = None

        # if there is no data look in the old names
        if (id == None):
            id, data = geneTable.findByOtherNames(entry["name"])

        if (data != None):
            strand = id.split(Table.ID_DELIMITER)[2]

        entry[TableGlobals.FIRST_STRAND_KEY] = strand
        entry[TableGlobals.SECOND_STRAND_KEY] = strand

    pdfTable = pdfLoader.createTable("dummy", result)
    pdfTable.dump("./final_format/s5_directed.table")

    #S7 handling
    result = pdfLoader.loadUnprocessed("./parsed_files/s7.table")

    for entry in result:

        strand = entry.pop("strand")

        if (strand == "f"):
            strand = Table.POS_STRAND

        elif (strand == "r"):
            strand = Table.NEG_STRAND

        entry[TableGlobals.FIRST_STRAND_KEY] = strand
        entry[TableGlobals.SECOND_STRAND_KEY] = entry[TableGlobals.FIRST_STRAND_KEY]

    pdfTable = pdfLoader.createTable("dummy", result)
    pdfTable.dump("./final_format/s7_directed.table")

    # S6 handling
    result = pdfLoader.loadUnprocessed("./parsed_files/s6.table")
    pdfLoader.createTable("dummy", result).dump("./final_format/s6_directed.table")


def printMatchAnalysis():
    """
    This function is used to print the matched directions for the know rna from
    the pdf so we can make sure the direction we matched them is valid and there
    was no problems with the name (since there are new and old ones and the article
    is from 2012).

    :return: None
    """

    loader = GeneTableLoader()
    result = loader.loadUnprocessed("genes.col")

    table = loader.createTable("geneDB", result)

    count = {}

    for entry in result:

        name = entry["NAME".lower()].lower()

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
                warnings += ("warning: %s is an old name of %s = %d\n" % (lower_name, entry["NAME".lower()], count[lower_name]))
            elif (old_count.has_key(lower_name)):
                old_count[lower_name] += 1
                old_warnings += ("warning: old multiplicity for %s = %d\n" % (lower_name, old_count[lower_name]))
            else:
                old_count[lower_name] = 1

    print "-" * 100
    print "counting instances of major names and old names"
    print "-" * 100
    print warnings
    print "-" * 100
    print "counting instances of old names between themselves"
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

    result = pdfLoader.loadUnprocessed("./parsed_files/s5.table")
    pdfTable = pdfLoader.createTable("s5", result)

    print "-" * 100
    print "print matches to majors"
    print "-" * 100


    other_name_matches = ""

    total_count = {}
    total_count.update(count)
    total_count.update(old_count)


    #table_name = "srna-name"
    table_name = "name"

    # find matching candidates to copy their directionality
    for entry in result:

        match = table.findByField("name", entry[table_name])


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
    print "matches to old names (in case of no match to major)"
    print "-" * 100
    print other_name_matches


if (__name__ == "__main__"):
    # printMatchAnalysis()
    addDirectionalityToPdfTables()