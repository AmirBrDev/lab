__author__ = 'amirbar'

from TableLoader import GeneTableLoader, PDFTableLoader, ZhangTableLoader, BilusicTableLoader, TssMasterTableLoader
from TableLoader import LybeckerTableLoader, LybeckerS2TableLoader, TableLoader
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

            first_entry_start, first_entry_end, dummy_a,\
                second_entry_start, second_entry_end, dummy_b = \
                pdfTable.findByField(table_name, entry[table_name])[0].split(";")

            startDiff = int(first_entry_start)
            endDiff = int(first_entry_end)

            first_entry_start, first_entry_end, dummy_a,\
                second_entry_start, second_entry_end, dummy_b = match[0].split(";")

            startDiff -= int (int(first_entry_start))
            endDiff -= int(first_entry_end)

            print "%s from PDF match to: %s diff is: %d;%d" % (entry[table_name], match[0], startDiff, endDiff)

            if(total_count[entry[table_name].lower()] > 1):
                print "%s has multiple instances" % entry[table_name]
        else:

            match = table.findByOtherNames(entry[table_name])

            other_name_matches += "%s match to other names: %s\n" % (entry[table_name], match[0])

            if (match[0] != None):
                first_entry_start, first_entry_end, dummy_a,\
                    second_entry_start, second_entry_end, dummy_b = \
                    pdfTable.findByField(table_name, entry[table_name])[0].split(";")

                startDiff = int(first_entry_start)
                endDiff = int(first_entry_end)

                first_entry_start, first_entry_end, dummy_a,\
                    second_entry_start, second_entry_end, dummy_b = match[0].split(";")

                startDiff -= int (int(first_entry_start))
                endDiff -= int(first_entry_end)

                other_name_matches += "the diff is %d;%d\n" % (startDiff, endDiff)

            if (total_count.has_key(entry[table_name])):

                other_name_matches += "%s appears %d in countings\n" % \
                                      (entry[table_name], total_count[entry[table_name]])

    print "-" * 100
    print "matches to old names (in case of no match to major)"
    print "-" * 100
    print other_name_matches


def lybecker_update(file_name,
                    show_warnings=True,
                    overlap_delimiter="/",
                    overlap_field="annotation of overlapping genes",
                    adjacent_field="adjacent genes",
                    loader_type=LybeckerS2TableLoader):

    geneLoader = GeneTableLoader()
    gene_table = geneLoader.createTable("genes", geneLoader.loadUnprocessed("./genes.col"))

    loader = loader_type()
    table = loader.createTable("lybecker", loader.load("lybecker/final/%s" % file_name))

    new_table_raw = []

    for id, info in table:

        dct = {}

        info.pop(Table.UNIQUE_ID_FIELD)
        start, end, strand = id.split(";")[:3]
        # print info["name"]
        start = int(start)
        end = int(end)

        result = gene_table.is_overlaps(start, end, "none")

        if show_warnings:
            print id
            overlaps = [gene for gene in info[overlap_field].split(overlap_delimiter) if gene != ""]
            print overlaps

            for index, gene in enumerate(result):
                print "%d: %s" % (index, gene[1][1]["name"])
                if gene[1][1]["name"] not in overlaps:
                    print "[warning] older name is being used"

            if len(overlaps) > len(result):
                print "[warning] missing overlapping gene"

            if len(overlaps) < len(result):
                print "[warning] extra overlapping gene"

        # Set the record location
        dct[TableGlobals.FIRST_START_BASE_KEY], dct[TableGlobals.FIRST_END_BASE_KEY], \
            dct[TableGlobals.FIRST_STRAND_KEY], dct[TableGlobals.SECOND_START_BASE_KEY], \
            dct[TableGlobals.SECOND_END_BASE_KEY], dct[TableGlobals.SECOND_STRAND_KEY] = id.split(Table.ID_DELIMITER)

        # Assume unknown strand
        dct[TableGlobals.FIRST_STRAND_KEY] = "none"
        dct[TableGlobals.SECOND_STRAND_KEY] = "none"

        is_valid = result[0][0]
        overlap_names = info[overlap_field].split(overlap_delimiter)

        # Check if major in overlapping names
        if not is_valid:
            is_valid = False

            # print "major names", [gene[1][1]["name"] for gene in result]

            for first in [gene[1][1]["name"] for gene in result]:
                for second in overlap_names:
                    if first in second:
                        is_valid = True
                        break
                if is_valid:
                    break

        # check if minor in overlapping names
        if not is_valid:
            is_valid = False

            other_names = []

            for gene in result:
                other_names.extend(gene[1][1]["other_names"])

            # print "old names", other_names

            for first in other_names:
                for second in overlap_names:
                    if first in second:
                        is_valid = True
                        break
                if is_valid:
                    break

        if "intergenic" == info["category"]:
            is_valid = False
            adjacent_genes = [val for val in info[adjacent_field].split(overlap_delimiter) if val != ""]

            first = gene_table.findByName(adjacent_genes[0])
            second = gene_table.findByName(adjacent_genes[1])

            if first == (None, None):
                first = gene_table.findByOtherNames(adjacent_genes[0])
            if second == (None, None):
                second = gene_table.findByOtherNames(adjacent_genes[1])

            if first != (None, None):
                representing = first
            elif second != (None, None):
                representing = second
            else:
                raise BaseException("No presenting gene found")

            strand = representing[0].split(Table.ID_DELIMITER)[2]

            if strand == TableGlobals.STRAND_POSITIVE:
                dct[TableGlobals.FIRST_STRAND_KEY] = TableGlobals.STRAND_NEGATIVE
                dct[TableGlobals.SECOND_STRAND_KEY] = TableGlobals.STRAND_NEGATIVE

            else:
                dct[TableGlobals.FIRST_STRAND_KEY] = TableGlobals.STRAND_POSITIVE
                dct[TableGlobals.SECOND_STRAND_KEY] = TableGlobals.STRAND_POSITIVE

        if is_valid:
            is_valid = "divergent" not in info["category"] and \
                       "convergent" not in info["category"]


        # Update the record name if gene was found
        if is_valid:
            info["name"] = "overlapping_"
        else:
            name_id = id.split(Table.ID_DELIMITER)[:2]
            name_id.append(dct[TableGlobals.FIRST_STRAND_KEY])
            info["name"] = "lybecker_%s_%s" % (info["category"], Table.ID_DELIMITER.join(name_id))

        pos_count = 0
        neg_count = 0

        if is_valid:

            # for each gene match
            for entry in result:

                # exact match add name
                info["name"] += "%s." % entry[1][1]["name"]
                strand = entry[1][0].split(Table.ID_DELIMITER)[2]

                if TableGlobals.STRAND_NEGATIVE == strand:
                    neg_count += 1

                if TableGlobals.STRAND_POSITIVE == strand:
                    pos_count += 1

            if neg_count == 0:
                dct[TableGlobals.FIRST_STRAND_KEY] = TableGlobals.STRAND_NEGATIVE
                dct[TableGlobals.SECOND_STRAND_KEY] = TableGlobals.STRAND_NEGATIVE

            elif pos_count == 0:
                dct[TableGlobals.FIRST_STRAND_KEY] = TableGlobals.STRAND_POSITIVE
                dct[TableGlobals.SECOND_STRAND_KEY] = TableGlobals.STRAND_POSITIVE
            else:
                print "no strand match: %s" % entry[1][1]["name"]

        # remove extra . from name if match found
        if is_valid:
            info["name"] = info["name"][:-1]

        dct.update(info)

        new_table_raw.append(dct)

        if show_warnings:
            print 20 * "*"

    # for row in new_table_raw:
    #     print row

    TableLoader().createTable("updated_lybecker", new_table_raw).dump("lybecker/final/updated_%s" % file_name)


if __name__ == "__main__":

    lybecker_update("lybecker_s2.table",
                show_warnings=True)
    # lybecker_update("sd01.table",
    #                 overlap_delimiter=",",
    #                 overlap_field="overlapping genes",
    #                 loader_type=LybeckerTableLoader,
    #                 adjacent_field="neighboring genes",
    #                 show_warnings=True)
