__author__ = 'amirbar'

import sys
import time

from TableLoader import TableLoader, OurTableLoader
from Table import Table
from TableFormatter import TableFormatter


class TableCrosser(object):
    """
    This class responsible for crossing the different tables we work with
    and print the result nicely
    """

    def __init__(self, our_tables_list, article_tables_list):

        loader = OurTableLoader()

        self._our_tables_list = \
            [loader.createTable(name, loader.loadUnprocessed(path)) for name, path in our_tables_list]

        loader = TableLoader()
        self._article_tables_list = [loader.createTable(name, loader.load(path)) for name, path in article_tables_list]

    def cross_tables(self, distance_treshold, workdir="./results", factor=10):

        table_keys = ["article name",
                      "article strand",
                      "match to",
                      "overlaps",
                      "distance",
                      "first_type",
                      "second_type",
                      "RNA1 EcoCyc ID",
                      "RNA2 EcoCyc ID",
                      "RNA1 name",
                      "RNA2 name",
                      "interactions",
                      "odds ratio",
                      "article table"]

        for curr_table in self._our_tables_list:

            print "*" * 100
            print curr_table.get_name()
            records = len(curr_table)
            curr_record = 0
            percentage = 0

            start = time.time()
            header = "\t".join(table_keys)

            fl = open("%s/%s.table" % (workdir, curr_table.get_name()), "wb")
            fl.write(header)
            fl.write("\n")

            for id, additional_data in curr_table:
                first_entry_start, first_entry_end, first_entry_strand,\
                    second_entry_start, second_entry_end, second_entry_strand = id.split(Table.ID_DELIMITER)

                curr_record += 1

                if curr_record % (records / factor) == 0:
                    percentage += 100 / factor
                    end = time.time()
                    expected = (float(end - start) / (float(percentage) / 100)) * (float(100 - percentage) / 100)
                    print "%d percent, time elapsed %ds, remaining: %ds" % (percentage, end - start, expected)

                for table in self._article_tables_list:
                    match_list = table.is_overlaps(int(first_entry_start),
                                                   int(first_entry_end),
                                                   first_entry_strand)

                    for res in match_list:

                        distance = res[2]

                        if (distance <= distance_treshold) or res[0]:

                            # gene = table.findById(res[1])
                            gene = res[1]

                            to_print = [gene[1]["name"],
                                        gene[0].split(";")[2],  # the strand
                                        "RNA1",
                                        str(res[0]),
                                        str(distance),
                                        additional_data["first_type"],
                                        additional_data["second_type"],
                                        additional_data["rna1 ecocyc id"],
                                        additional_data["rna2 ecocyc id"],
                                        additional_data["rna1 name"],
                                        additional_data["rna2 name"],
                                        additional_data["interactions"],
                                        additional_data["odds ratio"],
                                        table.get_name()]

                            line = "\t".join(to_print)
                            fl.write(line)
                            fl.write("\n")

                    # check for the second one only if its different
                    if (first_entry_start != second_entry_start or
                            first_entry_end != second_entry_end or first_entry_strand != second_entry_strand):

                        match_list = table.is_overlaps(int(second_entry_start),
                                                       int(second_entry_end),
                                                       second_entry_strand)

                        for res in match_list:
                            distance = res[2]

                            if (distance <= distance_treshold) or res[0]:

                                # gene = table.findById(res[1])
                                gene = res[1]

                                to_print = [gene[1]["name"],
                                            gene[0].split(";")[2],  # the strand
                                            "RNA2",
                                            str(res[0]),
                                            str(distance),
                                            additional_data["first_type"],
                                            additional_data["second_type"],
                                            additional_data["rna1 ecocyc id"],
                                            additional_data["rna2 ecocyc id"],
                                            additional_data["rna1 name"],
                                            additional_data["rna2 name"],
                                            additional_data["interactions"],
                                            additional_data["odds ratio"],
                                            table.get_name()]

                                line = "\t".join(to_print)
                                fl.write(line)
                                fl.write("\n")
            fl.close()

    def cross_tables_reverse(self, distance_treshold, workdir="./results", factor=10):

        table_keys = ["RNA1 name",
                      "RNA2 name",
                      "RNA1 EcoCyc ID",
                      "RNA2 EcoCyc ID",
                      "first_type",
                      "second_type",
                      "overlaps",
                      "article name",
                      "article strand",
                      "distance",
                      "match to",
                      "interactions",
                      "odds ratio",
                      "article table"]

        for curr_table in self._our_tables_list:

            print "*" * 100
            print curr_table.get_name()
            records = len(curr_table)
            curr_record = 0
            percentage = 0

            start = time.time()
            header = "\t".join(table_keys)

            fl = open("%s/%s.table" % (workdir, curr_table.get_name()), "wb")
            fl.write(header)
            fl.write("\n")

            for id, additional_data in curr_table:
                first_entry_start, first_entry_end, first_entry_strand,\
                    second_entry_start, second_entry_end, second_entry_strand = id.split(Table.ID_DELIMITER)

                curr_record += 1

                if curr_record % (records / factor) == 0:
                    percentage += 100 / factor
                    end = time.time()
                    expected = (float(end - start) / (float(percentage) / 100)) * (float(100 - percentage) / 100)
                    print "%d percent, time elapsed %ds, remaining: %ds" % (percentage, end - start, expected)

                for table in self._article_tables_list:
                    match_list = table.is_overlaps(int(first_entry_start),
                                                   int(first_entry_end),
                                                   first_entry_strand)

                    for res in match_list:

                        distance = res[2]

                        if (distance <= distance_treshold) or res[0]:

                            # gene = table.findById(res[1])
                            gene = res[1]

                            to_print = [additional_data["rna1 name"],
                                        additional_data["rna2 name"],
                                        additional_data["rna1 ecocyc id"],
                                        additional_data["rna2 ecocyc id"],
                                        additional_data["first_type"],
                                        additional_data["second_type"],
                                        str(res[0]),
                                        gene[1]["name"],
                                        gene[0].split(";")[2],  # the strand
                                        str(distance),
                                        "RNA1",
                                        additional_data["interactions"],
                                        additional_data["odds ratio"],
                                        table.get_name()]

                            line = "\t".join(to_print)
                            fl.write(line)
                            fl.write("\n")

                    # check for the second one only if its different
                    if (first_entry_start != second_entry_start or
                            first_entry_end != second_entry_end or first_entry_strand != second_entry_strand):

                        match_list = table.is_overlaps(int(second_entry_start),
                                                       int(second_entry_end),
                                                       second_entry_strand)

                        for res in match_list:
                            distance = res[2]

                            if (distance <= distance_treshold) or res[0]:

                                # gene = table.findById(res[1])
                                gene = res[1]

                                to_print = [additional_data["rna1 name"],
                                            additional_data["rna2 name"],
                                            additional_data["rna1 ecocyc id"],
                                            additional_data["rna2 ecocyc id"],
                                            additional_data["first_type"],
                                            additional_data["second_type"],
                                            str(res[0]),
                                            gene[1]["name"],
                                            gene[0].split(";")[2],  # the strand
                                            str(distance),
                                            "RNA2",
                                            additional_data["interactions"],
                                            additional_data["odds ratio"],
                                            table.get_name()]

                                line = "\t".join(to_print)
                                fl.write(line)
                                fl.write("\n")
            fl.close()


def cross_raghavan(file_name, reverse=False):

    print "Cross Raghavan running..."

    files = [file_name]

    our_tables = [(name, "our_files/%s" % name) for name in files]

    article_tables = [("s5", "final_format/s5_directed.table"),
                      ("s6", "final_format/s6_directed.table"),
                      ("s7", "final_format/s7_directed.table")]

    crosser = TableCrosser(our_tables, article_tables)

    if not reverse:
        crosser.cross_tables(1000, "results/raghavan")
    else:
        crosser.cross_tables_reverse(100000000, "results/raghavan/reverse")
        TableFormatter.split_reverse_rows(TableLoader, "results/raghavan/reverse/%s.table")


def cross_zhang(file_name, reverse=False):

    print "Cross Zhang running..."

    files = [file_name]

    our_tables = [(name, "our_files/%s" % name) for name in files]

    article_tables = [("zhang2013_s3_sheet_2008", "zhang/final/Table-s3-zhang-2013-sheet2008.table"),
                      ("zhang2013_s3_sheet_2009", "zhang/final/Table-s3-zhang-2013-sheet2009.table"),
                      ("zhang2013_s4_sheet_2008", "zhang/final/Table-s4-zhang-2013-sheet2008.table"),
                      ("zhang2013_s4_sheet_2009", "zhang/final/Table-s4-zhang-2013-sheet2009.table")]

    crosser = TableCrosser(our_tables, article_tables)

    if not reverse:
        crosser.cross_tables(1000, "results/zhang")
    else:
        crosser.cross_tables_reverse(100000000, "results/zhang/reverse")
        TableFormatter.split_reverse_rows(TableLoader, "results/zhang/reverse/%s.table")


def cross_bilusic(file_name, reverse=False):

    print "Cross Bilusic running..."

    files = [file_name]

    our_tables = [(name, "our_files/%s" % name) for name in files]

    article_tables = [("2014RNABIOL0069R_TableS1", "bilusic/final/2014RNABIOL0069R_TableS1.table"),
                      ("2014RNABIOL0069R_TableS2", "bilusic/final/2014RNABIOL0069R_TableS2.table"),
                      ("2014RNABIOL0069R_TableS3", "bilusic/final/2014RNABIOL0069R_TableS3.table"),
                      ("2014RNABIOL0069R_TableS4", "bilusic/final/2014RNABIOL0069R_TableS4.table")]

    crosser = TableCrosser(our_tables, article_tables)

    if not reverse:
        crosser.cross_tables(1000, "results/bilusic", 20)
    else:
        crosser.cross_tables_reverse(100000000, "results/bilusic/reverse", 20)
        TableFormatter.split_reverse_rows(TableLoader, "results/bilusic/reverse/%s.table")


def cross_tss(file_name, reverse=False):

    print "Cross Tss running..."

    files = [file_name]

    our_tables = [(name, "our_files/%s" % name) for name in files]

    article_tables = [("Tss master table", "tss/final/JB.02096-14_zjb999093409sd1-3.table")]

    crosser = TableCrosser(our_tables, article_tables)

    if not reverse:
        crosser.cross_tables(1000, "results/tss", 20)
    else:
        crosser.cross_tables_reverse(100000000, "results/tss/reverse", 20)
        TableFormatter.split_reverse_rows(TableLoader, "results/tss/reverse/%s.table")


def cross_lybecker(file_name, reverse=False):

    print "Cross Lybecker running..."

    files = [file_name]

    our_tables = [(name, "our_files/%s" % name) for name in files]

    article_tables = [("Lybecker sd01", "lybecker/final/updated_sd01.table"),
                      ("Lybecker s2", "lybecker/final/updated_lybecker_s2.table")]

    crosser = TableCrosser(our_tables, article_tables)

    if not reverse:
        crosser.cross_tables(1000, "results/lybecker", 20)

    else:
        crosser.cross_tables_reverse(100000000, "results/lybecker/reverse", 20)
        TableFormatter.split_reverse_rows(TableLoader, "results/lybecker/reverse/%s.table")


def usage():

    print """TableCrosser.py <TYPE> <FILE>
    where type can be:
    1. raghavan
    2. zhang

    and file is the name of one of our *.with-type files."""


def run():

    mode = sys.argv[1]
    file_name = sys.argv[2]

    try:
        if mode == "raghavan":
            cross_raghavan(file_name)

        elif mode == "zhang":
            cross_zhang(file_name)

        else:
            print "invalid mode."

    except IOError as ex:
        print ex.message


if __name__ == "__main__":

    expected_args = 3

    print sys.argv

    if len(sys.argv) != expected_args:
        usage()

    else:
        run()

    print "Done"
