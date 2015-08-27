__author__ = 'amirbar'

import sys
import os

from TableLoader import TableLoader, OurTableLoader
from Table import Table

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

    def crossTables(self, distance_treshold, workdir="./results"):

        table_keys = ["article name",
                      "article strand",
                      "match to",
                      "overlaps",
                      "distance",
                      "first_type",
                      "second_type",
                      "RNA1 EcoCyc ID",
                      "RNA1 name",
                      "RNA2 name",
                      "interactions",
                      "odds ratio",
                      "article table"]

        for curr_table in self._our_tables_list:

            print "*" * 100
            print curr_table._name
            records = len(curr_table._dctData)
            factor = 10
            curr_record = 0
            percentage = 0

            header = "\t".join(table_keys)

            fl = open("%s/%s.table" % (workdir, curr_table._name), "wb")
            fl.write(header)
            fl.write("\n")

            for id, additional_data in curr_table:
                first_entry_start, first_entry_end, first_entry_strand,\
                    second_entry_start, second_entry_end, second_entry_strand = id.split(Table.ID_DELIMITER)

                curr_record += 1

                if curr_record % (records / factor) == 0:
                    percentage += factor
                    print "%d percent" % percentage

                for table in self._article_tables_list:
                    match_list = table.is_overlaps(int(first_entry_start),
                                                   int(first_entry_end),
                                                   first_entry_strand)

                    for res in match_list:

                        distance = res[2]

                        if ((distance <= distance_treshold) or res[0]):

                            gene = table.findById(res[1])

                            to_print = [gene[1]["name"],
                                        gene[0].split(";")[2], # the strand
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
                                        table._name]

                            line = "\t".join(to_print)
                            fl.write(line)
                            fl.write("\n")

                    # check for the second one only if its different
                    if (first_entry_start != second_entry_start or
                        first_entry_end != second_entry_end or
                        first_entry_strand != second_entry_strand):

                        match_list = table.is_overlaps(int(second_entry_start),
                                                       int(second_entry_end),
                                                       second_entry_strand)

                        for res in match_list:
                            distance = res[2]

                            if (distance <= distance_treshold) or res[0]:

                                gene = table.findById(res[1])

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
                                            table._name]

                                line = "\t".join(to_print)
                                fl.write(line)
                                fl.write("\n")
            fl.close()


def crossRaghavan():

    files = [#"assign-type-to-all-chimeras-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_all_interactions.with-type"]
             #"assign-type-to-all-chimeras-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_all_interactions.with-type"]
             #"assign-type-to-all-chimeras-of-MG_hfq-WT101_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.with-type"]
             # "assign-type-to-all-chimeras-of-MG_hfq-wt202_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.with-type"]
             # "assign-type-to-all-chimeras-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_all_interactions.with-type"]
             # "assign-type-to-signif-chimeras-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_sig_interactions.with-type"]
             # "assign-type-to-signif-chimeras-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_sig_interactions.with-type"]
             # "assign-type-to-signif-chimeras-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_sig_interactions.with-type",
             # "assign-type-to-single-counts-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_single_counts.with-type"]
             # "assign-type-to-single-counts-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_single_counts.with-type"]
             # "assign-type-to-single-counts-of-MG_hfq-WT101_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.with-type"]
             # "assign-type-to-single-counts-of-MG_hfq-wt202_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.with-type"]
             "assign-type-to-single-counts-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_single_counts.with-type"]



    our_tables = [(name, "our_files/%s" % name) for name in files]

    # our_tables = [("assign-type-to-signif-chimeras-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_sig_interactions.with-type",
    #                "our_files/assign-type-to-signif-chimeras-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_sig_interactions.with-type")]

    article_tables = [("s5", "final_format/s5_directed.table"),
                      ("s6", "final_format/s6_directed.table"),
                      ("s7", "final_format/s7_directed.table")]

    crosser = TableCrosser(our_tables, article_tables)
    crosser.crossTables(1000)

def crossZhang(fileName):

    files = [fileName]

    # files = [#"assign-type-to-all-chimeras-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_all_interactions.with-type"]
             #"assign-type-to-all-chimeras-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_all_interactions.with-type"]
             #"assign-type-to-all-chimeras-of-MG_hfq-WT101_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.with-type"]
             # "assign-type-to-all-chimeras-of-MG_hfq-wt202_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.with-type"]
             # "assign-type-to-all-chimeras-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_all_interactions.with-type"]
             # "assign-type-to-signif-chimeras-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_sig_interactions.with-type"]
             # "assign-type-to-signif-chimeras-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_sig_interactions.with-type"]
             # "assign-type-to-signif-chimeras-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_sig_interactions.with-type",
             # "assign-type-to-single-counts-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_single_counts.with-type"]
             # "assign-type-to-single-counts-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_single_counts.with-type"]
             # "assign-type-to-single-counts-of-MG_hfq-WT101_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.with-type"]
             # "assign-type-to-single-counts-of-MG_hfq-wt202_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.with-type"]
             # "assign-type-to-single-counts-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_single_counts.with-type"]

    our_tables = [(name, "our_files/%s" % name) for name in files]

    article_tables = [("zhang2013_s3_sheet_2008", "zhang/final/Table-s3-zhang-2013-sheet2008.table"),
                      ("zhang2013_s3_sheet_2009", "zhang/final/Table-s3-zhang-2013-sheet2009.table"),
                      ("zhang2013_s4_sheet_2008", "zhang/final/Table-s4-zhang-2013-sheet2008.table"),
                      ("zhang2013_s4_sheet_2009", "zhang/final/Table-s4-zhang-2013-sheet2009.table")]

    crosser = TableCrosser(our_tables, article_tables)
    crosser.crossTables(1000, "results/zhang")

def usage():
    print """TableCrosser.py <TYPE> <FILE>
    where type can be:
    1. raghavan
    2. zhang

    and file is the name of one of our *.with-type files."""

def run():
    mode = sys.argv[1]
    fileName = sys.argv[2]

    try:
        if mode == "raghavan":
            crossRaghavan(fileName)

        elif mode == "zhang":
            crossZhang(fileName)

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




