__author__ = 'amirbar'

from Table import Table, GeneTable
from TableLoader import  GeneTableLoader, LybeckerS2TableLoader, LybeckerTableLoader, TableLoader
from Globals import TableGlobals


def extract_our_gens():
    files = ["assign-type-to-all-chimeras-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_all_interactions.with-type",
             "assign-type-to-all-chimeras-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_all_interactions.with-type",
             "assign-type-to-all-chimeras-of-MG_hfq-WT101_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.with-type",
             "assign-type-to-all-chimeras-of-MG_hfq-wt202_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.with-type",
             "assign-type-to-all-chimeras-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_all_interactions.with-type",
             "assign-type-to-signif-chimeras-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_sig_interactions.with-type",
             "assign-type-to-signif-chimeras-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_sig_interactions.with-type",
             "assign-type-to-signif-chimeras-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_sig_interactions.with-type",
             "assign-type-to-single-counts-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_single_counts.with-type",
             "assign-type-to-single-counts-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_single_counts.with-type",
             "assign-type-to-single-counts-of-MG_hfq-WT101_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.with-type",
             "assign-type-to-single-counts-of-MG_hfq-wt202_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.with-type",
             "assign-type-to-single-counts-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_single_counts.with-type"]


    names_to_types = {}

    loader = TableLoader()

    for file_name in files:

        i = 0

        print file_name
        print "*" * 100

        entries = loader.load("our_files/%s" % file_name)

        for entry in entries:
            # print entry
            # i+=1
            #
            # if i == 10:
            #     break
            name = entry["rna1 name"]

            if not names_to_types.has_key(name):
                names_to_types[name] = entry["first_type"]

            name = entry["rna2 name"]

            if not names_to_types.has_key(name):
                names_to_types[name] = entry["second_type"]


    for key, value in names_to_types.items():
        print "name: %s\t|\t" % key, "type: %s" % value

def extract_raghavan(file_name, type):
    names_to_types = {}

    loader = TableLoader()

    print file_name
    print "*" * 100

    entries = loader.load("final_format/%s" % file_name)

    for entry in entries:
        # print entry
        # i+=1
        #
        # if i == 10:
        #     break
        name = entry["name"]

        if not names_to_types.has_key(name):
            names_to_types[name] = type
        else:
            print "[warning] multiple entries of the same name"


    for key, value in names_to_types.items():
        print "name: %s\t|\t" % key, "type: %s" % value

def extract_zhang(file_name, type):
    names_to_types = {}

    loader = TableLoader()

    print file_name
    print "*" * 100

    entries = loader.load("zhang/final/%s" % file_name)

    for entry in entries:
        # print entry
        # i+=1
        #
        # if i == 10:
        #     break
        name = entry["name"]

        if "/" in name:
            print name, "might not be a srna type"

        if not names_to_types.has_key(name):
            names_to_types[name] = type
        else:
            print "[warning] multiple entries of the same name"


    for key, value in names_to_types.items():
        print "name: %s\t|\t" % key, "type: %s" % value


def extract_bilusic(file_name, type):
    names_to_types = {}

    loader = TableLoader()

    print file_name
    print "*" * 100

    entries = loader.load("bilusic/final/%s" % file_name)

    for entry in entries:
        # print entry
        # i+=1
        #
        # if i == 10:
        #     break
        name = entry["name"]

        if not names_to_types.has_key(name):
            names_to_types[name] = type

            if "5'utr" in name:
                names_to_types[name] = "5'utr"
            elif "3'utr" in name:
                names_to_types[name] = "3'utr"

        else:
            print "[warning] multiple entries of the same name"



    for key, value in names_to_types.items():
        print "name: %s\t|\t" % key, "type: %s" % value

if False:
    extract_bilusic("2014RNABIOL0069R_TableS2.table", "as")
    extract_bilusic("2014RNABIOL0069R_TableS3.table", "intra")
    extract_bilusic("2014RNABIOL0069R_TableS4.table", "igr")


if False:
    extract_zhang("Table-s3-zhang-2013-sheet2008.table", "srna")
    extract_zhang("Table-s3-zhang-2013-sheet2009.table", "srna")
    extract_zhang("Table-s4-zhang-2013-sheet2008.table", "srna")
    extract_zhang("Table-s4-zhang-2013-sheet2009.table", "srna")

if False:
    extract_raghavan("s5_directed.table", "srna")
    extract_raghavan("s6_directed.table", "igr")
    extract_raghavan("s7_directed.table", "srna")