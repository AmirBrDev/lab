__author__ = 'amirbar'

from Table import Table, GeneTable
from TableLoader import  GeneTableLoader, LybeckerS2TableLoader, LybeckerTableLoader, TableLoader, OurTableLoader
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


import MySQLdb
from Table import Table

def generate_table_old(table_path, name, is_our_table=False):

    if not is_our_table:
        loader = TableLoader()
        table = loader.createTable(name, loader.load(table_path))

    else:
        loader = OurTableLoader()
        table = loader.createTable(name, loader.loadUnprocessed(table_path))

    db=MySQLdb.connect(host="localhost",user="amirbar",db="amir")
    cur = db.cursor(MySQLdb.cursors.DictCursor)

    # Generate the keys for the table
    id_keys = [TableGlobals.FIRST_START_BASE_KEY,
                   TableGlobals.FIRST_END_BASE_KEY,
                   TableGlobals.FIRST_STRAND_KEY,
                   TableGlobals.SECOND_START_BASE_KEY,
                   TableGlobals.SECOND_END_BASE_KEY,
                   TableGlobals.SECOND_STRAND_KEY]

    for key in table._dctData.values()[0].keys():
        if key not in id_keys:
            id_keys.append(key)

    fields = ", ".join("%s VARCHAR(200)" % key.replace(" ", "_").replace("-", "_").replace("'", "").replace("/", "") for key in id_keys if key != "")
    print fields
    cur.execute("CREATE TABLE %s (%s)" % (table.get_name(), fields))


    table_as_list = []

    # Generate dictionary for each row according to the keys
    for key, value in table:
        id_values = key.split(Table.ID_DELIMITER)

        for id_key, id_val in zip(id_keys, id_values):

            value[id_key] = str(id_val)

        table_as_list.append(value)
        # print "-->", value

    # Go over the rows and add them to the db
    for row in table_as_list:

        values = ",".join("%s" % db.literal(str(row[key])) for key in id_keys if key != "")
        # print "INSERT INTO %s VALUES (%s)" % (table.get_name(), values)
        cur.execute("INSERT INTO %s VALUES (%s)" % (table.get_name(), values))

    db.commit()

def generate_table(table_path, name, is_our_table=False):

    if not is_our_table:
        loader = TableLoader()
        table = loader.createTable(name, loader.load(table_path))

    else:
        loader = OurTableLoader()
        table = loader.createTable(name, loader.loadUnprocessed(table_path))

    db=MySQLdb.connect(host="localhost",user="amirbar",db="amir")
    cur = db.cursor(MySQLdb.cursors.DictCursor)

    # Generate the keys for the table
    id_keys = [TableGlobals.FIRST_START_BASE_KEY,
                   TableGlobals.FIRST_END_BASE_KEY,
                   TableGlobals.FIRST_STRAND_KEY,
                   TableGlobals.SECOND_START_BASE_KEY,
                   TableGlobals.SECOND_END_BASE_KEY,
                   TableGlobals.SECOND_STRAND_KEY]

    for key in table._dctData.values()[0].keys():
        if key not in id_keys:
            id_keys.append(key)

    key_list = []

    for key in id_keys:
        if (key == "primary"):
            key_list.append("_primary")
        elif (key == "condition"):
            key_list.append("_condition")
        else:
            key_list.append(key)

    fields = ", ".join("%s VARCHAR(200)" % key.replace(" ", "_").replace("-", "_").replace("'", "").replace("/", "") for key in key_list if (key != "" and key != 'sequence -50 nt upstream + tss (51nt)') )
    print fields
    cur.execute("CREATE TABLE %s (%s)" % (table.get_name(), fields))


    table_as_list = []

    # Generate dictionary for each row according to the keys
    for key, value in table:
        id_values = key.split(Table.ID_DELIMITER)

        for id_key, id_val in zip(id_keys, id_values):

            value[id_key] = str(id_val)

        table_as_list.append(value)
        # print "-->", value

    # Go over the rows and add them to the db
    for row in table_as_list:

        values = ",".join("%s" % db.literal(str(row[key])) for key in id_keys if (key != "" and key != 'sequence -50 nt upstream + tss (51nt)'))
        # print "INSERT INTO %s VALUES (%s)" % (table.get_name(), values)
        cur.execute("INSERT INTO %s VALUES (%s)" % (table.get_name(), values))

    db.commit()

def generate_gene_table(table_path, name):
    loader = GeneTableLoader()
    table = loader.createTable(name, loader.loadUnprocessed(table_path))

    rows = []
    for key, value in table:

        start, end, strand = key.split(Table.ID_DELIMITER)[:3]

        row = [value["unique_id"],
               value["name"],
               value["blattner-id"],
               value["unique-id"],
               int(start),
               int(end),
               strand]

        synonyms = [""] * 4

        for index, name in enumerate(value["other_names"]):
            synonyms[index] = name

        row.extend(synonyms)

        rows.append(row)

    # for row in rows:
    #     print row

    db=MySQLdb.connect(host="localhost",user="amirbar",db="amir")
    cur = db.cursor(MySQLdb.cursors.DictCursor)

    key_list = ["unique_id",
                "name",
                "blattner_id",
                "gene_unique_id",
                "start",
                "end",
                "strand",
                "synonym_1",
                "synonym_2",
                "synonym_3",
                "synonym_4"]

    fields = ", ".join("%s VARCHAR(200)" % key for key in key_list)
    cur.execute("CREATE TABLE %s (%s)" % (table.get_name(), fields))

    for row in rows:
        values = ",".join("%s" % db.literal(str(val)) for val in row )
        cur.execute("INSERT INTO %s VALUES (%s)" % (table.get_name(), values))

    db.commit()

# generate_gene_table("genes.col", "genes")

# generate_table("final_format/s5_directed.table", "raghavan_s5")
# generate_table("final_format/s6_directed.table", "raghavan_s6")
# generate_table("final_format/s7_directed.table", "raghavan_s7")

# generate_table("bilusic/final/2014RNABIOL0069R_TableS1.table", "bilusic_s1")
# generate_table("bilusic/final/2014RNABIOL0069R_TableS2.table", "bilusic_s2")
# generate_table("bilusic/final/2014RNABIOL0069R_TableS3.table", "bilusic_s3")
# generate_table("bilusic/final/2014RNABIOL0069R_TableS4.table", "bilusic_s4")

# generate_table("zhang/final/Table-s3-zhang-2013-sheet2008.table", "zhang_s3_2013_sheet2008")
# generate_table("zhang/final/Table-s3-zhang-2013-sheet2009.table", "zhang_s3_2013_sheet2009")
# generate_table("zhang/final/Table-s4-zhang-2013-sheet2008.table", "zhang_s4_2013_sheet2008")
# generate_table("zhang/final/Table-s4-zhang-2013-sheet2009.table", "zhang_s4_2013_sheet2009")

# generate_table("lybecker/final/updated_sd01.table", "lybecker_s1")
# generate_table("lybecker/final/lybecker_s2.table", "lybecker_s2")

# generate_table("tss/final/JB.02096-14_zjb999093409sd1-3.table", "tss")

def dump_our_to_sql():

    files = [("assign-type-to-all-chimeras-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_all_interactions.with-type", "all_chimeras_of_iron_limitation_cl"),
             ("assign-type-to-all-chimeras-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_all_interactions.with-type", "all_chimeras_of_log_phase_cl"),
             ("assign-type-to-all-chimeras-of-MG_hfq-WT101_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.with-type", "all_chimeras_of_mg_hfq_wt101"),
             ("assign-type-to-all-chimeras-of-MG_hfq-wt202_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.with-type", "all_chimeras_of_mg_hfq_wt202_cl_stationary"),
             ("assign-type-to-all-chimeras-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_all_interactions.with-type", "all_chimeras_of_stationary_cl"),
             ("assign-type-to-signif-chimeras-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_sig_interactions.with-type", "signif_chimeras_of_iron_limitation_cl"),
             ("assign-type-to-signif-chimeras-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_sig_interactions.with-type","signif_chimeras_of_log_phase_cl"),
             ("assign-type-to-signif-chimeras-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_sig_interactions.with-type", "signif_chimeras_of_stationary_cl"),
             ("assign-type-to-single-counts-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_single_counts.with-type", "single_of_iron_limitation_cl"),
             ("assign-type-to-single-counts-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_single_counts.with-type", "single_of_log_phase_cl"),
             ("assign-type-to-single-counts-of-MG_hfq-WT101_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.with-type", "single_of_mg_hfq_wt101"),
             ("assign-type-to-single-counts-of-MG_hfq-wt202_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.with-type", "single_of_mg_hfq_wt202_cl_stationary"),
             ("assign-type-to-single-counts-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_single_counts.with-type", "single_of_stationary_cl")]

    for file_name, table_name in files:
        print table_name
        print file_name
        print "*" * 100
        generate_table("our_files/%s" % file_name, table_name, is_our_table=True)



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