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


# generate_table("final_format/s5_directed.table", "raghavan_s5")
# generate_table("final_format/s6_directed.table", "raghavan_s6")
# generate_table("final_format/s7_directed.table", "raghavan_s7")
# generate_table("final_format/s8_directed.table", "raghavan_s8")
# generate_table("final_format/2_directed.table", "raghavan_2")

# generate_table("bilusic/final/2014RNABIOL0069R_TableS1.table", "bilusic_s1")
# generate_table("bilusic/final/2014RNABIOL0069R_TableS2.table", "bilusic_s2")
# generate_table("bilusic/final/2014RNABIOL0069R_TableS3.table", "bilusic_s3")
# generate_table("bilusic/final/2014RNABIOL0069R_TableS3_1.table", "bilusic_s3_1")
# generate_table("bilusic/final/2014RNABIOL0069R_TableS3_2.table", "bilusic_s3_2")
# generate_table("bilusic/final/2014RNABIOL0069R_TableS4.table", "bilusic_s4")
# generate_table("bilusic/final/2014RNABIOL0069R_TableS4_1.table", "bilusic_s4_1")
# generate_table("bilusic/final/2014RNABIOL0069R_TableS4_2.table", "bilusic_s4_2")

# generate_table("zhang/final/Table-s3-zhang-2013-sheet2008.table", "zhang_s3_2013_sheet2008")
# generate_table("zhang/final/Table-s3-zhang-2013-sheet2009.table", "zhang_s3_2013_sheet2009")
# generate_table("zhang/final/Table-s4-zhang-2013-sheet2008.table", "zhang_s4_2013_sheet2008")
# generate_table("zhang/final/Table-s4-zhang-2013-sheet2009.table", "zhang_s4_2013_sheet2009")

# generate_table("lybecker/final/updated_sd01.table", "lybecker_s1")
# generate_table("lybecker/final/updated_lybecker_s2.table", "lybecker_s2")

# generate_table("tss/final/JB.02096-14_zjb999093409sd1-3_with_dup.table", "thomason_full")
# generate_table("tss/final/thomason_primary.table", "thomason_primary")
# generate_table("tss/final/thomason_secondary.table", "thomason_secondary")
# generate_table("tss/final/thomason_internal.table", "thomason_internal")
# generate_table("tss/final/thomason_antisense.table", "thomason_antisense")
# generate_table("tss/final/thomason_putative_asrna.table", "thomason_putative_asrna")

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

    db = MySQLdb.connect(host="localhost",user="amirbar",db="amir")
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

def add_field_to_our_tables(table_name, source_file):

    loader = OurTableLoader()
    result = loader.loadUnprocessed(source_file)

    db = MySQLdb.connect(host="localhost",user="amirbar",db="amir")
    cur = db.cursor(MySQLdb.cursors.DictCursor)

    cur.execute("select * from %s" % table_name)

    for entry in cur.fetchallDict():

        print "*" * 30
        print entry["rna1_name"], entry["rna2_name"]

        i = 0

        for row in result:
            if row["start_1"] == entry["start_1"] and row["start_2"] == entry["start_2"] and \
                row["end_1"] == entry["end_1"] and row["end_2"] == entry["end_2"]:

                i += 1

                cur.execute("""UPDATE %(table_name)s
                SET rna1_poly_u=%(rna1_poly_u)s, rna2_poly_u=%(rna2_poly_u)s
                WHERE
                start_1=%(start_1)s and start_2=%(start_2)s and
                end_1=%(end_1)s and end_2=%(end_2)s""" % {"table_name": table_name,
                                                          "rna1_poly_u": db.literal(row["poly-u of rna1"]),
                                                          "rna2_poly_u": db.literal(row["poly-u of rna2"]),
                                                          "start_1": row["start_1"],
                                                          "start_2": row["start_2"],
                                                          "end_1": row["end_1"],
                                                          "end_2": row["end_2"]})

                matches = cur.fetchallDict()

        if i != 1:
            print "--> error!"

        print i

    db.commit()


# add_field_to_our_tables("signif_chimeras_of_iron_limitation_cl", "our_with_poly_u/Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_sig_interactions.txt_with_pU.txt")
#
# add_field_to_our_tables("signif_chimeras_of_stationary_cl", "our_with_poly_u/Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_sig_interactions.txt_with_pU.txt")
#
# add_field_to_our_tables("signif_chimeras_of_log_phase_cl", "our_with_poly_u/Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_sig_interactions.txt_with_pU.txt")


def generate_signif_chimeras_name_table():
    db = MySQLdb.connect(host="localhost",user="amirbar",db="amir")
    cur = db.cursor(MySQLdb.cursors.DictCursor)

    cur.execute("""CREATE TABLE signif_chimeras (name VARCHAR(200), type VARCHAR(200))""")

    tables = ["signif_chimeras_of_iron_limitation_cl",
              "signif_chimeras_of_log_phase_cl",
              "signif_chimeras_of_stationary_cl"]

    union_statement = ""

    for table_name in tables[1:]:
        union_statement += """ UNION
    SELECT rna1_name, first_type
    FROM %(table_name)s
    UNION
    SELECT rna2_name, second_type
    FROM %(table_name)s""" % {"table_name": table_name}

    query = """SELECT rna1_name, first_type
    FROM %(table_name)s
    UNION
    SELECT rna2_name, second_type
    FROM %(table_name)s
    %(union_statement)s
    ORDER BY first_type
    """ % {"table_name": tables[0],
           "union_statement": union_statement}

    cur.execute(query)

    dct = cur.fetchallDict()

    for row in dct:
        print row
        cur.execute("INSERT INTO signif_chimeras VALUES (%s, %s)" % \
                    (db.literal(row["rna1_name"]), db.literal(row["first_type"])))

    db.commit()

# generate_signif_chimeras_name_table()

def generate_csv_table(path, name):

    db = MySQLdb.connect(host="localhost",user="amirbar",db="amir")
    cur = db.cursor(MySQLdb.cursors.DictCursor)

    fl = open(path, "rb")

    header = fl.readline().lower().replace("\n", "").split("\t")

    key_list = []

    for key in header:
        if key == "condition":
            key_list.append("_condition")
        else:
            key_list.append(key)

    fields = ", ".join("%s VARCHAR(200)" % key.replace(" ", "_").replace("-", "_").replace("?", "") for key in key_list)
    print len(key_list)

    print "CREATE TABLE %s (%s)" % (name, fields)
    cur.execute("CREATE TABLE %s (%s)" % (name, fields))

    for row in fl.readlines():
        values = ", ".join(db.literal(str(val)) for val in row.lower().replace("\n", "").split("\t"))
        cur.execute("INSERT INTO %s VALUES (%s)" % (name, values))

    db.commit()

    fl.close()


# generate_csv_table("MEME_results.csv", "meme_results")


def add_type_to_meme_results():

    db = MySQLdb.connect(host="localhost",user="amirbar",db="amir")
    cur = db.cursor(MySQLdb.cursors.DictCursor)

    cur.execute("""SELECT meme.*, genes.type
    FROM meme_results as meme, signif_chimeras as genes
    WHERE meme.gene_name = genes.name""")

    rows = cur.fetchallDict()

    fl = open("MEME_results.csv", "rb")

    header = fl.readline().lower().replace("\n", "").split("\t")

    fl.close()

    key_list = []

    for key in header:
        if key == "condition":
            key_list.append("_condition")
        else:
            key_list.append(key)

    fields = [key.replace(" ", "_").replace("-", "_").replace("?", "") for key in key_list]

    fields.append("type")

    fl = open("MEME_results_with_types.csv", "wb")

    fl.write("%s\n" % "\t".join(fields))

    for row in rows:

        to_write = [row[key] for key in fields]

        print row
        fl.write("%s\n" % "\t".join(to_write))

    fl.close()


# add_type_to_meme_results()

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


def get_zhang_stats_by_name(name):

    db = MySQLdb.connect(host="localhost",user="amirbar",db="amir")
    cur = db.cursor()

    # query = """SELECT k31a_ip_k31a_total, r16a_ip_____r16a_total, q8a_ip_____q8a_total,
    # k31a_total_wt_total, r16a_total____wt_total, q8a_total___wt_total
    # FROM zhang_s4_2013_sheet2009 as zhang
    # WHERE zhang.name=%(name)s;""" % {"name": db.literal(name)}

    query = """SELECT k31a_ip____k31a_total, r16a_ip______r16a_total, q8a_ip________q8a_total,
    k31a_total__wt_total, r16a_total__wt_total, q8a_total__wt_total
    FROM zhang_s3_2013_sheet2009 as zhang
    WHERE zhang.name=%(name)s;""" % {"name": db.literal(name)}
    cur.execute(query)

    return cur.fetchall()


def generate_zhang_stats():
    loader = TableLoader()

    rows_as_dictionary = loader.load("output/table_s6_range_0.csv")

    row_list = []

    header = ["name",
              "il_rna2_percent",
              "stat_rna2_percent",
              "log_rna2_percent",
              "k31_ip", #  distal
              "r16a_ip", #  rim
              "q8a_ip", #  proximal
              "k31",
              "r16a",
              "q8a",
              "average_percent"]

    for row in rows_as_dictionary:

        new_row = [row["name"],
                   row["signif_chimeras_of_iron_limitation_cl.as_rna2_percentage"].replace("-", ""),
                   row["signif_chimeras_of_stationary_cl.as_rna2_percentage"].replace("-", ""),
                   row["signif_chimeras_of_log_phase_cl.as_rna2_percentage"].replace("-", "")]

        matches = get_zhang_stats_by_name(row["name"])

        if len(matches) > 1:
            print "warning too many results"

        elif len(matches) == 1:
            new_row.extend(val for val in matches[0])
            # print new_row

        else:
            new_row.extend([""] * 6)

        average_percent = 0.0
        fields = 0

        if row["signif_chimeras_of_iron_limitation_cl.as_rna2_percentage"] != "-":
            average_percent += float(row["signif_chimeras_of_iron_limitation_cl.as_rna2_percentage"])
            fields += 1

        if row["signif_chimeras_of_stationary_cl.as_rna2_percentage"] != "-":
            average_percent += float(row["signif_chimeras_of_stationary_cl.as_rna2_percentage"])
            fields += 1

        if row["signif_chimeras_of_log_phase_cl.as_rna2_percentage"] != "-":
            average_percent += float(row["signif_chimeras_of_log_phase_cl.as_rna2_percentage"])
            fields += 1

        average_percent /= fields

        new_row.append(average_percent)

        row_list.append(new_row)

    # for row in row_list:
    #     print row

    fl = open("zhang_stats.csv", "wb")

    fl.write("%s\n" % "\t".join(header))

    for row in row_list:
        fl.write("%s\n" % "\t".join(str(val) for val in row))

    fl.close()

# generate_zhang_stats()


def find_all_combos(table_name_list, cursor, allow_5utr_dup):

    combos = []
    union_statement = ""

    for table_name in table_name_list[1:]:
        union_statement += """ UNION
        SELECT rna1_name, rna2_name
        FROM %(table_name)s""" % {"table_name": table_name}

    query = """SELECT rna1_name, rna2_name
    FROM %(table_name)s
    %(union_statement)s""" % {"table_name": table_name_list[0],
                              "union_statement": union_statement}
    cursor.execute(query)

    row = cursor.fetchone()

    while row is not None:
        
        first_name = row["rna1_name"]
        second_name = row["rna2_name"]

        if not allow_5utr_dup:
            first_name = first_name.replace(".5utr", "").replace(".est5utr", "")
            second_name = second_name.replace(".5utr", "").replace(".est5utr", "")

        if (first_name, second_name) not in combos and \
           (second_name, first_name) not in combos:
            combos.append((first_name, second_name))

        row = cursor.fetchone()

    return combos


def find_all_srna_combos(table_name_list, cursor, allow_5utr_dup):

    combos = []
    union_statement = ""

    for table_name in table_name_list[1:]:
        union_statement += """ UNION
        SELECT rna1_name, rna2_name
        FROM %(table_name)s
        WHERE first_type = 'srna' OR second_type='srna'""" % {"table_name": table_name}

    query = """SELECT rna1_name, rna2_name
    FROM %(table_name)s
    WHERE first_type = 'srna' OR second_type='srna'
    %(union_statement)s""" % {"table_name": table_name_list[0],
                              "union_statement": union_statement}
    cursor.execute(query)

    row = cursor.fetchone()

    while row is not None:
        
        first_name = row["rna1_name"]
        second_name = row["rna2_name"]

        if not allow_5utr_dup:
            first_name = first_name.replace(".5utr", "").replace(".est5utr", "")
            second_name = second_name.replace(".5utr", "").replace(".est5utr", "")

        if (first_name, second_name) not in combos and \
           (second_name, first_name) not in combos:
            combos.append((first_name, second_name))

        row = cursor.fetchone()

    return combos


def test_find_all_combos():
    db = MySQLdb.connect(host="localhost",user="amirbar",db="amir")
    cursor = db.cursor(MySQLdb.cursors.DictCursor)

    tables = ["signif_chimeras_of_iron_limitation_cl",
              "signif_chimeras_of_log_phase_cl",
              "signif_chimeras_of_stationary_cl"]    

    print len(find_all_combos(tables, cursor, True))
    print len(find_all_combos(tables, cursor, False))
    print len(find_all_srna_combos(tables, cursor, True))
    print len(find_all_srna_combos(tables, cursor, False))

#test_find_all_combos()

def find_all_srna_combos_when_srna_second(table_name_list, cursor, allow_5utr_dup):

    combos = []
    union_statement = ""

    for table_name in table_name_list[1:]:
        union_statement += """ UNION
        SELECT rna1_name, rna2_name
        FROM %(table_name)s
        WHERE second_type = 'srna' 
        UNION
        SELECT rna2_name, rna1_name
        FROM %(table_name)s
        WHERE first_type = 'srna'""" % {"table_name": table_name}

    query = """SELECT rna1_name, rna2_name
    FROM %(table_name)s
    WHERE second_type = 'srna' 
    UNION
    SELECT rna2_name, rna1_name
    FROM %(table_name)s
    WHERE first_type = 'srna' 
    %(union_statement)s""" % {"table_name": table_name_list[0],
                              "union_statement": union_statement}
    cursor.execute(query)

    row = cursor.fetchone()

    while row is not None:
        
        first_name = row["rna1_name"]
        second_name = row["rna2_name"]

        if not allow_5utr_dup:
            first_name = first_name.replace(".5utr", "").replace(".est5utr", "")
            second_name = second_name.replace(".5utr", "").replace(".est5utr", "")

        if (first_name, second_name) not in combos and \
           (second_name, first_name) not in combos:
            combos.append((first_name, second_name))

        row = cursor.fetchone()

    return combos

def test_srna_combos_2nd():
    db = MySQLdb.connect(host="localhost",user="amirbar",db="amir")
    cursor = db.cursor(MySQLdb.cursors.DictCursor)

    tables = ["signif_chimeras_of_iron_limitation_cl",
              "signif_chimeras_of_log_phase_cl",
              "signif_chimeras_of_stationary_cl"]

    #combos = find_all_srna_combos_when_srna_second(tables, cursor, False)
    combos = find_all_combos(tables, cursor, False)

    for first, second in combos:
        print "\t".join([first, second])

#test_srna_combos_2nd()

def fix_thomason_csv_cols():
    result = []

    with open("tss/csv/JB.02096-14_zjb999093409sd1-3.csv", "rb") as fl:

        header = fl.readline().lower().replace("\n", "").split(Table.TABLE_DELIMITER)

        name_list = []

        for key in header:

            if key == "Sequence -50 nt upstream + TSS (51nt)".lower():
                name_list.append("_sequence")

            elif key == "primary":
                name_list.append("_primary")

            elif key == "secondary":
                name_list.append("_secondary")

            elif key == "condition":
                name_list.append("_condition")

            elif key == "putative srna":
                name_list.append("putative_srna")

            elif key == "putative asrna":
                name_list.append("putative_asrna")

            elif key == "Overlap with RegulonDB".lower():
                name_list.append("overlap_with_regulondb")

            else:
                name_list.append(key)

        for line in fl.readlines():

            line = line.lower()

            dct = {}

            value_list = line.split(Table.TABLE_DELIMITER)

            if (len(name_list) != len(value_list)):
                raise BaseException("Unmatched arguments count in table.")

            for key, val in zip(name_list, value_list):
                dct[key.replace("\n", "")] = val.replace("\n", "")

            result.append(dct)

    with open("tss/csv/JB.02096-14_zjb999093409sd1-3_fixed_cols.csv", "wb") as fl:
        header = ["start_1", "end_1", "strand_1", "start_2", "end_2", "strand_2"]
        header += name_list[2:]
        fl.write("\t".join(key for key in header) + "\n")

        for row in result:
            value_list = [row["pos"],
                          row["pos"],
                          row["strand"],
                          row["pos"],
                          row["pos"],
                          row["strand"]]

            for key in name_list[2:]:
                value_list.append(row[key])

            fl.write("\t".join(val for val in value_list) + "\n")

fix_thomason_csv_cols()

def upload_thomason_full_table():

    TABLE_NAME = "thomason_full"
    result = []


    with open("tss/csv/JB.02096-14_zjb999093409sd1-3_fixed_cols.csv", "rb") as fl:

        name_list = fl.readline().lower().replace("\n", "").split(Table.TABLE_DELIMITER)

        for line in fl.readlines():

            line = line.lower()

            dct = {}

            value_list = line.split(Table.TABLE_DELIMITER)

            if (len(name_list) != len(value_list)):
                raise BaseException("Unmatched arguments count in table.")

            for key, val in zip(name_list, value_list):
                dct[key.replace("\n", "")] = val.replace("\n", "")

            result.append(dct)

    db = MySQLdb.connect(host="localhost", user="amirbar", db="amir")
    cursor = db.cursor(MySQLdb.cursors.DictCursor)


    # Remove the current table
    query="""DROP TABLE IF EXISTS %s""" % TABLE_NAME
    cursor.execute(query)

    # Recreate it with the new fields
    fields = ", ".join("%s VARCHAR(200)" % key for key in name_list)
    print "CREATE TABLE %s (%s)" % (TABLE_NAME, fields)
    cursor.execute("CREATE TABLE %s (%s)" % (TABLE_NAME, fields))

    # Recreate it with the new fields
    for row in result:
        values = ",".join("%s" % db.literal(str(row[key])) for key in name_list)
        cursor.execute("INSERT INTO %s VALUES (%s)" % (TABLE_NAME, values))

    db.commit()

upload_thomason_full_table()