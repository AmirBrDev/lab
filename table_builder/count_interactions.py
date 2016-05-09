__author__ = 'amirbar'

import MySQLdb
import csv

def get_names_from_tables(table_name_list, cursor):

    union_statement = ""

    for table_name in table_name_list[1:]:
        union_statement += """ UNION
    SELECT rna1_name
    FROM %(table_name)s
    UNION
    SELECT rna2_name
    FROM %(table_name)s""" % {"table_name": table_name}

    query = """SELECT rna1_name
    FROM %(table_name)s
    UNION
    SELECT rna2_name
    FROM %(table_name)s
    %(union_statement)s
    """ % {"table_name": table_name_list[0],
           "union_statement": union_statement}

    cursor.execute(query)

    row = cursor.fetchone()

    names = []

    while row is not None:
        if row["rna1_name"] not in names:
            names.append(row["rna1_name"])

        row = cursor.fetchone()

    print "*" * 100
    print len(names)
    print "*" * 100

    return names


def get_counts_per_name_for_table(name_list, table_name, cursor):

    short_name_list = set(name.replace(".5utr", "").replace(".est5utr", "") for name in name_list)
    name_to_interactions =  {name : 0 for name in short_name_list}

    query = """SELECT rna1_name, rna2_name, interactions
    FROM %(table_name)s
    """ % {"table_name": table_name}

    cursor.execute(query)

    row = cursor.fetchone()

    while row is not None:

        rna1_name = row["rna1_name"].replace(".5utr", "").replace(".est5utr", "")
        rna2_name = row["rna2_name"].replace(".5utr", "").replace(".est5utr", "")

        name_to_interactions[rna1_name] += int(row["interactions"])
        name_to_interactions[rna2_name] += int(row["interactions"])

        row = cursor.fetchone()

    return name_to_interactions

def upload_table(table_file, table_name):

    DELIMITER = "\t"
    result = []

    with open(table_file, "rb") as fl:

        name_list = fl.readline().lower().replace("\n", "").split(DELIMITER)
        print name_list
        for line in fl.readlines():

            line = line.lower()

            dct = {}

            value_list = line.split(DELIMITER)

            if (len(name_list) != len(value_list)):
                raise BaseException("Unmatched arguments count in table.")

            for key, val in zip(name_list, value_list):
                dct[key.replace("\n", "")] = val.replace("\n", "")

            result.append(dct)


    db = MySQLdb.connect(host="stone.md.huji.ac.il", user="amirbar", db="RILseq_seperated_annotations")
    cursor = db.cursor(MySQLdb.cursors.DictCursor)


    # Remove the current table
    query="""DROP TABLE IF EXISTS %s""" % table_name
    cursor.execute(query)

    # Recreate it with the new fields
    fields = ", ".join("%s VARCHAR(300)" % key.replace(" ", "_").replace("-", "_").replace("'", "") for key in name_list)
    print "CREATE TABLE %s (%s)" % (table_name, fields)
    cursor.execute("CREATE TABLE %s (%s)" % (table_name, fields))

    # Recreate it with the new fields
    for row in result:
        values = ",".join("%s" % db.literal(str(row[key])) for key in name_list)
        cursor.execute("INSERT INTO %s VALUES (%s)" % (table_name, values))

    cursor.execute("ALTER TABLE %s CHANGE rna1_from start_1 int not null" % table_name)
    cursor.execute("ALTER TABLE %s CHANGE rna1_to end_1 int not null" % table_name)
    cursor.execute("ALTER TABLE %s CHANGE rna1_strand strand_1 VARCHAR(300) not null" % table_name)
    cursor.execute("ALTER TABLE %s CHANGE rna2_from start_2 int not null" % table_name)
    cursor.execute("ALTER TABLE %s CHANGE rna2_to end_2 int not null" % table_name)
    cursor.execute("ALTER TABLE %s CHANGE rna2_strand strand_2 VARCHAR(300) not null" % table_name)
    cursor.execute("ALTER TABLE %s CHANGE interactions interactions int not null" % table_name)
    cursor.execute("ALTER TABLE %s CHANGE other_interactions_of_rna1 other_interactions_of_rna1 int not null" % table_name)
    cursor.execute("ALTER TABLE %s CHANGE other_interactions_of_rna2 other_interactions_of_rna2 int not null" % table_name)
    cursor.execute("ALTER TABLE %s CHANGE total_other_interactions total_other_interactions int not null" % table_name)
    cursor.execute("ALTER TABLE %s CHANGE odds_ratio odds_ratio int not null" % table_name)
    cursor.execute("ALTER TABLE %s CHANGE fishers_exact_test_p_value fishers_exact_test_p_value int not null" % table_name)

    db.commit()

def upload_tables():

    #table_list = [("our_files/refactor/myflip/MG_hfq-FLAG101_A_T1_7m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "101_all_interactions"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG101_A_T1_7m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "101_only_singles"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG102_A_T1_15m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "102_all_interactions"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG102_A_T1_15m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "102_only_singles"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG103_A_T1_5m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "103_all_interactions"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG103_A_T1_5m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "103_only_singles"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG104_A_T1_3m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "104_all_interactions"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG104_A_T1_3m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "104_only_singles"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG107_TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "107_all_interactions"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG107_TAP_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "107_only_singles"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG108_IP_50ul_Beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "108_all_interactions"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG108_IP_50ul_Beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "108_only_singles"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG109_lysed_by_MM400_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "109_all_interactions"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG109_lysed_by_MM400_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "109_only_singles"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG207_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "207_all_interactions"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG207_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "207_only_singles"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG208_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "208_all_interactions"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG208_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "208_only_singles"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG209_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "209_all_interactions"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG209_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "209_only_singles"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG210_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "210_all_interactions"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG210_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "210_only_singles"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG305-CL-Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "305_all_interactions"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG305-CL-Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "305_only_singles"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "312_all_interactions"),
    #              ("our_files/refactor/myflip/MG_hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "312_only_singles"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG101-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "101_total_all_interactions"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG101-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "101_total_only_singles"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG102-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "102_total_all_interactions"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG102-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "102_total_only_singles"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG103-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "103_total_all_interactions"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG103-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "103_total_only_singles"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG104-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "104_total_all_interactions"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG104-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "104_total_only_singles"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG107-CL-Log-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "107_total_all_interactions"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG107-CL-Log-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "107_total_only_singles"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG108-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "108_total_all_interactions"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG108-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "108_total_only_singles"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG109-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "109_total_all_interactions"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG109-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "109_total_only_singles"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG207-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "207_total_all_interactions"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG207-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "207_total_only_singles"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG208-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "208_total_all_interactions"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG208-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "208_total_only_singles"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG211-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "211_total_all_interactions"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG211-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "211_total_only_singles"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG212-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "212_total_all_interactions"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG212-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "212_total_only_singles"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG305-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "305_total_all_interactions"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG305-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "305_total_only_singles"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "312_total_all_interactions"),
    #              ("our_files/refactor/myflip/Total_MG-hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "312_total_only_singles")]

    table_list = [("our_files/refactor/set_2/myflip/MG-hfq-FLAG401-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "401_all_interactions"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG408-Log_CL-100-bella-beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "408_all_interactions"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG401-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "401_only_singles"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG408-Log_CL-100-bella-beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "408_only_singles"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG402-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "402_all_interactions"),   
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG410-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "410_all_interactions"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG402-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "402_only_singles"),           
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG410-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "410_only_singles"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG404-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "404_all_interactions"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG411-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "411_all_interactions"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG404-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "404_only_singles"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG411-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "411_only_singles"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG405-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "405_all_interactions"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG413-Log_CL_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "413_all_interactions"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG405-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "405_only_singles"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG413-Log_CL_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "413_only_singles"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG407-Log_CL-100-bella-beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "407_all_interactions"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG414-Log_CL_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "414_all_interactions"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG407-Log_CL-100-bella-beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "407_only_singles"),
                  ("our_files/refactor/set_2/myflip/MG-hfq-FLAG414-Log_CL_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "414_only_singles")]

    for our_file, table_name in table_list:
        print "*" * 100
        print our_file
        upload_table(our_file, table_name)

# upload_tables()


def flip_files_igr():

    workdir = "/home/users/amirbar/lab/table_builder/our_files/refactor/set_2"

    #file_list = ["MG_hfq-FLAG101_A_T1_7m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "Total_MG-hfq-FLAG101-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "MG_hfq-FLAG101_A_T1_7m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #             "Total_MG-hfq-FLAG101-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #             "MG_hfq-FLAG102_A_T1_15m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "Total_MG-hfq-FLAG102-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "MG_hfq-FLAG102_A_T1_15m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #             "Total_MG-hfq-FLAG102-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #             "MG_hfq-FLAG103_A_T1_5m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "Total_MG-hfq-FLAG103-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "MG_hfq-FLAG103_A_T1_5m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #             "Total_MG-hfq-FLAG103-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #             "MG_hfq-FLAG104_A_T1_3m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "Total_MG-hfq-FLAG104-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "MG_hfq-FLAG104_A_T1_3m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #             "Total_MG-hfq-FLAG104-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #             "MG_hfq-FLAG107_TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "Total_MG-hfq-FLAG107-CL-Log-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "MG_hfq-FLAG107_TAP_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #             "Total_MG-hfq-FLAG107-CL-Log-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #             "MG_hfq-FLAG108_IP_50ul_Beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "Total_MG-hfq-FLAG108-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "MG_hfq-FLAG108_IP_50ul_Beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #             "Total_MG-hfq-FLAG108-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #             "MG_hfq-FLAG109_lysed_by_MM400_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "Total_MG-hfq-FLAG109-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "MG_hfq-FLAG109_lysed_by_MM400_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #             "Total_MG-hfq-FLAG109-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #             "MG_hfq-FLAG207_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "Total_MG-hfq-FLAG207-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "MG_hfq-FLAG207_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #             "Total_MG-hfq-FLAG207-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #             "MG_hfq-FLAG208_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "Total_MG-hfq-FLAG208-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "MG_hfq-FLAG208_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #             "Total_MG-hfq-FLAG208-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #             "MG_hfq-FLAG209_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "Total_MG-hfq-FLAG211-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "MG_hfq-FLAG209_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #             "Total_MG-hfq-FLAG211-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #             "MG_hfq-FLAG210_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "Total_MG-hfq-FLAG212-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "MG_hfq-FLAG210_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #             "Total_MG-hfq-FLAG212-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #             "MG_hfq-FLAG305-CL-Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "Total_MG-hfq-FLAG305-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "MG_hfq-FLAG305-CL-Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #             "Total_MG-hfq-FLAG305-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #             "MG_hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "Total_MG-hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #             "MG_hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #             "Total_MG-hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt"]

    file_list = ["MG-hfq-FLAG401-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                  "MG-hfq-FLAG408-Log_CL-100-bella-beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                  "MG-hfq-FLAG401-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
                  "MG-hfq-FLAG408-Log_CL-100-bella-beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
                  "MG-hfq-FLAG402-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",   
                  "MG-hfq-FLAG410-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                  "MG-hfq-FLAG402-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",           
                  "MG-hfq-FLAG410-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
                  "MG-hfq-FLAG404-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",              
                  "MG-hfq-FLAG411-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                  "MG-hfq-FLAG404-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",                 
                  "MG-hfq-FLAG411-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
                  "MG-hfq-FLAG405-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",              
                  "MG-hfq-FLAG413-Log_CL_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                  "MG-hfq-FLAG405-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",                 
                  "MG-hfq-FLAG413-Log_CL_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
                  "MG-hfq-FLAG407-Log_CL-100-bella-beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                  "MG-hfq-FLAG414-Log_CL_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                  "MG-hfq-FLAG407-Log_CL-100-bella-beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
                  "MG-hfq-FLAG414-Log_CL_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt"]

    for file_path in file_list:
        print file_path
        print "*" * 100

        with open("/".join([workdir, file_path]), "rb") as in_file:
            reader = csv.reader(in_file, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            with open("/".join([workdir, "myflip", file_path]), "wb") as out_file:
                writer = csv.writer(out_file, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

                for row in reader:
                    rna1_ecocyc_id = row[0]
                    rna2_ecocyc_id = row[1]
                    rna1_name = row[2]
                    rna2_name = row[3]
                    strand_1 = row[9]
                    strand_2 = row[13]

                    if ".igr" in rna1_name.lower() and strand_1 == "-":
                        name_parts = rna1_name.split(".")
                        row[2] = ".".join([name_parts[1], name_parts[0], name_parts[2]])

                        id_parts = rna1_ecocyc_id.split(".")
                        row[0] = ".".join([id_parts[1], id_parts[0], id_parts[2]])

                    if ".igr" in rna2_name.lower() and strand_2 == "-":
                        name_parts = rna2_name.split(".")
                        row[3] = ".".join([name_parts[1], name_parts[0], name_parts[2]])

                        id_parts = rna2_ecocyc_id.split(".")
                        row[1] = ".".join([id_parts[1], id_parts[0], id_parts[2]])

                    writer.writerow(row)

# flip_files_igr()


def generate_beautiful_name_table(file_list):

    lower_to_upper_names = {}

    for file_name in file_list:

        with open(file_name, "rb") as fl:
            fl.readline()
            for line in fl.readlines():
                args = line.split("\t")
                name_1, name_2 = args[2], args[3]

                name_parts = set(name_1.split(".") + name_2.split("."))

                for part in name_parts:
                    if part.lower() not in lower_to_upper_names.keys():
                        lower_to_upper_names[part.lower()] = part

    return lower_to_upper_names


def get_name_dictionary(our_file_list):
    result = {}

    for file_path in our_file_list:
        fl = open(file_path, "rb")

        fl.readline()
        for line in fl.readlines():
            args = line.split("\t")
            name_1, name_2 = args[4], args[5]

            if name_1.lower() not in result.keys():
                result[name_1.lower()] = name_1

            if name_2.lower() not in result.keys():
                result[name_2.lower()] = name_2

    return result


def test():

    # file_list = ["our_files/refactor/MG_hfq-FLAG101_A_T1_7m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/MG_hfq-FLAG101_A_T1_7m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #              "our_files/refactor/MG_hfq-FLAG102_A_T1_15m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/MG_hfq-FLAG102_A_T1_15m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #              "our_files/refactor/MG_hfq-FLAG103_A_T1_5m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/MG_hfq-FLAG103_A_T1_5m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #              "our_files/refactor/MG_hfq-FLAG104_A_T1_3m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/MG_hfq-FLAG104_A_T1_3m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #              "our_files/refactor/MG_hfq-FLAG107_TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/MG_hfq-FLAG107_TAP_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #              "our_files/refactor/MG_hfq-FLAG108_IP_50ul_Beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/MG_hfq-FLAG108_IP_50ul_Beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #              "our_files/refactor/MG_hfq-FLAG109_lysed_by_MM400_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/MG_hfq-FLAG109_lysed_by_MM400_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #              "our_files/refactor/MG_hfq-FLAG207_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/MG_hfq-FLAG207_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #              "our_files/refactor/MG_hfq-FLAG208_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/MG_hfq-FLAG208_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #              "our_files/refactor/MG_hfq-FLAG209_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/MG_hfq-FLAG209_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #              "our_files/refactor/MG_hfq-FLAG210_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/MG_hfq-FLAG210_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #              "our_files/refactor/MG_hfq-FLAG305-CL-Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/MG_hfq-FLAG305-CL-Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #              "our_files/refactor/MG_hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/MG_hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG101-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG101-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG102-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG102-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG103-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG103-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG104-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG104-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG107-CL-Log-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG107-CL-Log-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG108-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG108-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG109-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG109-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG207-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG207-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG208-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG208-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG211-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG211-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG212-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG212-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG305-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG305-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
    #              "our_files/refactor/Total_MG-hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt"]

    workdir = "/home/users/amirbar/lab/table_builder/our_files/refactor/set_2/"

    file_list = [workdir + "MG-hfq-FLAG401-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                 workdir + "MG-hfq-FLAG408-Log_CL-100-bella-beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                 workdir + "MG-hfq-FLAG401-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
                 workdir + "MG-hfq-FLAG408-Log_CL-100-bella-beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
                 workdir + "MG-hfq-FLAG402-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                 workdir + "MG-hfq-FLAG410-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                 workdir + "MG-hfq-FLAG402-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
                 workdir + "MG-hfq-FLAG410-Log_CL-100-beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
                 workdir + "MG-hfq-FLAG404-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                 workdir + "MG-hfq-FLAG411-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                 workdir + "MG-hfq-FLAG404-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
                 workdir + "MG-hfq-FLAG411-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
                 workdir + "MG-hfq-FLAG405-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                 workdir + "MG-hfq-FLAG413-Log_CL_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                 workdir + "MG-hfq-FLAG405-Log_CL-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
                 workdir + "MG-hfq-FLAG413-Log_CL_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
                 workdir + "MG-hfq-FLAG407-Log_CL-100-bella-beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                 workdir + "MG-hfq-FLAG414-Log_CL_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt",
                 workdir + "MG-hfq-FLAG407-Log_CL-100-bella-beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt",
                 workdir + "MG-hfq-FLAG414-Log_CL_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt"]

    # conditions_tables = {101: ["101_total_all_interactions", "101_total_only_singles", "101_all_interactions", "101_only_singles"],
    #                      102: ["102_total_all_interactions", "102_total_only_singles", "102_all_interactions", "102_only_singles"],
    #                      103: ["103_total_all_interactions", "103_total_only_singles", "103_all_interactions", "103_only_singles"],
    #                      104: ["104_total_all_interactions", "104_total_only_singles", "104_all_interactions", "104_only_singles"],
    #                      107: ["107_total_all_interactions", "107_total_only_singles", "107_all_interactions", "107_only_singles"],
    #                      108: ["108_total_all_interactions", "108_total_only_singles", "108_all_interactions", "108_only_singles"],
    #                      109: ["109_total_all_interactions", "109_total_only_singles", "109_all_interactions", "109_only_singles"],
    #                      207: ["207_total_all_interactions", "207_total_only_singles", "207_all_interactions", "207_only_singles"],
    #                      208: ["208_total_all_interactions", "208_total_only_singles", "208_all_interactions", "208_only_singles"],
	# 	             	   209: ["209_all_interactions", "209_only_singles"],
    #                      210: ["210_all_interactions", "210_only_singles"],
    #                      211: ["211_total_all_interactions", "211_total_only_singles"],
    #                      212: ["212_total_all_interactions", "212_total_only_singles"],
    #                      305: ["305_total_all_interactions", "305_total_only_singles", "305_all_interactions", "305_only_singles"],
    #                      312: ["312_total_all_interactions", "312_total_only_singles", "312_all_interactions", "312_only_singles"]}
    conditions_tables = {401: ["401_all_interactions", "401_only_singles"],
                         402: ["402_all_interactions", "402_only_singles"],
                         404: ["404_all_interactions", "404_only_singles"],
                         405: ["405_all_interactions", "405_only_singles"],
                         407: ["407_all_interactions", "407_only_singles"],
                         408: ["408_all_interactions", "408_only_singles"],
                         410: ["410_all_interactions", "410_only_singles"],
                         411: ["411_all_interactions", "411_only_singles"],
                         413: ["413_all_interactions", "413_only_singles"],
			             414: ["414_all_interactions", "414_only_singles"]}


    our_tables = [table_name for table_list in conditions_tables.values() for table_name in table_list ]


    db = MySQLdb.connect(host="stone.md.huji.ac.il", user="amirbar", db="RILseq_seperated_annotations")
    cursor = db.cursor(MySQLdb.cursors.DictCursor)

    name_list = get_names_from_tables(our_tables, cursor)

    short_name_list = set(name.replace(".5utr", "").replace(".est5utr", "") for name in name_list)

    conditions_counts = {key: {} for key in conditions_tables.keys()}

    # Go over the different conditions
    for condition in conditions_tables.keys():

        # Get the counts per table
        for table in conditions_tables[condition]:
            names_to_counts = get_counts_per_name_for_table(name_list, table, cursor)

            conditions_counts[condition][table] = names_to_counts

            with open("count_results/%s_count.csv" % table, 'wb') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                writer.writerow(["name", "count"])

                for name, count in names_to_counts.items():

                    # print name, count
                    writer.writerow([name, count])



    # header = ["name",
    #           "FLAG101",
    #           "FLAG101_total",
    #           "FLAG102",
    #           "FLAG102_total",
    #           "FLAG103",
    #           "FLAG103_total",
    #           "FLAG104",
    #           "FLAG104_total",
    #           "FLAG107",
    #           "FLAG107_total",
    #           "FLAG108",
    #           "FLAG108_total",
    #           "FLAG109",
    #           "FLAG109_total",
    #           "FLAG207",
    #           "FLAG207_total",
    #           "FLAG208",
    #           "FLAG208_total",
    #           "FLAG209",
    #           "FLAG211_total",
    #           "FLAG210",
    #           "FLAG212_total",
    #           "FLAG305",
    #           "FLAG305_total",
    #           "FLAG312",
    #           "FLAG312_total"]
	#
    # conditions_list = [("101_all_interactions", "101_only_singles"),
    #                    ("101_total_all_interactions", "101_total_only_singles"),
    #                    ("102_all_interactions", "102_only_singles"),
    #                    ("102_total_all_interactions", "102_total_only_singles"),
    #                    ("103_all_interactions", "103_only_singles"),
    #                    ("103_total_all_interactions", "103_total_only_singles"),
    #                    ("104_all_interactions", "104_only_singles"),
    #                    ("104_total_all_interactions", "104_total_only_singles"),
    #                    ("107_all_interactions", "107_only_singles"),
    #                    ("107_total_all_interactions", "107_total_only_singles"),
    #                    ("108_all_interactions", "108_only_singles"),
    #                    ("108_total_all_interactions", "108_total_only_singles"),
    #                    ("109_all_interactions", "109_only_singles"),
    #                    ("109_total_all_interactions", "109_total_only_singles"),
    #                    ("207_all_interactions", "207_only_singles"),
    #                    ("207_total_all_interactions", "207_total_only_singles"),
    #                    ("208_all_interactions", "208_only_singles"),
    #                    ("208_total_all_interactions", "208_total_only_singles"),
    #                    ("209_all_interactions", "209_only_singles"),
    #                    ("211_total_all_interactions", "211_total_only_singles"),
    #                    ("210_all_interactions", "210_only_singles"),
    #                    ("212_total_all_interactions", "212_total_only_singles"),
    #                    ("305_all_interactions", "305_only_singles"),
    #                    ("305_total_all_interactions", "305_total_only_singles"),
    #                    ("312_all_interactions", "312_only_singles"),
    #                    ("312_total_all_interactions", "312_total_only_singles")]

    header = ["name",
              "FLAG401",
              "FLAG402",
              "FLAG404",
              "FLAG405",
              "FLAG407",
              "FLAG408",
              "FLAG410",
              "FLAG411",
              "FLAG413",
              "FLAG414"]

    conditions_list = [("401_all_interactions", "401_only_singles"),
                       ("402_all_interactions", "402_only_singles"),
                       ("404_all_interactions", "404_only_singles"),
                       ("405_all_interactions", "405_only_singles"),
                       ("407_all_interactions", "407_only_singles"),
                       ("408_all_interactions", "408_only_singles"),
                       ("410_all_interactions", "410_only_singles"),
                       ("411_all_interactions", "411_only_singles"),
                       ("413_all_interactions", "413_only_singles"),
                       ("414_all_interactions", "414_only_singles")]


    print "-" * 50
    print "Getting Beautiful name list..."

    beautiful_name_convert = generate_beautiful_name_table(file_list)

    print "done."

    #row_list = [[".".join(beautiful_name_convert[part] for part in name.split("."))] for name in short_name_list]
    row_list = [[name] for name in short_name_list]

    # Fill the table condition by condition
    for condition_name_1, condition_name_2 in conditions_list:

        print "-" * 50
        print "Processing conditions %s, %s" % (condition_name_1, condition_name_2)

        # Fill according to the name order
        for row in row_list:

            # Sum both tables of the same condition, row[0] is the gene name
            row.append(conditions_counts[int(condition_name_1[:3])][condition_name_1][row[0]] +
                       conditions_counts[int(condition_name_2[:3])][condition_name_2][row[0]])

        print "done."

    print "-" * 50
    print "Converting to fractions"

    # Sum totals per column
    col_totals = [sum([row[index + 1] for row in row_list]) for index in range(len(conditions_list))]

    # Update each cell to be fraction of the column total
    for index, val in enumerate(col_totals):
        for row in row_list:
            row[index + 1] = float(row[index + 1]) / val

    for index in range(1, len(header), 2):

        header.append("_".join([header[index], "fraction"]))

        for row in row_list:
            try:
                row.append(row[index] / row[index + 1])

            except ZeroDivisionError:
                row.append("N/A")



    print "done."

    print "-" * 50
    print "Generating final table"

    with open("count_results/fractions/counts.csv", 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(header)

        for row in row_list:
            row[0] = ".".join(beautiful_name_convert[part] for part in row[0].split("."))
            writer.writerow(row)


    print "finished!"

test()

