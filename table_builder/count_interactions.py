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


    db = MySQLdb.connect(host="localhost", user="amirbar", db="RILseq_seperated_annotations")
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

    table_list = [("our_files/refactor/MG_hfq-FLAG101_A_T1_7m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "101_all_interactions"),
                  ("our_files/refactor/MG_hfq-FLAG101_A_T1_7m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "101_only_singles"),
                  ("our_files/refactor/MG_hfq-FLAG102_A_T1_15m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "102_all_interactions"),
                  ("our_files/refactor/MG_hfq-FLAG102_A_T1_15m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "102_only_singles"),
                  ("our_files/refactor/MG_hfq-FLAG103_A_T1_5m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "103_all_interactions"),
                  ("our_files/refactor/MG_hfq-FLAG103_A_T1_5m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "103_only_singles"),
                  ("our_files/refactor/MG_hfq-FLAG104_A_T1_3m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "104_all_interactions"),
                  ("our_files/refactor/MG_hfq-FLAG104_A_T1_3m_22C_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "104_only_singles"),
                  ("our_files/refactor/MG_hfq-FLAG107_TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "107_all_interactions"),
                  ("our_files/refactor/MG_hfq-FLAG107_TAP_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "107_only_singles"),
                  ("our_files/refactor/MG_hfq-FLAG108_IP_50ul_Beads_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "108_all_interactions"),
                  ("our_files/refactor/MG_hfq-FLAG108_IP_50ul_Beads_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "108_only_singles"),
                  ("our_files/refactor/MG_hfq-FLAG109_lysed_by_MM400_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "109_all_interactions"),
                  ("our_files/refactor/MG_hfq-FLAG109_lysed_by_MM400_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "109_only_singles"),
                  ("our_files/refactor/MG_hfq-FLAG207_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "207_all_interactions"),
                  ("our_files/refactor/MG_hfq-FLAG207_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "207_only_singles"),
                  ("our_files/refactor/MG_hfq-FLAG208_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "208_all_interactions"),
                  ("our_files/refactor/MG_hfq-FLAG208_CL_Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "208_only_singles"),
                  ("our_files/refactor/MG_hfq-FLAG209_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "209_all_interactions"),
                  ("our_files/refactor/MG_hfq-FLAG209_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "209_only_singles"),
                  ("our_files/refactor/MG_hfq-FLAG210_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "210_all_interactions"),
                  ("our_files/refactor/MG_hfq-FLAG210_CL_Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "210_only_singles"),
                  ("our_files/refactor/MG_hfq-FLAG305-CL-Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "305_all_interactions"),
                  ("our_files/refactor/MG_hfq-FLAG305-CL-Iron_limitation_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "305_only_singles"),
                  ("our_files/refactor/MG_hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "312_all_interactions"),
                  ("our_files/refactor/MG_hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_single_counts.txt", "312_only_singles")]
                  

#    table_list = [("our_files/refactor/Total_MG-hfq-FLAG101-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "101_total_all_interactions"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG101-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "101_total_only_singles"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG102-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "102_total_all_interactions"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG102-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "102_total_only_singles"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG103-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "103_total_all_interactions"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG103-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "103_total_only_singles"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG104-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "104_total_all_interactions"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG104-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "104_total_only_singles"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG107-CL-Log-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "107_total_all_interactions"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG107-CL-Log-TAP_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "107_total_only_singles"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG108-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "108_total_all_interactions"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG108-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "108_total_only_singles"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG109-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "109_total_all_interactions"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG109-CL-Log_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "109_total_only_singles"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG207-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "207_total_all_interactions"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG207-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "207_total_only_singles"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG208-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "208_total_all_interactions"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG208-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "208_total_only_singles"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG211-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "211_total_all_interactions"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG211-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "211_total_only_singles"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG212-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "212_total_all_interactions"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG212-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "212_total_only_singles"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG305-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "305_total_all_interactions"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG305-CL-Iron_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "305_total_only_singles"),
#                  ("our_files/refactor/Total_MG-hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt", "312_total_all_interactions"),
                  #("our_files/refactor/Total_MG-hfq-FLAG312-CL-Stationary_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt", "312_total_only_singles")]

    for our_file, table_name in table_list:
        print "*" * 100
        print our_file
        upload_table(our_file, table_name)

upload_tables()

def test():

    conditions_tables = {101: ["101_total_all_interactions", "101_total_only_singles"],
                         102: ["102_total_all_interactions", "102_total_only_singles"],
                         103: ["103_total_all_interactions", "103_total_only_singles"],
                         104: ["104_total_all_interactions", "104_total_only_singles"],
                         107: ["107_total_all_interactions", "107_total_only_singles"],
                         108: ["108_total_all_interactions", "108_total_only_singles"],
                         109: ["109_total_all_interactions", "109_total_only_singles"],
                         207: ["207_total_all_interactions", "207_total_only_singles"],
                         208: ["208_total_all_interactions", "208_total_only_singles"],
                         211: ["211_total_all_interactions", "211_total_only_singles"],
                         212: ["212_total_all_interactions", "212_total_only_singles"],
                         305: ["305_total_all_interactions", "305_total_only_singles"],
                         312: ["312_total_all_interactions", "312_total_only_singles"]}

    our_tables = [table_name for table_list in conditions_tables.values() for table_name in table_list ]


    db = MySQLdb.connect(host="localhost", user="amirbar", db="RILseq_seperated_annotations")
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

    header = ["name",
              "101_total_all_interactions", "101_total_only_singles",
              "102_total_all_interactions", "102_total_only_singles",
              "103_total_all_interactions", "103_total_only_singles",
              "104_total_all_interactions", "104_total_only_singles",
              "107_total_all_interactions", "107_total_only_singles",
              "108_total_all_interactions", "108_total_only_singles",
              "109_total_all_interactions", "109_total_only_singles",
              "207_total_all_interactions", "207_total_only_singles",
              "208_total_all_interactions", "208_total_only_singles",
              "211_total_all_interactions", "211_total_only_singles",
              "211_total_all_interactions", "212_total_only_singles",
              "305_total_all_interactions", "305_total_only_singles",
              "312_total_all_interactions", "312_total_only_singles"]

    row_list = [[name] for name in short_name_list]

    # Fill the table condition by condition
    for condition_name in header[1:]:

        # Fill according to the name order
        for row in row_list:
            row.append(conditions_counts[int(condition_name[:3])][condition_name][row[0]])


    with open("count_results/counts.csv", 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(header)

        for row in row_list:
            writer.writerow(row)


# test()

