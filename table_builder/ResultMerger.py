__author__ = 'amirbar'

import MySQLdb
from DBCrosser import select_cross

class Counter():
    first = 0
    second = 0

def get_ecocyc_id_by_name(name_list, cursor):

    ecocyc_ids = []

    for name in name_list:

        query = """SELECT DISTINCT gene_unique_id
        FROM ecocyc_genes
        WHERE name='%s'""" % name.split(".")[0]

        cursor.execute(query)

        results = cursor.fetchallDict()

        if len(results) != 1:
            print name
            raise BaseException("too many results for one gene")


        ecocyc_ids.append(results[0]["gene_unique_id"])

    return ecocyc_ids


def get_tables_names(table_name_list, cursor):

    union_statement = ""

    for table_name in table_name_list[1:]:
        union_statement += """ UNION
    SELECT rna1_name, first_type, rna1_ecocyc_id
    FROM %(table_name)s
    UNION
    SELECT rna2_name, second_type, rna2_ecocyc_id
    FROM %(table_name)s""" % {"table_name": table_name}

    query = """SELECT rna1_name, first_type, rna1_ecocyc_id
    FROM %(table_name)s
    UNION
    SELECT rna2_name, second_type, rna2_ecocyc_id
    FROM %(table_name)s
    %(union_statement)s
    ORDER BY first_type
    """ % {"table_name": table_name_list[0],
           "union_statement": union_statement}

    cursor.execute(query)

    row = cursor.fetchone()

    names = []
    ecocyc_ids = []

    while row is not None:
        if row["rna1_name"] not in names:
            names.append(row["rna1_name"])
            ecocyc_ids.append(row["rna1_ecocyc_id"])

        row = cursor.fetchone()

    print "*" * 100
    print len(names)
    print "*" * 100

    return names, ecocyc_ids


def update_names_found(row, names_found):

    name_list = []

    if row["match_first"]:
        name_list.append(row["rna1_name"])

    if row["match_second"]:
        name_list.append(row["rna2_name"])

    if not row["match_first"] and not row["match_second"]:
        if row["closest_1"] < row["closest_2"]:
            name_list.append(row["rna1_name"])
        else:
            name_list.append(row["rna2_name"])

    for name in name_list:
        if name not in names_found.keys():
            names_found[name] = True


def merge_results(our_tables_list, their_tables_list, treshold):

    db = MySQLdb.connect(host="localhost", user="amirbar", db="amir")
    cursor = db.cursor(MySQLdb.cursors.DictCursor)

    total_names, ecocyc_ids = get_tables_names(our_tables_list, cursor)
    # ecocyc_ids = get_ecocyc_id_by_name(total_names, cursor)
    values = [[name, ecocyc_id] for name, ecocyc_id in zip(total_names, ecocyc_ids)]

    header = ["name", "ecocyc_id"]

    # Go over their tables
    for their_table in their_tables_list:

        # Go over the conditions
        for our_table in our_tables_list:

            header.append("%s_%s" % (their_table, our_table))

            names_found = {}
            select_cross(our_table, their_table, treshold, cursor)

            # Go over the cross results and add the matched rows
            row = cursor.fetchone()

            while row is not None:
                update_names_found(row, names_found)

                row = cursor.fetchone()

            for index, name in enumerate(total_names):
                if name in names_found.keys():
                    values[index].append('+')
                else:
                    values[index].append('-')

    interactions = count_interactions(our_tables, total_names, cursor)

    header.extend(interactions[interactions.keys()[0]].keys())

    for row in values:

        name = row[0]
        row.extend(interactions[name].values())

    return header, values


def get_name_count(name, table, cursor):

    query = """SELECT DISTINCT
    rna1_name, rna2_name
    FROM %(table_name)s
    WHERE rna1_name='%(name)s'""" % {"table_name": table,
                                     "name": name}

    cursor.execute(query)
    first_count = len(cursor.fetchallDict())

    query = """SELECT
    rna1_name, rna2_name, interactions
    FROM %(table_name)s
    WHERE rna1_name='%(name)s'""" % {"table_name": table,
                                     "name": name}

    cursor.execute(query)
    all_entries = cursor.fetchallDict()
    result = len(all_entries)

    if result != first_count:
        Counter.first += 1
        print "[warning] overlapped %d rows for %s first interaction under %s" % \
              (result - first_count, name, table), first_count

    total_interactions_as_first = 0

    # count the total interactions
    for entry in all_entries:
        total_interactions_as_first += int(entry["interactions"])


    query = """SELECT DISTINCT
    rna1_name, rna2_name
    FROM %(table_name)s
    WHERE rna2_name='%(name)s'""" % {"table_name": table,
                                     "name": name}

    cursor.execute(query)
    second_count = len(cursor.fetchallDict())

    query = """SELECT
    rna1_name, rna2_name, interactions
    FROM %(table_name)s
    WHERE rna2_name='%(name)s'""" % {"table_name": table,
                                     "name": name}

    cursor.execute(query)
    all_entries = cursor.fetchallDict()
    result = len(all_entries)

    if result != second_count:
        Counter.second += 1
        print "[warning] overlapped %d rows for %s second interaction under %s" % \
              (result - second_count, name, table), second_count

    total_interactions_as_second = 0

    # count the total interactions
    for entry in all_entries:
        total_interactions_as_second += int(entry["interactions"])

    return first_count, second_count, total_interactions_as_first, total_interactions_as_second


def count_interactions(our_tables_list, name_list, cursor):

    interactions_per_name = {}

    # Go over each entry in the name list
    for name in name_list:

        interactions = {}

        # Check the interactions per condition
        for condition in our_tables_list:

            interactions["%s_first_count" % condition], interactions["%s_second_count" % condition],\
                interactions["%s_first_interactions" % condition], interactions["%s_second_interactions" % condition] = \
                get_name_count(name, condition, cursor)

            interactions["%s_total" % condition] = \
                interactions["%s_first_count" % condition] + interactions["%s_second_count" % condition]

        interactions_per_name[name] = interactions

    return interactions_per_name

our_tables = ["signif_chimeras_of_iron_limitation_cl",
              "signif_chimeras_of_log_phase_cl",
              "signif_chimeras_of_stationary_cl"]

conditions = ["iron", "log_phase", "stationary"]

their_tables = ["raghavan_s5",
                "raghavan_s6",
                "raghavan_s7",
                "lybecker_s1",
                "lybecker_s2",
                "bilusic_s1",
                "bilusic_s2",
                "bilusic_s3",
                "bilusic_s4",
                "zhang_s3_2013_sheet2008",
                "zhang_s3_2013_sheet2009",
                "zhang_s4_2013_sheet2008",
                "zhang_s4_2013_sheet2009",
                "tss"]

# db = MySQLdb.connect(host="localhost", user="amirbar", db="amir")
# cursor = db.cursor(MySQLdb.cursors.DictCursor)
#
#
# names = get_tables_names(["signif_chimeras_of_iron_limitation_cl"], cursor)
#
# print len(names)
#
# dct = {}
#
# for name in names:
#     if name in dct.keys():
#         print name
#
#     dct[name] = ""
#
# print len(dct)

# interactions = count_interactions(our_tables, names, cursor)
#
# for key, value in interactions.items():
#     print key, value
#
# print len(interactions)

header, final_table = merge_results(our_tables, their_tables, 50)

fl = open("test.table", "wb")

fl.write("%s\n" % "\t".join(header))

for row in final_table:
    fl.write("%s\n" % "\t".join(str(val) for val in row))

fl.close()

print "*" * 100
print Counter.first
print Counter.second

# print header
# for row in final_table:
#     print row

# db = MySQLdb.connect(host="localhost", user="amirbar", db="amir")
# cursor = db.cursor(MySQLdb.cursors.DictCursor)
#
# tables_1 = ["all_chimeras_of_iron_limitation_cl",
#                   "all_chimeras_of_log_phase_cl",
#                   "all_chimeras_of_mg_hfq_wt101",
#                   "all_chimeras_of_mg_hfq_wt202_cl_stationary",
#                   "all_chimeras_of_stationary_cl"]
#
# tables_2 = ["single_of_iron_limitation_cl",
#             "single_of_log_phase_cl",
#             "single_of_mg_hfq_wt101",
#             "single_of_mg_hfq_wt202_cl_stationary",
#             "single_of_stationary_cl"]
#
# tables_1.extend(tables_2)
#
# names, ids = get_tables_names(tables_1, cursor)
