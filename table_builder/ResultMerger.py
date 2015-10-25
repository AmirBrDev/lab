__author__ = 'amirbar'

import MySQLdb
from DBCrosser import select_cross, select_defenetive_cross
from TableLoader import TableLoader

our_tables = ["signif_chimeras_of_iron_limitation_cl",
              "signif_chimeras_of_log_phase_cl",
              "signif_chimeras_of_stationary_cl"]

conditions = ["iron", "log_phase", "stationary"]

their_tables = [("raghavan_s5", "mev >= 1"),
                ("raghavan_s6", "1=1"),
                ("raghavan_s7", "mev >= 1"),
                ("raghavan_2", "mev >= 1"),
                # "raghavan_s8",
                ("lybecker_s1", "1=1"),
                ("lybecker_s2", "1=1"),
                ("bilusic_s1", "1=1"),
                ("bilusic_s2", "1=1"),
                ("bilusic_s3_1", "1=1"),
                ("bilusic_s3_2", "1=1"),
                ("bilusic_s4_1", "1=1"),
                ("bilusic_s4_2", "1="),
                # "zhang_s3_2013_sheet2008",
                # "zhang_s3_2013_sheet2009",
                # "zhang_s4_2013_sheet2008",
                # "zhang_s4_2013_sheet2009",
                ("thomason", "1=1"),
                ("thomason_primary", "1=1"),
                ("thomason_secondary", "1=1"),
                ("thomason_internal", "1=1"),
                ("thomason_antisense", "1=1"),
                ("thomason_putative_asrna", "1=1")]

our_file_list = ["our_files/assign-type-to-signif-chimeras-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_sig_interactions.with-type",
                 "our_files/assign-type-to-signif-chimeras-of-Log_phase_CL_FLAG101-104_108_109_all_fragments_l25.txt_sig_interactions.with-type",
                 "our_files/assign-type-to-signif-chimeras-of-Stationary_CL_FLAG209_210_312_all_fragments_l25.txt_sig_interactions.with-type"]


class Counter():

    first = 0
    second = 0

    def __init__(self):
        pass


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
    types = []

    while row is not None:
        if row["rna1_name"] not in names:
            names.append(row["rna1_name"])
            ecocyc_ids.append(row["rna1_ecocyc_id"])
            types.append(row["first_type"])

        row = cursor.fetchone()

    print "*" * 100
    print len(names)
    print "*" * 100

    return names, ecocyc_ids, types


def get_tables_srna_names(table_name_list, cursor):

    union_statement = ""

    for table_name in table_name_list[1:]:
        union_statement += """ UNION
    SELECT rna1_name, first_type
    FROM %(table_name)s
    WHERE first_type='srna'
    UNION
    SELECT rna2_name, second_type
    FROM %(table_name)s
    WHERE second_type='srna'""" % {"table_name": table_name}

    query = """SELECT rna1_name, first_type
    FROM %(table_name)s
    WHERE first_type='srna'
    UNION
    SELECT rna2_name, second_type
    FROM %(table_name)s
    WHERE second_type='srna'
    %(union_statement)s""" % {"table_name": table_name_list[0],
                              "union_statement": union_statement}

    cursor.execute(query)

    row = cursor.fetchone()

    names = []

    while row is not None:
        if row["rna1_name"] not in names:
            names.append(row["rna1_name"])

        row = cursor.fetchone()

    return names


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

    total_names, ecocyc_ids, types = get_tables_names(our_tables_list, cursor)
    # ecocyc_ids = get_ecocyc_id_by_name(total_names, cursor)
    values = [[name, ecocyc_id, type] for name, ecocyc_id, type in zip(total_names, ecocyc_ids, types)]

    header = ["name", "ecocyc_id", "type"]

    # Generate other article hit table

    # Go over their tables
    for their_table in their_tables_list:

        # Go over the conditions
        for our_table in our_tables_list:

            header.append("%s_%s" % (their_table[0], our_table))

            names_found = {}

            # Change this to allow none strands (remove the comment from the first)
            # select_cross(our_table, their_table, treshold, cursor)
            select_defenetive_cross(our_table, their_table[0], treshold, cursor, their_table[1])

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

    # Generate interactions count and polyU length
    interactions = count_interactions(our_tables, total_names, cursor)

    header.extend(interactions[interactions.keys()[0]].keys())
    header.append("max_poly_u_length")

    for row in values:

        name = row[0]
        row.extend(interactions[name].values())
        row.append(get_poly_u_for_name(name, our_tables_list, cursor))

    # Generate target columns
    conditions_counts, ignore = get_target_counts(total_names, our_tables, cursor)

    for condition, names_to_counts in conditions_counts.items():

        header.append(condition)

        for row in values:

            name = row[0]

            if name in names_to_counts.keys():
                row.append(names_to_counts[name])
            else:
                row.append(0)

    total_targets = get_total_targets(total_names, our_tables, True, cursor)
    header.append("total_targets")

    for row in values:
        name = row[0]
        row.append(total_targets[name])

    # Generate meme, mast columns
    meme_mast_values = get_meme_mast_values(total_names, cursor)

    header.append("meme")
    header.append("mast")
    header.append("motif")
    header.append("meme_results")
    header.append("binding_site_state")

    for row in values:
        name = row[0]
        meme_result = ""
        mast_result = ""
        motif = ""
        binding_site = "-"

        is_reciprocal = False

        # Select best meme/mast values for each entry
        for meme_val, mast_val, motif_val, binding_site_val in meme_mast_values[name]:
            if meme_val is not None and \
               mast_val is not None and \
               meme_val <= 0.05 and mast_val <= 0.05:

                # Update first time
                if not is_reciprocal:
                    meme_result = "%.2e" % meme_val
                    mast_result = "%.2e" % mast_val
                    motif = motif_val
                    binding_site = binding_site_val
                    is_reciprocal = True

                # Update if same binding and improves results
                elif binding_site == binding_site_val:

                    if mast_val < float(mast_result):
                        meme_result = "%.2e" % meme_val
                        mast_result = "%.2e" % mast_val
                        motif = motif_val
                        binding_site = binding_site_val

                    elif mast_val == float(mast_result) and \
                         meme_val < float(meme_result):

                        meme_result = "%.2e" % meme_val
                        mast_result = "%.2e" % mast_val
                        motif = motif_val
                        binding_site = binding_site_val

                # Update if improves binding
                elif binding_site_val == "yes" or \
                     (binding_site_val == "no" and binding_site == "-"):

                    meme_result = "%.2e" % meme_val
                    mast_result = "%.2e" % mast_val
                    motif = motif_val
                    binding_site = binding_site_val


            elif not is_reciprocal:
                meme_result = ";".join([val for val in ["%.2e" % meme_val, meme_result] if val != ""])
                mast_result = ";".join([val for val in ["%.2e" % mast_val, mast_result] if val != ""])

        row.append(meme_result)
        row.append(mast_result)
        row.append(motif)

        if is_reciprocal:
            row.append("recip")
        else:
            row.append("")

        row.append(binding_site)

    # Generate tb-srnas targets
    header.append("tb_srna_targets")
    names_to_tb_srna_targets = find_tbs_srna_targets(total_names, our_tables, cursor)

    for row in values:
        name = row[0]

        if name in names_to_tb_srna_targets.keys():
            row.append(len(names_to_tb_srna_targets[name]))

    # Generate mrna/5utr targets
    header.append("mrna_5utr_targets")
    names_to_mrna_5utr_targets = find_mrna_5utr_targets(total_names, our_tables, cursor)

    for row in values:
        name = row[0]

        if name in names_to_mrna_5utr_targets.keys():
            row.append(len(names_to_mrna_5utr_targets[name]))

    # Generate igr/3utr targets
    header.append("igr_3utr_targets")
    names_to_igr_3utr_targets = find_igr_3utr_targets(total_names, our_tables, cursor)

    for row in values:
        name = row[0]

        if name in names_to_igr_3utr_targets.keys():
            row.append(len(names_to_igr_3utr_targets[name]))

    # Generate consistency - Just for my curiosity
    # header.append("total_as_first")
    # header.append("unique_as_first")
    # header.append("as_first_consistency_precentage")
    #
    # for row in values:
    #     name = row[0]
    #
    #     total, unique, percentage = get_consistency_as_first(our_tables_list, name, cursor)
    #     row.append(total)
    #     row.append(unique)
    #     row.append(percentage)

    return header, values


def get_poly_u_for_name(name, table_name_list, cursor):
    union_statement = ""

    for table_name in table_name_list[1:]:
        union_statement += """ UNION
        SELECT rna1_name, rna1_poly_u
        FROM %(table_name)s
        WHERE rna1_name = '%(name)s'
        UNION
        SELECT rna2_name, rna2_poly_u
        FROM %(table_name)s
        WHERE rna2_name = '%(name)s'""" % {"table_name": table_name,
                                           "name": name}

    query = """SELECT unified.rna1_name, MAX(unified.rna1_poly_u)
    FROM (SELECT rna1_name, rna1_poly_u
    FROM %(table_name)s
    WHERE rna1_name = '%(name)s'
    UNION
    SELECT rna2_name, rna2_poly_u
    FROM %(table_name)s
    WHERE rna2_name = '%(name)s'
    %(union_statement)s) as unified""" % {"table_name": table_name_list[0],
                                          "name": name,
                                          "union_statement": union_statement}

    cursor.execute(query)
    row = cursor.fetchone()

    if row is not None:
        return row["MAX(unified.rna1_poly_u)"]
    else:
        print "[Error] Couldn't find poly_u for %s" % name


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


def find_targets(names, table_name_list, cursor):

    names_to_targets = {}

    for name in names:

        union_statement = ""

        for table_name in table_name_list[1:]:
            union_statement += """ UNION
            SELECT rna1_name
            FROM %(table_name)s
            WHERE rna2_name = '%(name)s'
            UNION
            SELECT rna2_name
            FROM %(table_name)s
            WHERE rna1_name = '%(name)s'""" % {"table_name": table_name,
                                               "name": name}

        query = """SELECT rna1_name
        FROM %(table_name)s
        WHERE rna2_name = '%(name)s'
        UNION
        SELECT rna2_name
        FROM %(table_name)s
        WHERE rna1_name = '%(name)s'
        %(union_statement)s""" % {"table_name": table_name_list[0],
                                  "name": name,
                                  "union_statement": union_statement}

        cursor.execute(query)

        row = cursor.fetchone()

        names_to_targets[name] = []

        while row is not None:
            if row["rna1_name"] not in names_to_targets[name]:
                names_to_targets[name].append(row["rna1_name"])

            row = cursor.fetchone()

    return names_to_targets


def find_tbs_srna_targets(names, table_name_list, cursor):

    names_to_targets = {}

    for name in names:

        union_statement = ""

        for table_name in table_name_list[1:]:
            union_statement += """ UNION
            SELECT rna1_name
            FROM %(table_name)s
            WHERE rna2_name = '%(name)s' AND first_type='srna'
            UNION
            SELECT rna2_name
            FROM %(table_name)s
            WHERE rna1_name = '%(name)s' AND second_type='srna'""" % {"table_name": table_name,
                                                                      "name": name}

        query = """SELECT rna1_name
        FROM %(table_name)s
        WHERE rna2_name = '%(name)s' AND first_type='srna'
        UNION
        SELECT rna2_name
        FROM %(table_name)s
        WHERE rna1_name = '%(name)s' AND second_type='srna'
        %(union_statement)s""" % {"table_name": table_name_list[0],
                                  "name": name,
                                  "union_statement": union_statement}

        cursor.execute(query)

        row = cursor.fetchone()

        names_to_targets[name] = []

        while row is not None:

            to_add = row["rna1_name"].replace(".5utr", "").replace(".est5utr", "")

            if to_add not in names_to_targets[name]:
                names_to_targets[name].append(to_add)

            row = cursor.fetchone()

    return names_to_targets


def find_mrna_5utr_targets(names, table_name_list, cursor):

    names_to_targets = {}

    for name in names:

        union_statement = ""

        for table_name in table_name_list[1:]:
            union_statement += """ UNION
            SELECT rna1_name
            FROM %(table_name)s
            WHERE rna2_name = '%(name)s' AND (first_type='mrna' OR first_type='5utr')
            UNION
            SELECT rna2_name
            FROM %(table_name)s
            WHERE rna1_name = '%(name)s' AND (second_type='mrna' OR second_type='5utr')""" % {"table_name": table_name,
                                                                                              "name": name}

        query = """SELECT rna1_name
        FROM %(table_name)s
        WHERE rna2_name = '%(name)s' AND (first_type='mrna' OR first_type='5utr')
        UNION
        SELECT rna2_name
        FROM %(table_name)s
        WHERE rna1_name = '%(name)s' AND (second_type='mrna' OR second_type='5utr')
        %(union_statement)s""" % {"table_name": table_name_list[0],
                                  "name": name,
                                  "union_statement": union_statement}

        cursor.execute(query)

        row = cursor.fetchone()

        names_to_targets[name] = []

        while row is not None:

            to_add = row["rna1_name"].replace(".5utr", "").replace(".est5utr", "")

            if to_add not in names_to_targets[name]:
                names_to_targets[name].append(to_add)

            row = cursor.fetchone()

    return names_to_targets


def find_igr_3utr_targets(names, table_name_list, cursor):

    names_to_targets = {}

    for name in names:

        union_statement = ""

        for table_name in table_name_list[1:]:
            union_statement += """ UNION
            SELECT rna1_name
            FROM %(table_name)s
            WHERE rna2_name = '%(name)s' AND (first_type='igr' OR first_type='3utr')
            UNION
            SELECT rna2_name
            FROM %(table_name)s
            WHERE rna1_name = '%(name)s' AND (second_type='igr' OR second_type='3utr')""" % {"table_name": table_name,
                                                                                              "name": name}

        query = """SELECT rna1_name
        FROM %(table_name)s
        WHERE rna2_name = '%(name)s' AND (first_type='igr' OR first_type='3utr')
        UNION
        SELECT rna2_name
        FROM %(table_name)s
        WHERE rna1_name = '%(name)s' AND (second_type='igr' OR second_type='3utr')
        %(union_statement)s""" % {"table_name": table_name_list[0],
                                  "name": name,
                                  "union_statement": union_statement}

        cursor.execute(query)

        row = cursor.fetchone()

        names_to_targets[name] = []

        while row is not None:

            to_add = row["rna1_name"]

            if to_add not in names_to_targets[name]:
                names_to_targets[name].append(to_add)

            row = cursor.fetchone()

    return names_to_targets


def get_combinations(targets_dictionary):
    combinations = []

    for name, target_list in targets_dictionary.items():
        for target in target_list:

            target_name = target.replace(".5utr", "").replace(".est5utr", "")

            if (name, target_name) not in combinations and (target_name, name) not in combinations:
                combinations.append((name, target_name))

    return combinations


def get_total_targets(names, table_name_list, merge_5utr, cursor):
    targets_dictionary = find_targets(names, table_name_list, cursor)

    for name, target_list in targets_dictionary.items():

        if merge_5utr:
            targets_dictionary[name] = \
                len(set(target.replace(".5utr", "").replace(".est5utr", "") for target in target_list))

        else:
            targets_dictionary[name] = len(set(target_list))

        if targets_dictionary[name] == 0:
            print "[warning] no targets! -", name, targets_dictionary[name]

    return targets_dictionary


def get_target_counts(names, table_name_list, cursor):

    conditions = {}
    conditions_combinations = {}

    for table in table_name_list:
        targets_dictionary = find_targets(names, [table], cursor)

        names_to_counts = {}

        # get the length for each name ignoring the .5utr
        for name, target_list in targets_dictionary.items():
            names_to_counts[name] = \
                len(set(target.replace(".5utr", "").replace(".est5utr", "") for target in target_list))

        combinations = get_combinations(targets_dictionary)

        conditions["%s_targets" % table] = names_to_counts
        conditions_combinations["%s_targets" % table] = combinations

    return conditions, conditions_combinations


def test_counts(our_tables):

    db = MySQLdb.connect(host="localhost", user="amirbar", db="amir")
    cursor = db.cursor(MySQLdb.cursors.DictCursor)

    conditions_counts, condtion_combinations = \
        get_target_counts(get_tables_srna_names(our_tables, cursor), our_tables, cursor)

    for condition, names_to_counts in conditions_counts.items():

        print condition
        print "*" * 30

        for name, count in names_to_counts.items():
            print "%s: %d" % (name, count)

    all_combs = []

    print condtion_combinations.keys()
    arr = condtion_combinations.values()

    for comb in arr:
        all_combs.extend([pair for pair in comb if pair not in all_combs])

    buckets = [0] * 8

    for pair in all_combs:

        bucket = 0

        if pair in arr[0]:
            bucket += 4

        if pair in arr[1]:
            bucket += 2

        if pair in arr[2]:
            bucket += 1

        buckets[bucket] += 1

    for index, count in enumerate(buckets):

        print format(index, "03b"), count

    print buckets


# test_counts(our_tables)

def get_meme_mast_values(names, cursor):

    query = """SELECT meme.gene_name, meme.meme_e_value, meme.mast_e_value, meme.meme_motif, meme.match_known_bs
    FROM meme_results as meme, signif_chimeras
    WHERE meme.gene_name = signif_chimeras.name and meme.number_of_targets >= 5"""

    cursor.execute(query)
    result = cursor.fetchallDict()

    names_to_values = {name: [] for name in names}

    for row in result:
        names_to_values[row["gene_name"]].append((row["meme_e_value"],
                                                  row["mast_e_value"],
                                                  row["meme_motif"],
                                                  row["match_known_bs"]))

    return names_to_values


def test_meme_values():
    db = MySQLdb.connect(host="localhost", user="amirbar", db="amir")
    cursor = db.cursor(MySQLdb.cursors.DictCursor)

    our_tables = ["signif_chimeras_of_iron_limitation_cl",
                  "signif_chimeras_of_log_phase_cl",
                  "signif_chimeras_of_stationary_cl"]

    total_names = get_tables_names(our_tables, cursor)[0]

    result = get_meme_mast_values(total_names, cursor)

    for name, values in result.items():
        print name, values


# test_meme_values()

def test_get_poly_u_by_name(name, tables_list, cursor):
    print "%s poly u is %s" % (name, str(get_poly_u_for_name(name, tables_list, cursor)))


# test_get_poly_u_by_name("yiep", our_tables, cursor)

def test_interactions(our_tables_list):
    db = MySQLdb.connect(host="localhost", user="amirbar", db="amir")
    cursor = db.cursor(MySQLdb.cursors.DictCursor)

    names = get_tables_names(["signif_chimeras_of_iron_limitation_cl"], cursor)[0]

    print len(names)

    dct = {}

    for name in names:
        if name in dct.keys():
            print name

        dct[name] = ""

    print len(dct)

    interactions = count_interactions(our_tables_list, names, cursor)

    for key, value in interactions.items():
        print key, value

    print len(interactions)


# test_interactions(our_tables)

def test_merge_results(treshold):
    header, final_table = merge_results(our_tables, their_tables, treshold)

    fl = open("test.csv", "wb")

    fl.write("%s\n" % "\t".join(header))

    for row in final_table:
        fl.write("%s\n" % "\t".join(str(val) for val in row))

    fl.close()

    print "*" * 100
    print Counter.first
    print Counter.second


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


def get_ecocyc_id_dictionary(our_file_list):
    result = {}

    for file_path in our_file_list:
        fl = open(file_path, "rb")

        fl.readline()
        for line in fl.readlines():
            args = line.split("\t")
            id_1, id_2= args[2], args[3]

            if id_1.lower() not in result.keys():
                result[id_1.lower()] = id_1

            if id_2.lower() not in result.keys():
                result[id_2.lower()] = id_2

    return result


def get_consistency_as_first(table_name_list, name, cursor):

    union_statement = ""

    for table_name in table_name_list[1:]:
        union_statement += """ UNION
        SELECT rna2_name
        FROM %(table_name)s
        WHERE rna1_name = '%(name)s'""" % {"table_name": table_name,
                                           "name": name}

    query = """SELECT rna2_name
    FROM %(table_name)s
    WHERE rna1_name = '%(name)s'
    %(union_statement)s""" % {"table_name": table_name_list[0],
                              "name": name,
                              "union_statement": union_statement}

    cursor.execute(query)

    total = len(cursor.fetchallDict())

    union_statement = ""

    for table_name in table_name_list[1:]:
        union_statement += """ UNION
        SELECT rna2_name
        FROM %(table_name)s
        WHERE rna1_name = '%(name)s' AND
        rna2_name NOT IN
        (SELECT rna1_name FROM %(table_name)s WHERE rna2_name='%(name)s')""" % {"table_name": table_name,
                                                                                "name": name}

    query = """SELECT rna2_name
        FROM %(table_name)s
        WHERE rna1_name = '%(name)s' AND
        rna2_name NOT IN
        (SELECT rna1_name FROM %(table_name)s WHERE rna2_name='%(name)s')
        %(union_statement)s""" % {"table_name": table_name,
                                  "name": name,
                                  "union_statement": union_statement}

    cursor.execute(query)

    unique_first = len(cursor.fetchallDict())

    if total != 0:
        percentage = float(unique_first) / float(total)

    else:
        percentage = '-'

    return total, unique_first, percentage

def get_condition_by_header_name(field, conditions, conditions_beauty_names):
    
    for index, name in enumerate(conditions_beauty_names):
        if field.endswith(name):
            return conditions[index]

    return None

def format_final_table(path, our_tables):

    sets = {"Raghavan et al 2011": ["raghavan_s5", "raghavan_s6", "raghavan_s7", "raghavan_2"],
            "Lybecker et al 2014": ["lybecker_s1", "lybecker_s2"],
            "Bilusic et al 2014": ["bilusic_s1", "bilusic_s2", "bilusic_s3_1", "bilusic_s3_2", "bilusic_s4_1", "bilusic_s4_2"],
            # "zhang": ["zhang_s3_2013_sheet2008", "zhang_s3_2013_sheet2009", "zhang_s4_2013_sheet2008", "zhang_s4_2013_sheet2009"],
            "Thomason et al 2015": ["thomason", "thomason_primary", "thomason_secondary", "thomason_internal", "thomason_antisense", "thomason_putative_asrna"]}

    conditions = ["signif_chimeras_of_iron_limitation_cl",
                  "signif_chimeras_of_log_phase_cl",
                  "signif_chimeras_of_stationary_cl"]

    conditions_beauty_names = ["Iron limitation",
                               "Log",
                               "Stationary"]

    short_name = {"raghavan_s5": "R1",
                  "raghavan_s6": "R2",
                  "raghavan_s7": "R3",
                  "raghavan_2": "R4",
                  # "raghavan_s8": "R4",
                  "lybecker_s1": "L1",
                  "lybecker_s2": "L2",
                  "bilusic_s1": "B1",
                  "bilusic_s2": "B2",
                  "bilusic_s3_1": "B3_1",
                  "bilusic_s3_2": "B3_2",
                  "bilusic_s4_1": "B4_1",
                  "bilusic_s4_2": "B4_2",
                  # "zhang_s3_2013_sheet2008": "Z1",
                  # "zhang_s3_2013_sheet2009": "Z2",
                  # "zhang_s4_2013_sheet2008": "Z3",
                  # "zhang_s4_2013_sheet2009": "Z4",
                  "thomason": "T1",
                  "thomason_primary": "T1_1",
                  "thomason_secondary": "T1_2",
                  "thomason_internal": "T1_3",
                  "thomason_antisense": "T1_4",
                  "thomason_putative_asrna": "T1_5"}

    beauty_type_names = {"3utr": "3UTR",
                         "5utr": "5UTR",
                         "as": "AS",
                         "cis_as_with_trans_t": "cASt",
                         "igr": "IGR",
                         "mrna": "CDS",
                         "other-ncrna": "oRNA",
                         "srna": "TB-sRNA",
                         "trna": "tRNA",
                         "tu": "IGT"}

    loader = TableLoader()
    results = loader.load(path)

    header = ["Name",
              "EcoCyc id",
              "Type",
              "Total UI",
              "TB-sRNA UI",
              "CDS & 5'UTR UI",
              "3'UTR & IGR UI"]

    start_of_total_interactions = len(header)
    for cond_name in conditions_beauty_names:
        header.append("TNR %s" % cond_name)

    end_of_total_interactions = len(header)

    start_of_interactions = end_of_total_interactions

    for cond_name in conditions_beauty_names:
        header.append("Fraction as RNA2 %s" % cond_name)

    end_of_interactions = len(header)

    header.extend(["Longest poly U",
                   "MEME E-value",
                   "MAST E-value",
                   "Meme motif",
                   "Overlaps known binding site"])

    start_of_tables = len(header)

    for set_name in sets:
        header.append(set_name)

    end_of_tables = len(header)

    header.append("Paper number")

    final_rows = []

    # Go over the rows and fill according to the header
    for index, row in enumerate(results):
        row_values = [row["name"],
                      row["ecocyc_id"],
                      row["type"],
                      row["total_targets"],
                      row["tb_srna_targets"],
                      row["mrna_5utr_targets"],
                      row["igr_3utr_targets"]]

        # Total interactions
        for field in header[start_of_total_interactions: end_of_total_interactions]:
            cond_name = get_condition_by_header_name(field, conditions, conditions_beauty_names)

            first_count = int(row["%s_first_interactions" % cond_name])
            second_count = int(row["%s_second_interactions" % cond_name])

            row_values.append(first_count + second_count)

        # interactions percentage
        for field in header[start_of_interactions: end_of_interactions]:
            cond_name = get_condition_by_header_name(field, conditions, conditions_beauty_names)

            first_count = float(row["%s_first_interactions" % cond_name])
            second_count = float(row["%s_second_interactions" % cond_name])
            total_count = first_count + second_count

            if float(total_count) == 0:
                res = "-"
            else:
                res = second_count / total_count
                res = "%.2f" % res

            row_values.append(res)

        row_values.extend([row["max_poly_u_length"],
                           row["meme"].upper(),
                           row["mast"].upper(),
                           row["motif"],
                           row["binding_site_state"]])

        total_articles = 0

        # Go over the table hit fields and merge columns
        for set_name in header[start_of_tables: end_of_tables]:

            field_values = []

            for cond_name in conditions:
                for table in sets[set_name]:
                    if row["%s_%s" % (table, cond_name)] == "+":
                        field_values.append(short_name[table])
                    elif row["%s_%s" % (table, cond_name)] == "-":
                        continue
                    else:
                        print "[warning] invalid value for cell"

            if len(field_values) > 0:
                total_articles += 1

            row_values.append(";".join(list(set(field_values))))

        row_values.append(total_articles)

        final_rows.append(row_values)

    # go over the rows and fix name, id, type notations
    names = get_name_dictionary(our_tables)
    ecocyc_ids = get_ecocyc_id_dictionary(our_tables)

    for row in final_rows:
        row[0] = names[row[0]]
        row[1] = ecocyc_ids[row[1]]
        row[2] = beauty_type_names[row[2]]

        if row[1].split(".")[-1].isdigit():
            row[1] = ".".join(row[1].split(".")[:-1])

        if ".TU" in row[0]:
            row[0] = row[0].replace(".TU", ".IGT")

        if ".TU" in row[1]:
            row[1] = row[1].replace(".TU", ".IGT")

    fl = open("formatted_test.csv", "wb")

    fl.write("%s\n" % "\t".join(header))

    for row in final_rows:
        fl.write("%s\n" % "\t".join(str(val) for val in row))

#test_merge_results(0)
format_final_table("test.csv", our_file_list)
