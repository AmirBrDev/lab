__author__ = 'amirbar'

import MySQLdb
from TableFormatter import TableFormatter


def select_cross(our_table, their_table, treshold, cursor, specific_condition="1=1"):

    match_first = """(((their.start_1 BETWEEN our.start_1 and our.end_1) OR
    (their.end_1 BETWEEN our.start_1 and our.end_1) OR
    (our.start_1 BETWEEN their.start_1 and their.end_1) OR
    (our.end_1 BETWEEN their.start_1 and their.end_1)) AND
    (their.strand_1=our.strand_1 OR their.strand_1='none' OR our.strand_1='none'))"""

    match_second = """(((their.start_2 BETWEEN our.start_2 and our.end_2) OR
    (their.end_2 BETWEEN our.start_2 and our.end_2) OR
    (our.start_2 BETWEEN their.start_2 and their.end_2) OR
    (our.end_2 BETWEEN their.start_2 and their.end_2)) AND
    (their.strand_2=our.strand_2 OR their.strand_2='none' OR our.strand_2='none'))"""

    distance_1 = """(their.start_1 - our.start_1) as distance_1"""
    distance_2 = """(their.start_2 - our.start_2) as distance_2"""

    closest_1 = """LEAST(ABS(our.start_1 - their.end_1), ABS(their.start_1 - our.end_1))"""
    closest_2 = """LEAST(ABS(our.start_2 - their.end_2), ABS(their.start_2 - our.end_2))"""
    closest = "LEAST(%s, %s)" % (closest_1, closest_2)

    in_range_first = """(((their.start_1 BETWEEN (our.start_1 - %(treshold)s) and (our.end_1 + %(treshold)s)) OR
    (their.end_1 BETWEEN (our.start_1 - %(treshold)s) and (our.end_1 + %(treshold)s)) OR
    (our.start_1 BETWEEN (their.start_1 - %(treshold)s) and (their.end_1 + %(treshold)s)) OR
    (our.end_1 BETWEEN (their.start_1 - %(treshold)s) and (their.end_1 + %(treshold)s))) AND
    (their.strand_1=our.strand_1 OR their.strand_1='none' OR our.strand_1='none'))""" % {"treshold": treshold}

    in_range_second = """(((their.start_2 BETWEEN (our.start_2 - %(treshold)s) and (our.end_2 + %(treshold)s)) OR
    (their.end_2 BETWEEN (our.start_2 - %(treshold)s) and (our.end_2 + %(treshold)s)) OR
    (our.start_2 BETWEEN (their.start_2 - %(treshold)s) and (their.end_2 + %(treshold)s)) OR
    (our.end_2 BETWEEN (their.start_2 - %(treshold)s) and (their.end_2 + %(treshold)s))) AND
    (their.strand_2=our.strand_2 OR their.strand_2='none' OR our.strand_2='none'))""" % {"treshold": treshold}

    query = """SELECT
    their.name, their.start_1, their.end_1, their.strand_1, our.rna1_name, our.start_1, our.end_1, our.strand_1,
    our.rna2_name, our.start_2, our.end_2, our.strand_2, our.rna1_ecocyc_id, our.rna2_ecocyc_id, our.first_type,
    our.second_type, our.odds_ratio, our.interactions, %(match_first)s as match_first, %(match_second)s as match_second,
    %(distance_1)s, %(distance_2)s, %(closest)s as closest, %(closest_1)s as closest_1, %(closest_2)s as closest_2,
    '%(their_table)s' as table_name
    FROM
    %(their_table)s as their, %(our_table)s as our
    WHERE
    (%(in_range_first)s OR %(in_range_second)s) AND %(specific_condition)s;""" % {"our_table": our_table,
                                                                                  "their_table": their_table,
                                                                                  "treshold": treshold,
                                                                                  "match_first": match_first,
                                                                                  "match_second": match_second,
                                                                                  "in_range_first": in_range_first,
                                                                                  "in_range_second": in_range_second,
                                                                                  "distance_1": distance_1,
                                                                                  "distance_2": distance_2,
                                                                                  "closest": closest,
                                                                                  "closest_1": closest_1,
                                                                                  "closest_2": closest_2,
                                                                                  "specific_condition": specific_condition}

    cursor.execute(query)


def select_defenetive_cross(our_table, their_table, treshold, cursor, specific_condition="1=1"):

    match_first = """(
    (
    (their.start_1 BETWEEN our.start_1 and our.end_1) OR
    (their.end_1 BETWEEN our.start_1 and our.end_1) OR
    (our.start_1 BETWEEN their.start_1 and their.end_1) OR
    (our.end_1 BETWEEN their.start_1 and their.end_1)
    ) AND
    (their.strand_1=our.strand_1 AND their.strand_1!='none' AND our.strand_1!='none')
    )"""

    match_second = """(
    (
    (their.start_2 BETWEEN our.start_2 and our.end_2) OR
    (their.end_2 BETWEEN our.start_2 and our.end_2) OR
    (our.start_2 BETWEEN their.start_2 and their.end_2) OR
    (our.end_2 BETWEEN their.start_2 and their.end_2)
    ) AND
    (their.strand_2=our.strand_2 AND their.strand_2!='none' AND our.strand_2!='none')
    )"""

    distance_1 = """(their.start_1 - our.start_1) as distance_1"""
    distance_2 = """(their.start_2 - our.start_2) as distance_2"""

    closest_1 = """LEAST(ABS(our.start_1 - their.end_1), ABS(their.start_1 - our.end_1))"""
    closest_2 = """LEAST(ABS(our.start_2 - their.end_2), ABS(their.start_2 - our.end_2))"""
    closest = "LEAST(%s, %s)" % (closest_1, closest_2)

    in_range_first = """(
    (
    (their.start_1 BETWEEN (our.start_1 - %(treshold)s) and (our.end_1 + %(treshold)s)) OR
    (their.end_1 BETWEEN (our.start_1 - %(treshold)s) and (our.end_1 + %(treshold)s)) OR
    (our.start_1 BETWEEN (their.start_1 - %(treshold)s) and (their.end_1 + %(treshold)s)) OR
    (our.end_1 BETWEEN (their.start_1 - %(treshold)s) and (their.end_1 + %(treshold)s))
    ) AND
    (their.strand_1=our.strand_1 AND their.strand_1!='none' AND our.strand_1!='none')
    )""" % {"treshold": treshold}

    in_range_second = """(
    (
    (their.start_2 BETWEEN (our.start_2 - %(treshold)s) and (our.end_2 + %(treshold)s)) OR
    (their.end_2 BETWEEN (our.start_2 - %(treshold)s) and (our.end_2 + %(treshold)s)) OR
    (our.start_2 BETWEEN (their.start_2 - %(treshold)s) and (their.end_2 + %(treshold)s)) OR
    (our.end_2 BETWEEN (their.start_2 - %(treshold)s) and (their.end_2 + %(treshold)s))
    ) AND
    (their.strand_2=our.strand_2 AND their.strand_2!='none' AND our.strand_2!='none')
    )""" % {"treshold": treshold}

    query = """SELECT
    their.name, their.start_1, their.end_1, their.strand_1, our.rna1_name, our.start_1, our.end_1, our.strand_1,
    our.rna2_name, our.start_2, our.end_2, our.strand_2, our.rna1_ecocyc_id, our.rna2_ecocyc_id, our.first_type,
    our.second_type, our.odds_ratio, our.interactions, %(match_first)s as match_first, %(match_second)s as match_second,
    %(distance_1)s, %(distance_2)s, %(closest)s as closest, %(closest_1)s as closest_1, %(closest_2)s as closest_2,
    '%(their_table)s' as table_name
    FROM
    %(their_table)s as their, %(our_table)s as our
    WHERE
    (%(in_range_first)s OR %(in_range_second)s) AND %(specific_condition)s;""" % {"our_table": our_table,
                                                                                  "their_table": their_table,
                                                                                  "treshold": treshold,
                                                                                  "match_first": match_first,
                                                                                  "match_second": match_second,
                                                                                  "in_range_first": in_range_first,
                                                                                  "in_range_second": in_range_second,
                                                                                  "distance_1": distance_1,
                                                                                  "distance_2": distance_2,
                                                                                  "closest": closest,
                                                                                  "closest_1": closest_1,
                                                                                  "closest_2": closest_2,
                                                                                  "specific_condition": specific_condition}

    cursor.execute(query)


def cross_tables(our_table, their_table, format_func, treshold, output):

    db = MySQLdb.connect(host="localhost", user="amirbar", db="amir")
    cur = db.cursor(MySQLdb.cursors.DictCursor)

    select_defenetive_cross(our_table, their_table[0], treshold, cur, specific_condition=their_table[1])

    format_func(cur, output)


def format_our_their(db_cursor, output):
    table_keys = ["rna1_name",
                  "rna2_name",
                  "rna1_ecocyc_id",
                  "rna2_ecocyc_id",
                  "first_type",
                  "second_type",
                  #result
                  "name",  # our
                  #"our.strand_1",  # our
                  #"our.strand_2",  # our
                  # distance
                  # match to
                  "interactions",
                  "odds_ratio",
                  "match_first",
                  "match_second"]

    fl = open(output, "wb")

    fl.write("%s\n" % "\t".join(table_keys))

    # keep going over the results until done
    row = db_cursor.fetchone()

    # print row.keys()

    while row is not None:

        values = [str(row[key]) for key in table_keys]
        fl.write("%s\n" % "\t".join(values))
        print values

        row = db_cursor.fetchone()


def format_their_our(db_cursor, output):

    table_keys = ["name",
                  "start_1",
                  "end_1",
                  "strand_1",
                  "rna1_name",
                  "rna2_name",
                  "rna1_ecocyc_id",
                  "rna2_ecocyc_id",
                  "first_type",
                  "second_type",
                  "interactions",
                  "odds_ratio",
                  "match_first",
                  "match_second"]

    fl = open(output, "wb")

    fl.write("%s\n" % "\t".join(table_keys))

    # keep going over the results until done
    row = db_cursor.fetchone()

    print row.keys()

    while row is not None:

        values = [str(row[key]) for key in table_keys]
        fl.write("%s\n" % "\t".join(values))
        print values

        row = db_cursor.fetchone()

    fl.close()


def format_our_single(db_cursor, output):

    header = ["our_name",
              "ecocyc_id",
              "type",
              "our_start",
              "our_end",
              "result",
              "their_name",
              "their_start",
              "their_end",
              "their_strand",
              "distance",
              "interactions",
              "odds_ratio"]

    fl = open(output, "wb")

    fl.write("%s\n" % "\t".join(header))

    # keep going over the results until done
    row = db_cursor.fetchone()

    if row is not None:
        print row.keys()

    names_to_data = {}

    while row is not None:

        values_list = []

        if row["match_first"]:
            values = [row["rna1_name"],
                      row["rna1_ecocyc_id"],
                      row["first_type"],
                      row["our.start_1"],
                      row["our.end_1"],
                      "true",
                      row["name"],
                      row["start_1"],
                      row["end_1"],
                      row["strand_1"],
                      row["distance_1"],
                      row["interactions"],
                      row["odds_ratio"]]
            values_list.append(values)

        if row["match_second"]:
            values = [row["rna2_name"],
                      row["rna2_ecocyc_id"],
                      row["second_type"],
                      row["start_2"],
                      row["end_2"],
                      "true",
                      row["name"],
                      row["start_1"],
                      row["end_1"],
                      row["strand_2"],
                      row["distance_2"],
                      row["interactions"],
                      row["odds_ratio"]]
            values_list.append(values)

        if not row["match_first"] and not row["match_second"]:
            if row["closest_1"] < row["closest_2"]:
                values = [row["rna1_name"],
                          row["rna1_ecocyc_id"],
                          row["first_type"],
                          row["our.start_1"],
                          row["our.end_1"],
                          "false",
                          row["name"],
                          row["start_1"],
                          row["end_1"],
                          row["strand_1"],
                          row["closest_1"],
                          row["interactions"],
                          row["odds_ratio"]]
            else:
                values = [row["rna2_name"],
                          row["rna2_ecocyc_id"],
                          row["second_type"],
                          row["start_2"],
                          row["end_2"],
                          "false",
                          row["name"],
                          row["start_1"],
                          row["end_1"],
                          row["strand_2"],
                          row["closest_2"],
                          row["interactions"],
                          row["odds_ratio"]]

            values_list.append(values)

        for curr_values in values_list:
            TableFormatter._add_names_to_data(names_to_data, curr_values)
            fl.write("%s\n" % "\t".join(str(val) for val in curr_values))
        # print values

        row = db_cursor.fetchone()

    fl.close()

    TableFormatter._sort_names_to_data(names_to_data)

    fl = open("%s.sorted" % output, "wb")

    fl.write("%s\n" % "\t".join(header))

    for curr in TableFormatter._sort_names_to_data(names_to_data):
        fl.write("%s\n" % "\t".join(str(val) for val in curr))

    fl.close()
