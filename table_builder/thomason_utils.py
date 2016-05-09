__author__ = 'amirbar'

THOMASON_CROSS_TABLE = "thomason_refactor.csv"

import MySQLdb
from Table import Table

class SignifTableUpdater():
    def __init__(self, table_name_list):
        self.db = MySQLdb.connect(host="localhost", user="amirbar", db="article_refactor_24_3_2016")
        self.cursor = self.db.cursor(MySQLdb.cursors.DictCursor)
        self.table_name_list = table_name_list

    def get_names_from_tables(self):
        """
        Extract a tuple of lists containing names, ecocyc_ids and types from our tables.

        :return: see above
        """

        union_statement = ""

        for table_name in self.table_name_list[1:]:
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
        """ % {"table_name": self.table_name_list[0],
               "union_statement": union_statement}

        self.cursor.execute(query)

        row = self.cursor.fetchone()

        names = []
        ecocyc_ids = []
        types = []

        while row is not None:
            if row["rna1_name"] not in names:
                names.append(row["rna1_name"])
                ecocyc_ids.append(row["rna1_ecocyc_id"])
                types.append(row["first_type"])

            row = self.cursor.fetchone()

        return names, ecocyc_ids, types

    def get_min_start_max_end(self, name):

        union_statement = ""

        for table_name in self.table_name_list[1:]:
            union_statement += """ UNION
            SELECT rna1_name, min(start_1), max(end_1), strand_1
            FROM %(table_name)s
            WHERE rna1_name='%(name)s'
            UNION
            SELECT rna2_name, min(start_2), max(end_2), strand_2
            FROM %(table_name)s
            WHERE rna2_name='%(name)s'""" % {"table_name": table_name,
                                           "name": name}

            query = """SELECT rna1_name, min(start_1), max(end_1), strand_1
            FROM %(table_name)s
            WHERE rna1_name='%(name)s'
            UNION
            SELECT rna2_name, min(start_2), max(end_2), strand_2
            FROM %(table_name)s
            WHERE rna2_name='%(name)s'
            %(union_statement)s
            """ % {"table_name": self.table_name_list[0],
                   "name": name,
                   "union_statement": union_statement}

        self.cursor.execute(query)

        row = self.cursor.fetchone()

        max_end = None
        min_start = None
        strand = None

        while row is not None:

            if row["min(start_1)"] is None or row["max(end_1)"] is None:
                row = self.cursor.fetchone()
                continue

            if max_end is None or min_start is None:
                min_start = int(row["min(start_1)"])
                max_end = int(row["max(end_1)"])
                strand = row["strand_1"]

            if int(row["min(start_1)"]) < min_start:
                min_start = int(row["min(start_1)"])

            if int(row["max(end_1)"]) > max_end:
                max_end = int(row["max(end_1)"])

            row = self.cursor.fetchone()

        return min_start, max_end, strand

    def add_min_max_limits(self):
        names_list = self.get_names_from_tables()[0]

        for name in names_list:
            min_start, max_end, strand = self.get_min_start_max_end(name)

            query = """UPDATE signif_chimeras
            SET min_start='%(min_start)s', max_end='%(max_end)s', strand='%(strand)s'
            WHERE name='%(name)s'""" % {"name": name,
                                        "min_start": min_start,
                                        "max_end": max_end,
                                        "strand": strand}
            self.cursor.execute(query)

        self.db.commit()

def test_update_signif_table():
    our_tables = ["signif_chimeras_of_iron_limitation_cl",
                      "signif_chimeras_of_log_phase_cl",
                      "signif_chimeras_of_stationary_cl"]

    SignifTableUpdater(our_tables).add_min_max_limits()

# test_update_signif_table()

def get_bounding_genes(name):

    db = MySQLdb.connect(host="localhost", user="amirbar", db="article_refactor_24_3_2016")
    cursor = db.cursor(MySQLdb.cursors.DictCursor)

    query="""SELECT *
    FROM
        (SELECT our.name as our_name, eco.name as eco_name, (our.min_start-eco.start) as dist
         FROM signif_chimeras as our, ecocyc_genes as eco
         WHERE eco.strand=our.strand and our.name='%s' order by dist) as lower
    WHERE dist > 0 limit 1""" % name

    cursor.execute(query)

    row = cursor.fetchone()

    if row is  None:
        print "[ERROR]bounding not found for %s" % name
        start_bound ="none"

    else:
        start_bound = row["eco_name"]

    query="""SELECT *
    FROM
        (SELECT our.name as our_name, eco.name as eco_name, (eco.end-our.max_end) as dist
         FROM signif_chimeras as our, ecocyc_genes as eco
         WHERE eco.strand=our.strand and our.name='%s' order by dist) as lower
    WHERE dist > 0 limit 1""" % name

    cursor.execute(query)

    row = cursor.fetchone()

    if row is  None:
        print "[ERROR]bounding not found for %s" % name
        end_bound ="none"

    else:
        end_bound = row["eco_name"]

    return start_bound, end_bound

def get_locus_name(name):

    db = MySQLdb.connect(host="localhost", user="amirbar", db="article_refactor_24_3_2016")
    cursor = db.cursor(MySQLdb.cursors.DictCursor)

    query="""SELECT name
        FROM ecocyc_genes
        WHERE blattner_id='%s'""" % name

    cursor.execute(query)

    row = cursor.fetchone()

    if row is None:
        if name == "orphan":
            return "orphan"

        print "[ERROR]No name matching this id '%s'" % name

        return "not_found"

    return row["name"]


def upload_thomason_cross_table():

    TABLE_NAME = "thomason_cross"
    result = []

    with open(THOMASON_CROSS_TABLE, "rb") as fl:

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


    db = MySQLdb.connect(host="localhost", user="amirbar", db="article_refactor_24_3_2016")
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

    cursor.execute("ALTER TABLE %s CHANGE antisense antisense int not null" % TABLE_NAME)
    cursor.execute("ALTER TABLE %s CHANGE _primary _primary int not null" % TABLE_NAME)
    cursor.execute("ALTER TABLE %s CHANGE _secondary _secondary int not null" % TABLE_NAME)
    cursor.execute("ALTER TABLE %s CHANGE internal internal int not null" % TABLE_NAME)
    cursor.execute("ALTER TABLE %s CHANGE putative_asrna putative_asrna int not null" % TABLE_NAME)
    cursor.execute("ALTER TABLE %s CHANGE tss_position tss_position int not null" % TABLE_NAME)
    cursor.execute("ALTER TABLE %s CHANGE min_start min_start int not null" % TABLE_NAME)
    cursor.execute("ALTER TABLE %s CHANGE max_end max_end int not null" % TABLE_NAME)
    cursor.execute("ALTER TABLE %s CHANGE is_inside is_inside int not null" % TABLE_NAME)


    db.commit()


def test_cross_tss(threshold):
    print "running"

    db = MySQLdb.connect(host="localhost", user="amirbar", db="article_refactor_24_3_2016")
    cursor = db.cursor(MySQLdb.cursors.DictCursor)

    query="""SELECT signif_chimeras.name, type, min_start, max_end, start_1, strand,
    (start_1>=min_start and start_1<=max_end) as inside, thomason_full.locus_tag as tss_locus, product,
    antisense, _primary, _secondary, internal, putative_asrna, product

    FROM signif_chimeras, thomason_full
    WHERE strand=strand_1 and (min_start - %(threshold)s)<=start_1 and start_1<=(max_end + %(threshold)s)
    ORDER BY name""" % {"threshold": threshold}

    cursor.execute(query)

    rows = []

    row = cursor.fetchone()

    while row is not None:
        rows.append(row)
        row = cursor.fetchone()

    header = ["name",
              "type",
              "strand",
              "min_start",
              "max_end",
              "is_inside",
              "tss_position",
              "_primary",
              "_secondary",
              "internal",
              "antisense",
              "putative_asrna",
              "distance_from_start",
              "tss_locus",
              "tss_locus_name",
              # "bound_start",
              # "bound_end",
              "product"]

    with open(THOMASON_CROSS_TABLE, "wb") as fl:
        fl.write("\t".join(header) + "\n")

        for row in rows:

            values = [row["name"],
                      row["type"],
                      row["strand"],
                      row["min_start"],
                      row["max_end"],
                      row["inside"],
                      row["start_1"],
                      row["_primary"],
                      row["_secondary"],
                      row["internal"],
                      row["antisense"],
                      row["putative_asrna"]]

            if row["strand"] == "+":
                values.append(int(row["start_1"]) - int(row["min_start"]))

            else:
                values.append(int(row["max_end"]) - int(row["start_1"]))

            # bound_start, bound_end = get_bounding_genes(row["name"])

            values.append(row["tss_locus"])
            values.append(get_locus_name(row["tss_locus"]))
            # values.append(bound_start)
            # values.append(bound_end)
            values.append(row["product"])

            fl.write("\t".join(["%s" % val for val in values]) + "\n")

test_cross_tss(100)
upload_thomason_cross_table()
