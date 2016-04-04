__author__ = 'amirbar'

import MySQLdb
from Table import Table

MCDOWELL_CROSS_TABLE = "mcdowell_cross.csv"

def generate_mcdowell_cross_(threshold):
    print "running"

    db = MySQLdb.connect(host="localhost", user="amirbar", db="article_refactor_24_3_2016")
    cursor = db.cursor(MySQLdb.cursors.DictCursor)

    query="""SELECT signif_chimeras.name, type, min_start, max_end, start_1, strand,
    (start_1>=min_start and start_1<=max_end) as inside

    FROM signif_chimeras, mcdowell

    WHERE strand=strand_1 and (min_start - %(threshold)s)<=start_1 and start_1<=(max_end + %(threshold)s) and
    m_geq_5='yes' and m_geq_3_4='yes'

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
			  "mcdowell_pos",
              "distance_from_start"]

    with open(MCDOWELL_CROSS_TABLE, "wb") as fl:
        fl.write("\t".join(header) + "\n")

        for row in rows:

            values = [row["name"],
                      row["type"],
                      row["strand"],
                      row["min_start"],
                      row["max_end"],
                      row["inside"],
                      row["start_1"]]

            if row["strand"] == "+":
                values.append(int(row["start_1"]) - int(row["min_start"]))

            else:
                values.append(int(row["max_end"]) - int(row["start_1"]))


            fl.write("\t".join(["%s" % val for val in values]) + "\n")

def upload_mcdowell_cross_table():

    TABLE_NAME = "mcdowell_cross"
    result = []

    with open(MCDOWELL_CROSS_TABLE, "rb") as fl:

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

    cursor.execute("ALTER TABLE %s CHANGE min_start min_start int not null" % TABLE_NAME)
    cursor.execute("ALTER TABLE %s CHANGE max_end max_end int not null" % TABLE_NAME)
    cursor.execute("ALTER TABLE %s CHANGE is_inside is_inside int not null" % TABLE_NAME)


    db.commit()


generate_mcdowell_cross_(100)
upload_mcdowell_cross_table()