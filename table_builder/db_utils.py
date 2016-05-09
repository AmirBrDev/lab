import MySQLdb

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


    db.commit()
