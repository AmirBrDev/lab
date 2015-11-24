__author__ = 'amirbar'

import MySQLdb

FIELD_DELIMITER = "\t"
GENE_NAME_DELIMITER = ","


def handle_tf_operon_file():

    table_name = "regulon_tf_gene"
    header = ["tf_name", "gene_name", "regulatory_effect", "operon_name", "operon_size"]

    rows = []

    with open("network_tf_operon.txt", "rb") as fl:

        for row in fl.readlines():

            row = row.lower()

            # Row structure is:
            # [TF_NAME]-del-[OPERON_NAME][gene_name_list delimited by ,]-del-[effect]-del-....
            row_values = row.split(FIELD_DELIMITER)

            # get gene names and operon
            gene_list_index = row_values[1].find("[")
            operon_name = row_values[1][:gene_list_index]
            genes = [val.strip() for val in row_values[1][gene_list_index+1:-1].split(GENE_NAME_DELIMITER)]

            for gene_name in genes:
                rows.append([row_values[0], gene_name, row_values[2], operon_name, len(genes)])

    for row in rows:
        print row

    db=MySQLdb.connect(host="localhost",user="amirbar",db="amir")
    cur = db.cursor(MySQLdb.cursors.DictCursor)

    fields = ", ".join("%s VARCHAR(200)" % key for key in header)
    cur.execute("CREATE TABLE %s (%s)" % (table_name, fields))

    for row in rows:
        values = ",".join("%s" % db.literal(str(val)) for val in row )
        cur.execute("INSERT INTO %s VALUES (%s)" % (table_name, values))

    db.commit()


def handle_ppi_file():
    table_name = "ppi_interactions_as_genes"
    header = ["gene_1", "gene_2"]

    rows = []

    with open("PPI_network_from_Rajagopala2014_including_literature.sif", "rb") as fl:

        for row in fl.readlines():

            row = row.lower().replace("\n", "")

            # Row structure is:
            # [first name] [second name]
            rows.append(row.split(FIELD_DELIMITER))

    for row in rows:
        print row

    db=MySQLdb.connect(host="localhost",user="amirbar",db="amir")
    cur = db.cursor(MySQLdb.cursors.DictCursor)

    fields = ", ".join("%s VARCHAR(200)" % key for key in header)
    cur.execute("CREATE TABLE %s (%s)" % (table_name, fields))

    for row in rows:
        values = ",".join("%s" % db.literal(str(val)) for val in row )
        cur.execute("INSERT INTO %s VALUES (%s)" % (table_name, values))

    db.commit()

handle_ppi_file()
# handle_tf_operon_file()