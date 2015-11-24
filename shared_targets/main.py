__author__ = 'amirbar'

import MySQLdb


def get_partial_gene(name):
    pass

class DBHandler(object):
    """
    Supplies useful methods to access db data without logging every time
    """

    def __init__(self, host, user, db, table_name_list):
        self.db = MySQLdb.connect(host=host, user=user, db=db)
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

    def get_tf_names(self):
        query = """SELECT tf_name
                FROM regulon_tf_gene"""

        self.cursor.execute(query)

        row = self.cursor.fetchone()

        names = []

        while row is not None:
            if row["tf_name"] not in names:
                names.append(row["tf_name"])

            row = self.cursor.fetchone()

        return names

    def get_ppi_genes_names(self):
         query = """SELECT gene_1
                FROM ppi_interactions_as_genes
                UNION
                SELECT gene_2
                FROM ppi_interactions_as_genes"""

         self.cursor.execute(query)

         row = self.cursor.fetchone()

         names = []

         while row is not None:
            if row["gene_1"] not in names:
                names.append(row["gene_1"])

            row = self.cursor.fetchone()

         return names

    def find_targets(self, names):
        """
        Extracts all the interacting RNAs for each name in the names list given
        :param names: a list of names to find their targets
        :return: dictionary where the keys are the given names and the values are lists of targets names
        """

        names_to_targets = {}

        for name in names:

            union_statement = ""

            for table_name in self.table_name_list[1:]:
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
            %(union_statement)s""" % {"table_name": self.table_name_list[0],
                                      "name": name,
                                      "union_statement": union_statement}

            self.cursor.execute(query)

            row = self.cursor.fetchone()

            names_to_targets[name] = []

            while row is not None:
                if row["rna1_name"] not in names_to_targets[name]:
                    names_to_targets[name].append(row["rna1_name"])

                row = self.cursor.fetchone()

        return names_to_targets

    def find_our_targets(self, name):

        union_statement = ""

        for table_name in self.table_name_list[1:]:
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
        %(union_statement)s""" % {"table_name": self.table_name_list[0],
                                  "name": name,
                                  "union_statement": union_statement}

        self.cursor.execute(query)

        row = self.cursor.fetchone()

        targets = []

        while row is not None:
            if row["rna1_name"] not in targets:
                targets.append(row["rna1_name"])

            row = self.cursor.fetchone()

        return targets

    def find_tf_targets(self, tf_name):
        """
        Finds the targets names of a given tf
        :param tf_name: the tf name
        :return: a list containing names of target genes
        """
        query = """SELECT gene_name
                FROM regulon_tf_gene
                WHERE tf_name = '%s'""" % tf_name

        self.cursor.execute(query)

        row = self.cursor.fetchone()

        targets = []

        while row is not None:
            if row["gene_name"] not in targets:
                targets.append(row["gene_name"])

            row = self.cursor.fetchone()

        return targets

    def find_ppi_genes_targets(self, ppi_gene_name):
        query = """SELECT gene_1
                FROM ppi_interactions_as_genes
                WHERE gene_2 = '%(ppi_gene_name)s'
                UNION
                SELECT gene_2
                FROM ppi_interactions_as_genes
                WHERE gene_1 = '%(ppi_gene_name)s'""" % {"ppi_gene_name": ppi_gene_name}

        self.cursor.execute(query)

        row = self.cursor.fetchone()

        targets = []

        while row is not None:
            if row["gene_1"] not in targets:
                targets.append(row["gene_1"])

            row = self.cursor.fetchone()

        return targets

class Runner(object):
    """
    Runs the program
    """

    def __init__(self):
        our_tables = ["signif_chimeras_of_iron_limitation_cl",
                      "signif_chimeras_of_log_phase_cl",
                      "signif_chimeras_of_stationary_cl"]

        self.db = DBHandler(host="localhost", user="amirbar", db="amir", table_name_list=our_tables)

        # Generate target table for our genes
        self.total_names, ecocyc_ids, types = self.db.get_names_from_tables()
        self.total_tf_names = self.db.get_tf_names()
        self.total_ppi_names = self.db.get_ppi_genes_names()

        self.targets = self.db.find_targets(self.total_names)
        self.filtered_targets = {}

        for name, target_list in self.targets.items():

            filtered_names = self.filter_name(name)

            if len(self.filter_name(name)) > 0:

                for filtered_name in filtered_names:

                    if self.filtered_targets.has_key(name):
                        self.filtered_targets[filtered_name].union(self.filter_targets(target_list))

                    else:
                        self.filtered_targets[filtered_name] = self.filter_targets(target_list)

    def filter_genes(self, names):
        result = []

        for name in names:
            result += self.filter_name(name)

        return set(result)

    def filter_targets(self, targets_names):
        result = []

        for name in targets_names:
            result += self.filter_name(name, True)

        return set(result)

    def filter_name(self, name, keep_cds=False):
        NAME_DELIMITER = "."
        result = []

        if "3utr" in name:
            result.append(name.split(NAME_DELIMITER)[0])

        elif "tu" in name:
            result += name.split(NAME_DELIMITER)[:2]

        elif "igr" in name:
            result += name.split(NAME_DELIMITER)[:2]

        elif "5utr" in name and keep_cds:
            result.append(name.split(NAME_DELIMITER)[0])

        elif keep_cds and len(name.split(NAME_DELIMITER)) == 1:
            result += [name]

        return result

    def run(self,
            candidate_name_list,
            candidate_target_retrieve_func,
            interaction_name_list,
            interaction_target_retrieve_func):
        """
        Prints all the pairs of regulatory elements which has common target.

        :param candidate_name_list: the group we want to check from our data
        :param candidate_target_retrieve_func: function to retrieve the targets of our candidate
        :param interaction_name_list: the group we want to compare with
        :param interaction_target_retrieve_func: function to retrieve the targets of our compared group
        :return: 0 on success anything else on failure.
        """       

        # filtered_gene_name_set = self.filter_genes(candidate_name_list)
        # print filtered_gene_name_set

        for our_name in candidate_name_list:
            for interaction_name in interaction_name_list:

                if interaction_name in self.filter_name(our_name, True):
                    interaction_targets = interaction_target_retrieve_func(interaction_name)
                    candidate_targets = self.filter_targets(candidate_target_retrieve_func(our_name))

                    # self.filtered_targets[our_name]
                    intersection = set(interaction_targets).intersection(candidate_targets)

                    if len(intersection) > 0:
                        print our_name, interaction_name, intersection

        return 0


def our_targets_finder(name):

    mapping = {"1.3utr": ["2.5.igr", "3.3utr"],
               "1": ["1", "2"],
               "2.4.tu": ["5","6", "7"]}

    return mapping[name]

def their_targets_finder(name):
    mapping = {"1": ["3", "20", "5", "2", "1"],
               "2": ["22", "6"]}

    return mapping[name]

if __name__ == "__main__":
    print "running"

    runner = Runner()

    candidate_name_list=["1.3utr", "2.4.tu"]
    their_names = ["2", "1"]

    # exit(runner.run(runner.total_names,
    #                 runner.db.find_our_targets,
    #                 runner.total_ppi_names,
    #                 runner.db.find_ppi_genes_targets))
    # exit(runner.run(runner.total_names,
    #                 runner.db.find_our_targets,
    #                 runner.total_tf_names,
    #                 runner.db.find_tf_targets))

    exit(runner.run(candidate_name_list,
                    our_targets_finder,
                    their_names,
                    their_targets_finder))