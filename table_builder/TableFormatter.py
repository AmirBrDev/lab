__author__ = 'amirbar'

from TableLoader import TableLoader


class TableFormatter(object):

    def __init__(self):
        pass

    @staticmethod
    def _add_names_to_data(names_to_data, row):

        # if new entry
        if not row[0] in names_to_data.keys():
            names_to_data[row[0]] = row[1:]

        # otherwise save the closest one
        else:

            curr_overlaps = names_to_data[row[0]][2] == "true"
            candidate_overlaps = row[3] == "true"

            try:
                curr_distance = abs(int(names_to_data[row[0]][5]))
            except ValueError:
                curr_distance = 999999999

            try:
                candidate_distance = abs(int(row[6]))
            except ValueError:
                candidate_distance = 999999999

            if curr_overlaps and candidate_overlaps and candidate_distance < curr_distance:
                names_to_data[row[0]] = row[1:]

            elif not curr_overlaps and candidate_overlaps:
                names_to_data[row[0]] = row[1:]

            elif (not curr_overlaps) and \
                 (not candidate_overlaps) and \
                 (candidate_distance < curr_distance):
                names_to_data[row[0]] = row[1:]

    @staticmethod
    def _get_key(row):

        score = 0

        # give score according to found
        if row[3] == "false":
            score += 100

        # give score according to type
        if row[2] == "3utr":
            score += 1

        elif row[2] == "5utr":
            score += 2

        elif row[2] == "as":
            score += 3

        elif row[2] == "as_with_cis_t":
            score += 4

        elif row[2] == "cis_as_with_trans_t":
            score += 5

        elif row[2] == "other-ncrna":
            score += 6

        elif row[2] == "srna":
            score += 7

        elif row[2] == "trna":
            score += 8

        elif row[2] == "tu":
            score += 9

        elif row[2] == "mrna":
            score += 10

        elif row[2] == "igr":
            score += 11

        else:
            raise BaseException("Unsupported type '%s'" % row[2])

        return score

    @staticmethod
    def _sort_names_to_data(names_to_data):

        to_sort = []

        for key, value in names_to_data.items():
            to_add = [key]
            to_add.extend(value)

            to_sort.append(to_add)

        return sorted(to_sort, key=TableFormatter._get_key)

    @staticmethod
    def split_reverse_rows(loader_type, file_path):

        loader = loader_type()

        table_rows = loader.load(file_path)

        names_to_data = {}

        header = ["RNA name",
                  "EcoCyc ID",
                  "type",
                  "overlaps",
                  "article name",
                  "article strand",
                  "distance",
                  "interactions",
                  "odds ratio",
                  "article table"]

        print "finished reading '%s'" % file_path

        for row in table_rows:
            # match_name = row["match to"].lower()
            # print row
            match_name = row["match to"].lower()

            # Handle first
            to_print = [row["rna1 name"],
                        row["rna1 ecocyc id"],
                        row["first_type"]]

            if match_name == "rna1":
                to_print.extend([row["overlaps"],
                                 row["article name"],
                                 row["article strand"],
                                 row["distance"],  # should be distance
                                 row["article table"]])

            else:
                to_print.extend([""] * 5)

            TableFormatter._add_names_to_data(names_to_data, to_print)

            # fl.write("%s\n" % "\t".join(to_print))

            # Handle second
            to_print = [row["rna2 name"],
                        row["rna2 ecocyc id"],
                        row["second_type"]]

            if match_name == "rna2":
                to_print.extend([row["overlaps"],
                                 row["article name"],
                                 row["article strand"],
                                 row["distance"],  # should be distance
                                 row["article table"]])

            else:
                to_print.extend([""] * 5)

            TableFormatter._add_names_to_data(names_to_data, to_print)

            # fl.write("%s\n" % )

        fl = open("%s.split" % file_path, "wb")

        fl.write("%s\n" % "\t".join(header))

        for curr in TableFormatter._sort_names_to_data(names_to_data):
            fl.write("%s\n" % "\t".join(curr))

        fl.close()

if __name__ == "__main__":

    TableFormatter.split_reverse_rows(TableLoader, "results/raghavan/reverse/assign-type-to-all-chimeras-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_all_interactions.with-type.table")

    # header = ["RNA name",
    #           "EcoCyc ID",
    #           "type",
    #           "overlaps",
    #           "article name",
    #           "article strand",
    #           "distance",
    #           "interactions",
    #           "odds ratio",
    #           "article table"]
    #
    # rows = [["aaa", "", "as", "false", "b1", "-", "100", "s5"],
    #         ["aaa", "", "as", "false", "b2", "-", "8", "s6"],
    #         ["aaa", "", "as", "false", "b3", "-", "12", "s7"],
    #         ["ccc", "", "tu", "false", "b1", "-", "100", "s5"],
    #         ["aaa", "", "as", "true", "b4", "-", "-13", "s7"],
    #         ["aaa", "", "as", "false", "b5", "-", "7", "s5"],
    #         ["ddd", "", "trna", "false", "b1", "-", "100", "s5"]]
    #
    # dct = {}
    #
    # for row in rows:
    #     TableFormatter._add_names_to_data(dct, row)
    #
    # print dct
    #
    # for row in TableFormatter._sort_names_to_data(dct):
    #     print row
