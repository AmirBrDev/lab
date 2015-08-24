__author__ = 'amirbar'

class TableLoader(object):

    def __init__(self):
        pass

    def load(self, path):
        result = []

        fl = open(path, "rb")

        name_list = fl.readline().split("\t")

        for line in fl.readlines():

            value_list = line.split("\t")

            dct = {}

            for key, val in zip(name_list, value_list):
                dct[key] = val

            result.append(dct)

        fl.close()

        return result

    def loadGeneTable(self, path):

        ignore = ["REPLICON",
                  "SWISS-PROT-ID",
                  "GENE-CLASS\n",
                  "PRODUCT-NAME",
                  "GENE-CLASS"]

        result = []

        fl = open(path, "rb")

        name_list = fl.readline().split("\t")

        for line in fl.readlines():

            value_list = line.split("\t")

            dct = {"other_names" : []}

            for key, val in zip(name_list, value_list):

                if (key == "SYNONYMS"):
                    if (val != ""):
                        dct["other_names"].append(val)
                elif (key in ignore):
                    continue
                else:
                    dct[key] = val

            result.append(dct)

        fl.close()

        for entry in result:
            if (int(entry['START-BASE']) < int(entry['END-BASE'])):
                entry["strand"] = "F"
            else:
                entry["strand"] = "R"


        return result


    def load_our_tables(self, path):
        ignore = ["RNA1 description",
                  "RNA2 description",
                  "odds ratio",
                  "total other interactions",
                  "interactions",
                  "other interactions of RNA1",
                  "other interactinos of RNA2",
                  "Fisher's exact test p-value"]

        result = []

        fl = open(path, "rb")

        dct = {}

        name_list = fl.readline().split("\t")

        for line in fl.readlines():

            value_list = line.split("\t")

            for key, val in zip(name_list, value_list):
                if (key in ignore):
                    continue
                else:
                    dct[key] = val

            result.append(dct)

        fl.close()


        return result

if (__name__ == "__main__"):


    #result = TableLoader().loadGeneTable("genes.col")
    result = TableLoader().load_our_tables("./our_files/assign-type-to-all-chimeras-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_all_interactions.with-type")

    i = 0
    for entry in result:
        print entry

        if (i == 10):
            break

        i += 1


