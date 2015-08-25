__author__ = 'amirbar'

from Table import Table

class TableLoader(object):

    START_BASE_KEY = "from"
    END_BASE_KEY = "to"
    STRAND_KEY = "strand"

    def load(self, path):
        raise NotImplementedError

    def createTable(self, name, row_list):
        result = Table(name)

        for row in row_list:

            result.insert(row.pop(TableLoader.START_BASE_KEY),
                          row.pop(TableLoader.END_BASE_KEY),
                          row.pop(TableLoader.STRAND_KEY),
                          **row)

        return result

class PDFTableLoader(TableLoader):

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

class GeneTableLoader(TableLoader):
    def load(self, path):

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

            entry[TableLoader.START_BASE_KEY] = entry.pop('START-BASE')
            entry[TableLoader.END_BASE_KEY] = entry.pop('END-BASE')


        return result

class OurTableLoader(TableLoader):
    def load(self, path):
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


    loader = GeneTableLoader()
    result = loader.load("genes.col")
    #result = TableLoader().load_our_tables("./our_files/assign-type-to-all-chimeras-of-Iron_limitation_CL_FLAG207_208_305_all_fragments_l25.txt_all_interactions.with-type")


    table = loader.createTable("geneDB", result)


    res, id = table.is_overlaps(510860, 510866, "F")
    print res, table.findById(id)
