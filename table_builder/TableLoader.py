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

if (__name__ == "__main__"):
    TableLoader().load("parsed_files/s5.table")