__author__ = 'amirbar'

class ExtractedConverter(object):
    def __init__(self):
        pass

    def convertS5(self, path):

        # Parse into dictionary
        new_format_list = []
        fl = open(path, "rb")

        for line in fl.readlines():

            lst = eval(line)

            dct = {'left_end' : int(lst[0]),
                   'right_end' : int(lst[1]),
                   'length' : lst[2],
                   'sRNA-name' : lst[3],
                   'MEV' : lst[4],
                   'PRM' : lst[5]
                   }

            new_format_list.append(dct)

        fl.close()

        self.writeToFile("./s5.table", new_format_list)

    def convertS6(self, path):

        # Parse into dictionary
        new_format_list = []
        fl = open(path, "rb")

        for line in fl.readlines():

            lst = eval(line)

            dct = {'left_end' : int(lst[0]),
                   'right_end' : int(lst[1]),
                   'orientation' : lst[2]
                   }

            new_format_list.append(dct)

        fl.close()

        self.writeToFile("./s6.table", new_format_list)

    def convertS7(self, path):

        # Parse into dictionary
        new_format_list = []
        fl = open(path, "rb")

        for line in fl.readlines():

            lst = eval(line)

            dct = { 'name' : lst[0],
                    'left_end' : int(lst[1]),
                    'right_end' : int(lst[2]),
                    'length' : int(lst[3]),
                    'strand' : lst[4],
                    'MEV' : lst[5],
                    'PRM' : lst[6],
                    'source' : lst[7]
                   }

            new_format_list.append(dct)

        fl.close()

        self.writeToFile("./s7.table", new_format_list)

    def writeToFile(self, path, table_as_list):
        fl = open(path, "wb")

        header = "\t".join(key for key in table_as_list[0].keys())
        fl.write("%s\n" % header)

        for dct in table_as_list:
            row = "\t".join(str(val) for val in dct.values())
            fl.write("%s\n" % row)

        fl.close()

if (__name__ == "__main__"):

    converter = ExtractedConverter()

    converter.convertS5("./unparsed_files/s5_parsed.txt")
    converter.convertS6("./unparsed_files/s6_parsed.txt")
    converter.convertS7("./unparsed_files/s7_parsed.txt")