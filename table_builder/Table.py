__author__ = 'amirbar'

from Globals import TableGlobals

class Table(object):
    """
    This class represents a table defined by us or other sources.
    The table contains different kind of data sorted according to unique
    key - the loci start, end and directionality
    """

    class TableIterator(object):
        """
        This class is used to iterate over the table rows
        """

        def __init__(self, table):
            self._rows = table._dctData
            self._iterator = self._rows.iteritems()

        def next(self):
            return self._iterator.next()

        def __iter__(self):
            return self._iterator

    TABLE_DELIMITER = "\t"
    ID_DELIMITER = ";"
    UNIQUE_ID_FIELD = "unique_id"
    NAME_FIELD = "name"
    POS_STRAND = "+"
    NEG_STRAND = "-"

    def __init__(self, name):
        """
        Initialize a table with a name
        :param name: the name for the table
        :return: new table
        """

        self._name = name
        self._dctData = {}

    def __iter__(self):
        return Table.TableIterator(self)

    def __len__(self):
        return len(self._dctData)

    def get_name(self):
        return self._name

    def dump(self, path):
        """
        Writes the table to a file
        :param path: the path to write to
        :return: None
        """

        if (len(self._dctData) == 0):
            return

        id_keys = [TableGlobals.FIRST_START_BASE_KEY,
                   TableGlobals.FIRST_END_BASE_KEY,
                   TableGlobals.FIRST_STRAND_KEY,
                   TableGlobals.SECOND_START_BASE_KEY,
                   TableGlobals.SECOND_END_BASE_KEY,
                   TableGlobals.SECOND_STRAND_KEY]

        for key in id_keys:
            self._dctData.values()[0][key] = None

        header = Table.TABLE_DELIMITER.join(key for key in self._dctData.values()[0])

        fl = open(path, "wb")

        fl.write("%s\n" % header)

        table_as_list = []

        for key, value in self._dctData.items():

            id_values = key.split(Table.ID_DELIMITER)

            for id_key, id_val in zip(id_keys, id_values):

                value[id_key] = id_val

            table_as_list.append(value)

        for dct in table_as_list:
            row = "\t".join(str(val) for val in dct.values())

            fl.write("%s\n" % row)

        fl.close()

    def insert(self, first, second, append_dup = False, **kwargs):
        """
        Insert a new entry to the table
        :param first: a tuple (start, end, strand)
        :param second: a tuple (start, end, strand)
        :param kwargs: a dictionary containing additional info
        :return: None
        """
        additional_data = {}
        additional_data.update(kwargs)
        additional_data[Table.UNIQUE_ID_FIELD] = len(self._dctData)

        # reorder the start and end if they were given in wrong order (sometimes on the negative strand
        # it is written opposite)
        start, end, strand = first

        if (int(start) > int(end)):
            first = (end, start, strand)

        start, end, strand = second

        if (int(start) > int(end)):
            second = (end, start, strand)

        # Generate unique id
        id = ("%s;%s;%s;" % first) + ("%s;%s;%s" % second)

        if (self._dctData.has_key(id)):

            if (not append_dup):
                print "[Table:insert]warning multiple entries has the unique id: %s" % id

            else:

                for key, value in additional_data.items():

                    if key != Table.UNIQUE_ID_FIELD:
                        self._dctData[id][key] += ";%s" % value

        if not append_dup or not self._dctData.has_key(id):
            self._dctData[id] = additional_data

    def is_overlaps(self, start, end, strand_val):
        """
        Checks whether a given segment overlaps with
        any of the entries in the table
        :param start: the start position
        :param end: the end position
        :param strand_val: the strand string (e.g "positive" or "negative")
        :return: a tuple containing True or False whether it has overlapping segment
                 and the id of the closest segment found (None if there isn't one)
        """
        closest_segment = None
        closest_segment_distance = None
        distance = None

        matches = []

        for key, value in self._dctData.items():

            first_entry_start, first_entry_end, first_entry_strand,\
                second_entry_start, second_entry_end, second_entry_strand = key.split(Table.ID_DELIMITER)

            first_entry_start = int(first_entry_start)
            first_entry_end = int(first_entry_end)
            second_entry_start = int(second_entry_start)
            second_entry_end = int(second_entry_end)

            #print entry_start, entry_end, entry_strand

            strand_values = []

            if (strand_val == "none"):
                strand_values.append(Table.POS_STRAND)
                strand_values.append(Table.NEG_STRAND)
            else:
                strand_values.append(strand_val)


            # in case we don't sure of the strand check for both
            for strand in strand_values:

                # We only care about matches on the same strand
                if (strand != first_entry_strand and strand != second_entry_strand):
                    continue

                if ((start >= first_entry_start and start <= first_entry_end) or
                    (end >= first_entry_start and end <= first_entry_end) or
                    (first_entry_start >= start and first_entry_start <= end) or
                    (first_entry_end >= start and first_entry_end <= end)) and strand == first_entry_strand:

                    matches.append((True, (key, value), start - first_entry_start))

                elif ((start >= second_entry_start and start <= second_entry_end) or
                    (end >= second_entry_start and end <= second_entry_end) or
                    (second_entry_start >= start and second_entry_start <= end) or
                    (second_entry_end >= start and second_entry_end <= end)) and strand == second_entry_strand:

                    matches.append((True, (key, value), start - second_entry_start))

                first_distance = min(abs(start - first_entry_start),
                                     abs(end - first_entry_end),
                                     abs(start - first_entry_end),
                                     abs(end - first_entry_start))

                second_distance = min(abs(start - second_entry_start),
                                      abs(end - second_entry_end),
                                      abs(start - second_entry_end),
                                      abs(end - second_entry_start))

                distance = min(first_distance, second_distance)

                if (closest_segment == None):
                    closest_segment = (key, value)
                    closest_segment_distance = distance

                elif (distance < closest_segment_distance):
                    closest_segment = (key, value)
                    closest_segment_distance = distance

        if ((len(matches) == 0) and (distance is not None)):

            matches.append((False, closest_segment, closest_segment_distance))

        return matches

    def findByField(self, field_name, search_value):
        for key, value in self._dctData.items():
            if (value[field_name] == search_value):
                return (key, value)

        return (None, None)

    def findById(self, id):
        return self.findByField(Table.UNIQUE_ID_FIELD, id)

    def findByName(self, name):
        return self.findByField(Table.NAME_FIELD, name)


class GeneTable(Table):

    OTHER_NAMES_FIELD = "other_names"

    def findByOtherNames(self, name):
        for key, value in self._dctData.items():
            if (name in value[GeneTable.OTHER_NAMES_FIELD]):
                return (key, value)

        return (None, None)

if __name__ == "__main__":

    table = Table("Known sRNA")

    dct = {"name" : "test1",
           "name2" : "test2"}

    table.insert((10, 12, "+"), (13, 100, "+"), **dct)
    table.insert((6, 9, "+"), (300, 400, "-"), **dct)

    print table._dctData

    print table.is_overlaps(1, 7, "-")

    print "-" * 100
    print "iterator checks"
    print "-" * 100

    for row in table:
        print row