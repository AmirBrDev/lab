__author__ = 'amirbar'

class Table(object):
    """
    This class represents a table defined by us or other sources.
    The table contains different kind of data sorted according to unique
    key - the loci start, end and directionality
    """

    UNIQUE_ID_FIELD = "unique_id"
    NAME_FIELD = "NAME"
    OTHER_NAMES_FIELD = "other_names"
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

    def insert(self, first, second, **kwargs):
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

        id = ("%s;%s;%s;" % first) + ("%s;%s;%s" % second)

        if (self._dctData.has_key(id)):
            print "[Table:insert]warning multiple entries has the unique id: %s" % id

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

        matches = []

        for key, value in self._dctData.items():

            first_entry_start, first_entry_end, first_entry_strand,\
                second_entry_start, second_entry_end, second_entry_strand = key.split(";")

            first_entry_start = int(first_entry_start)
            first_entry_end = int(first_entry_end)
            second_entry_start = int(second_entry_start)
            second_entry_end = int(second_entry_end)

            #print entry_start, entry_end, entry_strand

            strand_values = []

            if (strand_val == None):
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
                    (first_entry_end >= start and first_entry_end <= end)):

                    matches.append((True, value[Table.UNIQUE_ID_FIELD]))

                elif ((start >= second_entry_start and start <= second_entry_end) or
                    (end >= second_entry_start and end <= second_entry_end) or
                    (second_entry_start >= start and second_entry_start <= end) or
                    (second_entry_end >= start and second_entry_end <= end)):

                    matches.append((True, value[Table.UNIQUE_ID_FIELD]))

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
                    closest_segment = value[Table.UNIQUE_ID_FIELD]
                    closest_segment_distance = distance

                elif (distance < closest_segment_distance):
                    closest_segment = value[Table.UNIQUE_ID_FIELD]
                    closest_segment_distance = distance

        if (len(matches) == 0):

            matches.append((False, closest_segment))

        return matches

    def findByField(self, field_name, search_value):
        for key, value in self._dctData.items():

            compared = value[field_name]

            if (type(search_value) is str):
                search_value = search_value.lower()
                compared = value[field_name].lower()

            if (compared == search_value):
                return (key, value)

        return (None, None)

    def findById(self, id):
        return self.findByField(Table.UNIQUE_ID_FIELD, id)

    def findByName(self, name):
        return self.findByField(Table.NAME_FIELD, name)

    def findByOtherNames(self, name):
        for key, value in self._dctData.items():
            if (name in value[Table.OTHER_NAMES_FIELD]):
                return (key, value)

        return (None, None)

if __name__ == "__main__":

    table = Table("Known sRNA")

    dct = {"name" : "test1",
           "name2" : "test2"}

    table.insert((10, 12, "+"), (13, 100, "+"), **dct)
    table.insert((6, 9, "+"), (6, 9, "+"), **dct)

    print table._dctData

    print table.is_overlaps(1, 7, "+")