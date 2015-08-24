__author__ = 'amirbar'

class Table(object):
    """
    This class represents a table defined by us or other sources.
    The table contains different kind of data sorted according to unique
    key - the loci start, end and directionality
    """

    UNIQUE_ID_FIELD = "unique_id"

    def __init__(self, name):
        """
        Initialize a table with a name
        :param name: the name for the table
        :return: new table
        """

        self._name = name
        self._dctData = {}

    def insert(self, start, end, strand, **kwargs):
        """
        Insert a new entry to the table
        :param start: the start position
        :param end: the end position
        :param strand: the strand string (e.g "positive" or "negative")
        :param kwargs: a dictionary containing additional info
        :return: None
        """
        additional_data = {}
        additional_data.update(kwargs)
        additional_data[Table.UNIQUE_ID_FIELD] = len(self._dctData)

        self._dctData["%s;%s;%s" % (str(start), str(end), strand)] = additional_data

    def is_overlaps(self, start, end, strand):
        """
        Checks whether a given segment overlaps with
        any of the entries in the table
        :param start: the start position
        :param end: the end position
        :param strand: the strand string (e.g "positive" or "negative")
        :return: a tuple containing True or False whether it has overlapping segment
                 and the id of the closest segment found (None if there isn't one)
        """
        closest_segment = None
        closest_segment_distance = None

        for key, value in self._dctData.items():
            entry_start, entry_end, entry_strand = key.split(";")
            entry_start = int(entry_start)
            entry_end = int(entry_end)

            print entry_start, entry_end, entry_strand

            # We only care about matches on the same strand
            if (strand != entry_strand):
                continue

            if ((start >= entry_start and start <= entry_end) or
                (end >= entry_start and end <= entry_end) or
                (entry_start >= start and entry_start <= end) or
                (entry_end >= start and entry_end <= end)):
                return (True, value[Table.UNIQUE_ID_FIELD])

            distance = min(abs(start - entry_start),
                           abs(end - entry_end),
                           abs(start - entry_end),
                           abs(end - entry_start))

            if (closest_segment == None):
                closest_segment = value[Table.UNIQUE_ID_FIELD]
                closest_segment_distance = distance

            elif (distance < closest_segment_distance):
                closest_segment = value[Table.UNIQUE_ID_FIELD]
                closest_segment_distance = distance


        return (False, closest_segment)

if __name__ == "__main__":

    print min(3, 2)
    table = Table("Known sRNA")

    dct = {"name" : "test1",
           "name2" : "test2"}

    table.insert(10, 100, "positive", **dct)
    table.insert(6, 9, "positive", **dct)

    print table._dctData

    print table.is_overlaps(1, 7, "positive")