__author__ = 'amirbar'

import unittest
from table_builder.TableLoader import TableLoader
from table_builder.Table import Table
import os
import copy

class table_test_case(unittest.TestCase):

    def  setUp(self):
        self.loader = TableLoader()

    def tearDown(self):
        self.loader = None

    def test_dump(self):

        rows = self.loader.load("test_files/test_structure_4")
        expected_rows = copy.deepcopy(rows)
        expected_rows[0][Table.UNIQUE_ID_FIELD] = "0"

        table = self.loader.createTable("test", rows)


        file_name = "test_output_files/table_dump"

        if os.path.exists(file_name):
            os.remove(file_name)

        table.dump(file_name)

        actual_rows = self.loader.load(file_name)

        self.assertSequenceEqual(actual_rows, expected_rows)

    def test_insert(self):

        table = Table("dummy")
        table.insert((1, 2, "+"), (11, 12, "+"))

        self.assertTrue(len(table) == 1)

        table.insert((21, 22, "-"), (31, 32, "-"))
        self.assertTrue(len(table) == 2)

        self.assertTrue(table._dctData.has_key("1;2;+;11;12;+"))
        self.assertTrue(table._dctData["1;2;+;11;12;+"] == {Table.UNIQUE_ID_FIELD : 0})
        self.assertTrue(table._dctData.has_key("21;22;-;31;32;-"))
        self.assertTrue(table._dctData["21;22;-;31;32;-"] == {Table.UNIQUE_ID_FIELD : 1})

    def test_overlap_unmatched(self):

        table = Table("dummy")
        table.insert((1, 2, "+"), (11, 12, "+"))
        table.insert((21, 22, "-"), (31, 32, "-"))
        table.insert((100, 101, "-"), (100, 101, "-"))


        result = table.is_overlaps(3, 4, "+")
        self.assertEqual(len(result), 1)
        self.assertFalse(result[0][0])
        self.assertEqual(result[0][2], 1)

        result = table.is_overlaps(3, 4, "-")
        self.assertEqual(len(result), 1)
        self.assertFalse(result[0][0])
        self.assertEqual(result[0][2], 17)

        result = table.is_overlaps(0, 0, "+")
        self.assertEqual(len(result), 1)
        self.assertFalse(result[0][0])
        self.assertEqual(result[0][2], 1)

        result = table.is_overlaps(13, 14, "+")
        self.assertEqual(len(result), 1)
        self.assertFalse(result[0][0])
        self.assertEqual(result[0][2], 1)

        result = table.is_overlaps(1, 2, "-")
        self.assertEqual(len(result), 1)
        self.assertFalse(result[0][0])
        self.assertEqual(result[0][2], 19)

    def test_overlap_match(self):
        table = Table("dummy")
        table.insert((1, 2, "+"), (11, 12, "+"))
        table.insert((21, 22, "-"), (31, 32, "-"))

        result = table.is_overlaps(1, 12, "+")
        self.assertEqual(len(result), 1)
        self.assertTrue(result[0][0])
        self.assertEqual(result[0][2], 0)

        table.insert((12, 13, "+"), (14, 15, "+"))
        result = table.is_overlaps(1, 13, "+")
        self.assertEqual(len(result), 2)
        self.assertTrue(result[0][0])
        self.assertEqual(result[0][2], 0)
        self.assertTrue(result[1][0])
        self.assertEqual(result[1][2], -11)

if (__name__ == '__main__'):
    unittest.main()