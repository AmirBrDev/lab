__author__ = 'amirbar'

import unittest

from table_builder.TableLoader import TableLoader
from table_builder.BuilderExceptions.StructureException import StructureException
from table_builder.Globals import TableGlobals
from table_builder.Table import Table

class loader_test_case(unittest.TestCase):

    def  setUp(self):
        self.loader = TableLoader()

    def tearDown(self):
        self.loader = None

    def test_simple_load(self):

        result = self.loader.load("test_files/test_structure_1")

        expected_dict = {"field1" : "val1",
                         "field2" : "val2"}

        expected_keys = expected_dict.keys()
        expected_values = expected_dict.values()

        actual_keys = result[0].keys()
        actual_values = result[0].values()

        self.assertSequenceEqual(actual_keys, expected_keys)
        self.assertSequenceEqual(actual_values, expected_values)

    def test_bad_structure(self):

        self.assertRaises(StructureException, lambda: self.loader.load("test_files/test_structure_2"))
        self.assertRaises(StructureException, lambda: self.loader.load("test_files/test_structure_3"))

    def test_table_create(self):

        row_list = self.loader.load("test_files/test_structure_4")
        table = self.loader.createTable("test", row_list)

        expected_keys = "10;19;+;20;29;+"
        expected_values = {"info" : "more info",
                           Table.UNIQUE_ID_FIELD : 0}

        self.assertTrue(len(table), 1)
        self.assertTrue(table._dctData.has_key(expected_keys))
        self.assertTrue(expected_values == table._dctData[expected_keys])



if (__name__ == '__main__'):
    unittest.main()