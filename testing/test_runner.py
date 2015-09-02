__author__ = 'amirbar'

import unittest
import test_table
import test_table_loader

if __name__ == '__main__':
    test_classes_to_run = [test_table.table_test_case,
                           test_table_loader.loader_test_case]

    loader = unittest.TestLoader()

    suites_list = []

    for test_class in test_classes_to_run:
        suite = loader.loadTestsFromTestCase(test_class)
        suites_list.append(suite)

    big_suite = unittest.TestSuite(suites_list)

    runner = unittest.TextTestRunner()
    results = runner.run(big_suite)