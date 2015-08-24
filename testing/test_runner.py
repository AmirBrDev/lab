__author__ = 'amirbar'

import unittest
import test_suite1
import test_suite2

if __name__ == '__main__':
    test_classes_to_run = [test_suite1.some_test_module,
                           test_suite2.another_test_case]

    loader = unittest.TestLoader()

    suites_list = []

    for test_class in test_classes_to_run:
        suite = loader.loadTestsFromTestCase(test_class)
        suites_list.append(suite)

    big_suite = unittest.TestSuite(suites_list)

    runner = unittest.TextTestRunner()
    results = runner.run(big_suite)