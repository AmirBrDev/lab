from twisted.internet.test.test_gireactor import skip

__author__ = 'amirbar'

import unittest

class my_test_case(unittest.TestCase):

    def  setUp(self):
        x = 1

    def tearDown(self):
        x = 0


class some_test_module(my_test_case):

    def test_example_success(self):
        self.assertTrue(True)

    def test_example_failure(self):
        self.assertTrue(False)


if (__name__ == '__main__'):
    unittest.main()