__author__ = 'amirbar'


class StructureException(BaseException):
    def __init__(self, msg):

        self.message = "[Structure Exception] %s" % msg