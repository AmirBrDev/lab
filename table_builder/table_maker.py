__author__ = 'amirbar'

from TableLoader import GeneTableLoader, PDFTableLoader, ZhangTableLoader, BilusicTableLoader, TssMasterTableLoader
from TableLoader import LybeckerTableLoader, LybeckerS2TableLoader
from Globals import TableGlobals
from Table import Table


def makeFinalZhangTables():

    tables = ["Table-s3-zhang-2013-sheet2008",
              "Table-s3-zhang-2013-sheet2009",
              "Table-s4-zhang-2013-sheet2008",
              "Table-s4-zhang-2013-sheet2009"]

    loader = ZhangTableLoader()

    for table in tables:
        result = loader.loadUnprocessed("./zhang/csv/%s.csv" % table)
        loader.createTable("dummy", result).dump("./zhang/final/%s.table" % table)

def makeFinalBilusicTables():

    tables = ["2014RNABIOL0069R_TableS1",
              "2014RNABIOL0069R_TableS2",
              "2014RNABIOL0069R_TableS3",
              "2014RNABIOL0069R_TableS4"]

    loader = BilusicTableLoader()

    for table in tables:
        result = loader.loadUnprocessed("./bilusic/csv/%s.csv" % table)
        loader.createTable("dummy", result).dump("./bilusic/final/%s.table" % table)

def makeFinalTssTables():
    tables = ["JB.02096-14_zjb999093409sd1-3"]

    loader = TssMasterTableLoader()

    for table in tables:
        result = loader.loadUnprocessed("./tss/csv/%s.csv" % table)
        loader.createTable("dummy", result, True).dump("./tss/final/%s.table" % table)

def makeFinalLybeckerTables():
    tables = ["sd01"]

    loader = LybeckerTableLoader()

    for table in tables:
        result = loader.loadUnprocessed("./lybecker/csv/%s.csv" % table)
        loader.createTable("dummy", result).dump("./lybecker/final/%s.table" % table)

    loader = LybeckerS2TableLoader()

    result = loader.loadUnprocessed("./lybecker/csv/lybecker_s2.csv")
    loader.createTable("dummy", result).dump("./lybecker/final/lybecker_s2.table")


if (__name__ == "__main__"):
    # printMatchAnalysis()
    #addDirectionalityToPdfTables()
    #makeFinalZhangTables()
    # makeFinalBilusicTables()
    # makeFinalTssTables()
    makeFinalLybeckerTables()