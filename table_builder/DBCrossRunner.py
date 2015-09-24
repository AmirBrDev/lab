__author__ = 'amirbar'

from DBCrosser import cross_tables, format_our_single, format_their_our

our_tables = ["signif_chimeras_of_iron_limitation_cl",
              "signif_chimeras_of_log_phase_cl",
              "signif_chimeras_of_stationary_cl"]

raghavan_tables = [("raghavan_s5", "mev >= 1"),
                   ("raghavan_s6", "1=1"),
                   ("raghavan_s7", "mev >= 1"),
                   ("raghavan_s8", "1=1"),
                   ("raghavan_2", "mev >= 1")]

bilusic_tables = [("bilusic_s1", "1=1"),
                  ("bilusic_s2", "1=1"),
                  ("bilusic_s3_1", "1=1"),
                  ("bilusic_s3_2", "1=1"),
                  ("bilusic_s4_1", "1=1"),
                  ("bilusic_s4_2", "1=1")]

zhang_tables = [("zhang_s3_2013_sheet2008", "1=1"),
                ("zhang_s3_2013_sheet2009", "1=1"),
                ("zhang_s4_2013_sheet2008", "1=1"),
                ("zhang_s4_2013_sheet2009", "1=1")]

lybecker_tables = [("lybecker_s1", "1=1"),
                   ("lybecker_s2", "1=1")]

tss_tables = [("thomason", "1=1")]

def run_set(name, tables, format_func):

    print "#" * 100

    print "running set '%s'\n" % name


    # Go over our tables
    for current_our_table in our_tables:

        print "*" * 100
        print "processing our table: %s\n" % current_our_table


        # for each of them go over their tables in the current set
        for current_their_table in tables:

            print "*" * 50
            print "processing their table: %s" % current_their_table[0]

            # generate crossing table and unique table by our name
            cross_tables(current_our_table,
                         current_their_table,
                         format_func,
                         0,
                         "result/%s/%s_to_%s.0.table" % (name, current_our_table, current_their_table[0]))

format_func = format_our_single
# format_func = format_their_our

run_set("raghavan", raghavan_tables, format_func)
# run_set("bilusic", bilusic_tables, format_func)
# run_set("zhang", zhang_tables, format_func)
# run_set("lybecker", lybecker_tables, format_func)
# run_set("tss", tss_tables, format_func)