__author__ = 'amirbar'

from DBCrosser import cross_tables, format_our_single

our_tables = ["signif_chimeras_of_iron_limitation_cl",
              "signif_chimeras_of_log_phase_cl",
              "signif_chimeras_of_stationary_cl"]

raghavan_tables = ["raghavan_s5",
                   "raghavan_s6",
                   "raghavan_s7"]

bilusic_tables = ["bilusic_s1",
                  "bilusic_s2",
                  "bilusic_s3",
                  "bilusic_s4"]

zhang_tables = ["zhang_s3_2013_sheet2008",
                "zhang_s3_2013_sheet2009",
                "zhang_s4_2013_sheet2008",
                "zhang_s4_2013_sheet2009"]

lybecker_tables = ["lybecker_s1",
                   "lybecker_s2"]

tss_tables = ["tss"]

def run_set(name, tables):

    print "#" * 100

    print "running set '%s'\n" % name


    # Go over our tables
    for current_our_table in our_tables:

        print "*" * 100
        print "processing our table: %s\n" % current_our_table


        # for each of them go over their tables in the current set
        for current_their_table in tables:

            print "*" * 50
            print "processing their table: %s" % current_their_table

            # generate crossing table and unique table by our name
            cross_tables(current_our_table,
                         current_their_table,
                         format_our_single,
                         50,
                         "result/%s/%s_to_%s.table" % (name, current_our_table, current_their_table))

run_set("raghavan", raghavan_tables)
run_set("bilusic", bilusic_tables)
run_set("zhang", zhang_tables)
run_set("lybecker", lybecker_tables)
run_set("tss", tss_tables)