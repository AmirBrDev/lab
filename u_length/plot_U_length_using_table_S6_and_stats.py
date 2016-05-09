#!/usr/bin/env python

"""
Use Table S6 of the paper with the entities, their clasification and U longest
stretch to plot the length of the U of the different clasifications and compare
to the annotated genes terminators.
"""
from EcoCyc.parse_genes_dat import *
import sys
import argparse
from Bio import SeqIO
from scipy.stats import mannwhitneyu
from numpy import mean
import matplotlib.pyplot as plt
from collections import defaultdict
import random

def process_command_line(argv):
    """
    Return a 2-tuple: (settings object, args list).
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object, replace the description
    parser = argparse.ArgumentParser(
        description='Plot poly-U lengths',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-g', '--genome', default='/home/users/assafp/deep-data/reference_genomes/E_coli/genome.fa',
        help='The genome sequence.')
    parser.add_argument(
        '-m', '--min_ints', type=int, default=5,
        help='Minimal interactions to include an entry.')
    parser.add_argument(
        '-n', '--num_of_genes', type=int, default=5,
        help='Minimal number of genes in entity to include in plot.')
    parser.add_argument(
        'table',
        help='Table S6 of the paper.')

        

    settings = parser.parse_args(argv)

    return settings

def cons_T(seq):
    maxs = 0
    pos = 0
    conss = 0
    while (pos < len(seq)):
        if seq[pos] == 'T':
            conss += 1
        else:
            maxs = max(conss, maxs)
            conss = 0
        pos += 1
    maxs = max(maxs, conss)
    return maxs

def read_S6_table(tfile, min_ints):
    """
    Read the table (tab delimited) and return a dict with the classification
    as key and a list with U length as values.
    Assume classification in third column and U=length in last (7)
    Number of interactions are in column 4
    Arguments:
    - `tfile`: The table file
    - `min_ints`: Minimal number of interactions to incldue
    """
    class_len = defaultdict(list)
    with open(tfile) as tin:
        for line in tin:
            spl = line.strip().split()
            if int(spl[3]) >= min_ints:
                class_len[spl[2]].append(int(spl[13]))
    return class_len




def main(argv=None):
    settings = process_command_line(argv)
    uid_pos, uid_names, uid_tudata, sRNAs_list, other_RNAs_list =\
        read_genes_data()
    terms, ttp = read_terminators_data(get_type=True)
    gen = SeqIO.read(settings.genome, 'fasta').seq
    # Read the mRNA terminators
    rest_ters = {}
    for sr in set(uid_pos.keys()) - set(sRNAs_list):    
        for tud in uid_tudata[sr]:
            if tud in ttp and ttp[tud]=="Rho-Independent-Terminators":
                rest_ters[tud] = gen[terms[tud][0]-25:terms[tud][1]+26]
                if uid_pos[sr][3] == '-':
                    rest_ters[tud] = rest_ters[tud].reverse_complement()
    mRNAs_tlens = [cons_T(s) for s in rest_ters.values()]
    clens = read_S6_table(settings.table, settings.min_ints)
    print "mean mRNA terminator poly-U length: %g"%mean(mRNAs_tlens)
    plt.figure()
    plt.scatter([random.random() for i in range(len(mRNAs_tlens))], mRNAs_tlens)
    plt.plot([0,1], [mean(mRNAs_tlens)]*2, 'r')
    spos = 3
    tick_labs = ['mRNAs']
    tick_poss = [0.5]

    ccat = 'CDS'
    clist = clens[ccat]
    lengths_cds = clist
    if len(clist) >= settings.num_of_genes:
        plt.scatter([random.random()+spos for i in range(len(clist))], clist)
        plt.plot([spos, spos + 1], [mean(clist)]*2, 'r')
        tick_labs.append(ccat)
        tick_poss.append(spos+0.5)
        spos += 3

    ccat = '5UTR'
    clist = clens[ccat]
    lengths_5utr = clist
    if len(clist) >= settings.num_of_genes:
        plt.scatter([random.random()+spos for i in range(len(clist))], clist)
        plt.plot([spos, spos + 1], [mean(clist)]*2, 'r')
        tick_labs.append(ccat)
        tick_poss.append(spos+0.5)
        spos += 3

    ccat = 'IGR'
    clist = clens[ccat]
    lengths_igr = clist
    if len(clist) >= settings.num_of_genes:
        plt.scatter([random.random()+spos for i in range(len(clist))], clist)
        plt.plot([spos, spos + 1], [mean(clist)]*2, 'r')
        tick_labs.append(ccat)
        tick_poss.append(spos+0.5)
        spos += 3

    ccat = '3UTR'
    clist = clens[ccat]
    lengths_3utr = clist
    if len(clist) >= settings.num_of_genes:
        plt.scatter([random.random()+spos for i in range(len(clist))], clist)
        plt.plot([spos, spos + 1], [mean(clist)]*2, 'r')
        tick_labs.append(ccat)
        tick_poss.append(spos+0.5)
        spos += 3

    ccat = 'sRNA'
    clist = clens[ccat]
    lengths_srna = clist
    if len(clist) >= settings.num_of_genes:
        plt.scatter([random.random()+spos for i in range(len(clist))], clist)
        plt.plot([spos, spos + 1], [mean(clist)]*2, 'r')
        tick_labs.append(ccat)
        tick_poss.append(spos+0.5)
        spos += 3

    max_len = max(len(mRNAs_tlens),
                  len(lengths_cds), 
                  len(lengths_5utr), 
                  len(lengths_igr),
                  len(lengths_3utr), 
                  len(lengths_srna))

    print "\t".join(["mRNAs", "CDS", "5UTR", "IGR", "3UTR", "sRNA"])

    for i in range(max_len):

        vals = []
        if i < len(mRNAs_tlens):
            vals.append(str(mRNAs_tlens[i]))
        else:
            vals.append("")

        if i < len(lengths_cds):
            vals.append(str(lengths_cds[i]))
        else:
            vals.append("")

        if i < len(lengths_5utr):
            vals.append(str(lengths_5utr[i]))
        else:
            vals.append("")

        if i < len(lengths_igr):
            vals.append(str(lengths_igr[i]))
        else:
            vals.append("")

        if i < len(lengths_3utr):
            vals.append(str(lengths_3utr[i]))
        else:
            vals.append("")

        if i < len(lengths_srna):
            vals.append(str(lengths_srna[i]))
        else:
            vals.append("")

        print "\t".join(vals)

    # for ccat, clist in clens.items():
    #     if len(clist) < settings.num_of_genes:
    #         continue
    #     print "mean %s poly-U length: %g"%(ccat, mean(clist))
    #     m, p = mannwhitneyu(clist, mRNAs_tlens)
    #     print "Mann-Whitney-U stat and p-value: %g %g"%(m, p)
    #     plt.scatter([random.random()+spos for i in range(len(clist))], clist)
    #     plt.plot([spos, spos + 1], [mean(clist)]*2, 'r')
    #     tick_labs.append(ccat)
    #     tick_poss.append(spos+0.5)
    #     spos += 3
    plt.xlim([-1, spos-1])
    plt.xticks(tick_poss, tick_labs)
    plt.ylabel('Longest U track')
    plt.savefig('longest_U_stretch_scatters.eps')
    plt.savefig('longest_U_stretch_scatters.tif')
    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
