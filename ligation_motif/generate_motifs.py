__author__ = 'amirbar'

import os
from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet import Gapped, IUPAC
import csv

GAP = "-"
ALPHABET = Gapped(IUPAC.unambiguous_dna)

def format_input(chimera_frags_list):

    result = []

    for chimera in chimera_frags_list:
        first_frag, second_frag = chimera
        result.append((Seq(first_frag, ALPHABET), Seq(second_frag, ALPHABET)))

    return  result


def format_file_input(path_to_file):

    FRAG_1_SEQ_INDEX = 6
    FRAG_2_SEQ_INDEX = 8

    fragment_as_seq_list = []

    with open(path_to_file, "rb") as csvfile:
        reader = csv.reader(csvfile, delimiter='\t', quotechar='|')

        for row in reader:
            fragment_as_seq_list((Seq(FRAG_1_SEQ_INDEX, ALPHABET), Seq(row[FRAG_2_SEQ_INDEX], ALPHABET)))


def get_max_len_in_chimera(fragment_as_seq_list):

    first_frag_list, second_frag_list = zip(*fragment_as_seq_list)

    return len(max(first_frag_list, key=len)), len(max(second_frag_list, key=len))

def pad_sequences(chimera_frags_list):

    max_first, max_second = get_max_len_in_chimera(chimera_frags_list)

    result = []

    for chimera in chimera_frags_list:
        first_frag, second_frag = chimera

        first_frag = (max_first - len(first_frag)) * GAP + first_frag
        second_frag += (max_second - len(second_frag)) * GAP

        result.append((first_frag, second_frag))

    return result

def merge_frgaments(padded_fragments_list):

    result = []

    for chimera in padded_fragments_list:
        first_frag, second_frag = chimera
        result.append(Seq(str(first_frag + second_frag), ALPHABET))

    return result

def create_sequnce_file(filename, sequence_list):

    with open(filename, "wb") as fl:

        for sequence in sequence_list:
            fl.write(">\n")
            fl.write(str(sequence) + "\n")


def run(chimera_frags_list):

    # fragment_as_seq_list = format_input(chimera_frags_list)
    fragment_as_seq_list = format_file_input("reads_fusion.txt")

    padded_fragments_list = pad_sequences(fragment_as_seq_list)

    merged_fragments_list = merge_frgaments(padded_fragments_list)

    chimera_motifs = motifs.create(merged_fragments_list, ALPHABET)

    print chimera_motifs
    print chimera_motifs.counts

    create_sequnce_file("logo_data", merged_fragments_list)

    os.system("~/.local/bin/weblogo --format PNG < logo_data  > logo.png")







sequences = [("TACAAA", "TACGC"),
             ("TACAC", "TACCC"),
             ("TACCC", "TACCC"),
             ("TACCC", "AATGC"),
             ("AATGC", "AATGCA")]

run(sequences)