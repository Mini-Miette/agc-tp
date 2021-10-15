#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

from collections import Counter
import argparse
import sys
import os
import numpy as np
import math
import gzip
import statistics
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Bénédicte Noblet et Laura Xénard"
__copyright__ = "Universite de Paris"
__credits__ = ["Bénédicte Noblet et Laura Xénard"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Bénédicte Noblet et Laura Xénard"
__email__ = "laura.xenard@protonmail.com"
__status__ = "Developpement"


parse = False


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default=400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default=10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default=100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default=8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    """Read gz fasta file and get sequences one at a time.

    Parameters
    ----------
    fastq_file: string
        Path to gz.fastq file, identifier starts with ">" while
        sequences may be writtent on one or several lines.

    Returns
    -------
    iterator
        An iterator operating on reads
    """
    # Version Bénédicte mais ne prend pas la dernière séquence.
    with gzip.open(amplicon_file, 'rt') as my_file:

        activeone = False
        sequence = ''
        for line in my_file:
            if str(line).startswith(">"):
                if activeone and len(sequence) >= minseqlen:
                    yield sequence
                else:
                    activeone = True
                line = next(my_file, None)
                sequence = str(line).strip()
            else:
                sequence += str(line).strip()
        if activeone and len(sequence) >= minseqlen:
            yield sequence

# ===================================================================
#         for line in my_file:
#             sequence = ''
#             while line[0] != '>':
#                 sequence += line
#                 line = next(my_file, None)
#             if len(sequence) >= minseqlen:
#                 print(sequence)
#                 yield sequence
# ===================================================================

# ===================================================================
#         for line in my_file:
#             if line[0] != '>' and len(line) >= minseqlen:
#                 yield line
#
# ===================================================================

# ===================================================================
# # Version Bénédicte mais trop de séquences (59)
#         activeone = False
#         sequence = ''
#         for line in my_file:
#             if str(line).startswith(">"):
#                 if activeone:
#                     activeone = False
#                 else:
#                     line = next(my_file, None)
#                     sequence = str(line).strip()
#                     activeone = True
#             else:
#                 if activeone:
#                     sequence += str(line).strip()
#                 else:
#                     if len(sequence) >= minseqlen:
#                         yield sequence
# ===================================================================


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """Dereplicate an amplicon file.

    Args:
        amplicon_file (str): Path of the amplicon file.
        minseqlen (int): Minimum sequence length accepted.
        mincount (int): Minimun sequence occurrency accepted.

    Yields:
        key (str): Sequence.
        int: Number of occurrences of the sequence.

    """
    # Building the dictionary.
    seq_dict = {}
    fasta_reader = read_fasta(amplicon_file, minseqlen)
    for seq in fasta_reader:
        if seq in seq_dict:
            seq_dict[seq] += 1
        else:
            seq_dict[seq] = 1
    # print(seq_dict)

    # Yielding (sequence, count) in decreasing count order.
    for key in sorted(seq_dict, key=seq_dict.get, reverse=True):
        if seq_dict[key] >= mincount:
            yield (key, seq_dict[key])
        else:
            break  # No need to go throught the rest it's ordered.


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """Retrieve chunks (fragments) of given size for input sequence.
    """
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError("Sequence length ({}) is too short to be splitted in 4"
                         " chunk of size {}".format(len_seq, chunk_size))
    return [sequence[i:i+chunk_size]
            for i in range(0, len_seq, chunk_size)
            if i+chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]


def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    pass


def search_mates(kmer_dict, sequence, kmer_size):
    pass


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Yield a generator of non chimeric sequences.

    Args:
        amplicon_file (str): Path of the amplicon file.
        minseqlen (int): Minimum sequence length accepted.
        mincount (int): Minimun sequence occurrency accepted.
        chunk_size (int): Size of the sequence chunks.
        kmer_size (int): Size of the k-mers to be considered.

    Yields:
        None.

    """
    matrix_file = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                               "MATCH"))
    reader = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    non_chimeras = []  # Non chimeric sequences.
    # The first 2 sequences are declared non chimeric by default. So we jump
    # directly to the third one.
    non_chimeras.append(next(reader, None))
    non_chimeras.append(next(reader, None))
    # We are going to build a k-mer dictionary making an inventory of all
    # the non chimeric sequences k-mers.
    # The dictionary needs to be initialized with those first sequences,
    # then it will be updated each time a new non chimeric sequence is
    # identified.
    kmer_dict = get_unique_kmer({}, non_chimeras[0], 0, kmer_size)
    kmer_dict = get_unique_kmer(kmer_dict, non_chimeras[1], 1, kmer_size)

    for current_seq, current_count in reader:

        parents = []  # Parent sequences to the current sequence.

        # Splitting the sequence into non overlaping segments.
        current_segments = get_chunks(current_seq, chunk_size)

        # We want to compare each segment to the corresponding segment of
        # the nc (non chimeric) sequences.
        for nc_seq, nc_count in non_chimeras:
            nc_segments = get_chunks(nc_seq, chunk_size)
            for current_seg, nc_seg in zip(current_segments, nc_segments):

                # We want to compare the k-mers of each segment so it's
                # easier to first build a Counter to count the number
                # of occurrencies of each k-mer in a segment.
                current_seg_kmers_cnt = Counter()
                kmer_reader = cut_kmer(current_seg, kmer_size)
                for kmer in kmer_reader:
                    current_seg_kmers_cnt[kmer] += 1
                # We do the same for the non chimeric segment.
                nc_seg_kmers_cnt = Counter()
                kmer_reader = cut_kmer(nc_seg, kmer_size)
                for kmer in kmer_reader:
                    nc_seg_kmers_cnt[kmer] += 1

                # How many identical k-mers do we have between the 2 segments?
                # We look at the intersection between the 2 Counters.
                shared_kmers = [key for key in current_seg_kmers_cnt
                                if key in nc_seg_kmers_cnt]

                if len(shared_kmers) > 1:
                    # The segments share several k-mers. We have found a
                    # parent sequence.
                    parents.append((nc_seq, nc_count))
                    break  # No need to check the other segments.
                elif (len(shared_kmers) == 1
                      and current_seg_kmers_cnt[shared_kmers[0]] > 1
                      and nc_seg_kmers_cnt[shared_kmers[0]] > 1):
                    # There's only one shared k-mer and it appears at least
                    # twice in both segments. We have found a parent sequence.
                    parents.append((nc_seq, nc_count))
                    break  # No need to check the other segments.

        if len(parents) >= 2:
            # We use a np.array to hold all the identidy % between the
            # current sequence and each parent sequence.
            # The rows hold the sequences in the same order than the parents
            # list. The columns are for the segments.
            nb_segments = math.ceil(len(current_seq) / chunk_size)
            identities = np.zeros((len(parents), nb_segments))

            for i in range(len(parents)):
                parent_segments = get_chunks(parents[i], chunk_size)
                aligns = [nw.global_align(current_seg, parent_seg,
                                          gap_open=-1, gap_extend=-1,
                                          matrix=matrix_file)
                          for current_seg, parent_seg
                          in zip(current_segments, parent_segments)]
                identities[i, :] = [get_identity(aligns)]

                print(identities)

        # TODO: partie 4.
        # Vérifier si le cas où il n'y pas le même nombre de segments
        # apparaît.




    # Must yield (sequence, count)


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount,
                                chunk_size, kmer_size):
    """Perform an abundance greedy clustering on a sequences file.

    Args:
        amplicon_file (str): Path to gz.fastq file.
        minseqlen (int): Minimum sequence length accepted.
        mincount (int): Minimun sequence occurrency accepted.
        chunk_size (int): Size of the chunks to be considered.
        kmer_size (int): Size of the k-mers to be considered.

    Returns:
        otu_list (list(str, int)): List of OTUs and their respective number
                                   of occurences.

    """
    reader = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    matrix_file = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                               "MATCH"))
    otu_list = []

    # We initialise with the sequence with the most abundance.
    otu_list.append(next(reader, None))

    for element in reader:
        # print(element)
        identities = []
        for otu in otu_list:

            # First we align the sequences.
            align = nw.global_align(otu[0], element[0], gap_open=-1,
                                    gap_extend=-1,
                                    matrix=matrix_file)
            # Then we compute the identity.
            identities.append(get_identity(align))

        # If the sequence is less than 97% identical to all the others,
        # we have found a new OTU.
        # print(identities)
        if all([identity <= 97 for identity in identities]):
            #print('on ajoute')
            otu_list.append(element)

    # print(otu_list)
    return otu_list


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_OTU(OTU_list, output_file):
    """Write OTUs in a file following fastas format.

    Args:
        OTU_list (list(str)): lists of identified OTU.
        output_file (str): path of out file.

    Yields:
        None.
    """
    with open(output_file, 'wt') as my_out_file:
        for i, sequence in enumerate(OTU_list):
            my_out_file.write(f">OTU_{i+1} occurrence:{sequence[1]}\n")
            my_out_file.write(fill(sequence[0])+"\n")

# ==============================================================
# Main program
# ==============================================================


def main():
    """
    Main program function
    """
    if parse:
        # Get arguments
        args = get_arguments()
        amplicon_file = args.amplicon_file
        minseqlen = args.minseqlen
        mincount = args.mincount
        chunk_size = args.chunk_size
        kmer_size = args.kmer_size
        output_file = args.output_file
        # Votre programme ici

    else:
        amplicon_file = '/home/laura/Documents/M2BI/Omiques/agc-tp/tests/test_sequences.fasta.gz'
        minseqlen = 400
        mincount = 10
        chunk_size = 100
        kmer_size = 8
        output_file = '/home/laura/Documents/M2BI/Omiques/agc-tp/output/test_output.txt'

# =============================================================================
#     reader = read_fasta(amplicon_file, 1)
#     cpt = 0
#     for seq in reader:
#         # print(seq)
#         cpt += 1
#         # print('--')
#     print(cpt)
# =============================================================================

    #abundance_greedy_clustering(amplicon_file, 200, 1, chunk_size, kmer_size)
    chimera_removal(amplicon_file, 100, 1, 10, 4)

# =============================================================================
#     dereplication_reader = dereplication_fulllength(amplicon_file, 200, 3)
#     derep_1 = next(dereplication_reader)
#     derep_2 = next(dereplication_reader)
#     # Should be the most abundant sequence: seq4 counted 5 times
#     assert(derep_1[0] == "ACTACGGGGCGCAGCAGTAGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGTATGAAGAAGGTTTTCGGATCGTAAAGTACTGTTGTTAGAGAAGAACAAGGATAAGAGTAACTGCTTGTCCCTTGACGGTATCTAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAGTTAGTGGGCGTAAAGCGCGCGCAGGCGGTCTTTTAAGTCTGATGTCAAAGCCCCCGGCTTAACCGGGGAGGGTCATTGGAAACTGGAAGACTGGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGATATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAA")
#     assert(derep_1[1] == 5)
#     # Should be the second most abundant sequence: seq3 counted 4 times
#     assert(derep_2[0] == "TAGGGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTCTTCGGATCGTAAAACTCTGTTATTAGGGAAGAACATATGTGTAAGTAACTGTGCACATCTTGACGGTACCTAATCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTACAGCGCG")
#     assert(derep_2[1] == 4)
#     try:
#         derep_3 = next(dereplication_reader)
#         # derep_3 should be empty
#         assert(len(derep_3) == 0)
#     except StopIteration:
#         # Congrats only two sequences to detect
#         assert(True)
# =============================================================================


if __name__ == '__main__':
    main()
