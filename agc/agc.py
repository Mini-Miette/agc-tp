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
import time
import argparse
import sys
import os
import math
import gzip
import statistics
import numpy as np
from tqdm import tqdm
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
            msg = f"{path} is a directory"
        else:
            msg = f"{path} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__,
                                     usage=f"{sys.argv[0]} -h")
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
    global NB_SEQUENCES

    # Building the dictionary.
    seq_dict = {}
    fasta_reader = read_fasta(amplicon_file, minseqlen)
    cpt_fasta_reader = 0
    for seq in fasta_reader:
        cpt_fasta_reader += 1
        if seq in seq_dict:
            seq_dict[seq] += 1
        else:
            seq_dict[seq] = 1
    print('cpt_fasta_reader:', cpt_fasta_reader)

    NB_SEQUENCES = 0
    for key in sorted(seq_dict, key=seq_dict.get, reverse=True):
        if seq_dict[key] >= mincount:
            NB_SEQUENCES += 1
    print('NB_SEQUENCES:', NB_SEQUENCES)

    # Yielding (sequence, count) in decreasing count order.
    for key in sorted(seq_dict, key=seq_dict.get, reverse=True):
        if seq_dict[key] >= mincount:
            yield (key, seq_dict[key])
        else:
            break  # No need to go throught the rest since it's ordered.


def get_unique(ids):
    """Create a empty dictionnary with non-redondant ids as keys.
    """
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    """Retrieve a list of all elements included in both lists.
    """
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """Retrieve at least 4 chunks (fragments) of given size for input
    sequence.
    """
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError(f"Sequence length ({len_seq}) is too short to be "
                         f"splitted in 4 chunks of size {chunk_size}")
    return [sequence[i:i+chunk_size]
            for i in range(0, len_seq, chunk_size)
            if i+chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]


def get_identity(alignment_list):
    """
    Prend en entrée une liste de séquences alignées au format
    ["SE-QUENCE1", "SE-QUENCE2"] et retourne le pourcentage
    d'identité entre les deux.
    """
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """Update kmer_dict with sequences kmers.
    DO NOT record multiple occurences (id_seq is added once).

    Kmer_dict keys are kmers sequences while values storing a list
    of sequences indices (id_seq) if this kmer can be found in
    matching id_seq sequence.
    """
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer in kmer_dict:
            if id_seq in kmer_dict[kmer]:
                continue
            kmer_dict[kmer].append(id_seq)
        else:
            kmer_dict[kmer] = [id_seq]
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    """Identify 2 best matching reads for sequence.
    Returns a list of 2 id_seq corresponding to reads that match
    best with current studied sequence (here, a chunk).
    """
    all_found = []
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer in kmer_dict:
            all_found += kmer_dict[kmer]
    best_mates = Counter(all_found).most_common(2)
    #print('best_mates:', best_mates)

    return [seq_id for seq_id, count in best_mates]


def detect_chimera(perc_identity_matrix):
    """Detecting chimera sequences.
    Si l’écart type moyen des pourcentages d’identité des segments est supérieur à 5
    et que 2 segments minimum de notre séquence montrent une similarité différente
    à un des deux parents, nous identifierons cette séquence comme chimérique.
    """
    # Standard deviations mean is expected above 5.0.
    std_devs = [statistics.stdev(values) for values in perc_identity_matrix]
    ident_std_enough = statistics.mean(std_devs) > 5

    # Neither [0,0,0,0] nor [1,1,1,1] are wanted for below list,
    # not chimeral reads.
    fragment_max_pos = [0 if values[0] == max(values) else 1
                                       for values in perc_identity_matrix]
    diff_segments = not (sum(fragment_max_pos) == 0
                         or sum(fragment_max_pos) == 4)

    # Final answer requires both boolean to be True.
    return ident_std_enough and diff_segments


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

    non_chimerics = []  # Non chimeric sequences.
    # The first 2 sequences are declared non chimeric by default. So we jump
    # directly to the third one.
    non_chimerics.append(next(reader, None))
    non_chimerics.append(next(reader, None))
    #print(len(non_chimerics))
    # We are going to build a k-mer dictionary making an inventory of all
    # the non chimeric sequences k-mers.
    # The dictionary needs to be initialized with those first sequences,
    # then it will be updated each time a new non chimeric sequence is
    # identified.
    kmer_dict = get_unique_kmer({}, non_chimerics[0][0], 0, kmer_size)
    kmer_dict = get_unique_kmer(kmer_dict, non_chimerics[1][0], 1, kmer_size)

    for current_seq, current_count in tqdm(reader, unit='seq',
                                           total=NB_SEQUENCES):
        #print(len(non_chimerics))

        parents = []  # Parent sequences to the current sequence.

        # Splitting the sequence into non overlaping segments.
        current_segments = get_chunks(current_seq, chunk_size)

        # Searching for the nc sequence that best matches the segment k-mers.
        all_mates_seg = []
        for current_seg in current_segments:
            mates_seg = search_mates(kmer_dict, current_seg, kmer_size)
            all_mates_seg.append(mates_seg)
        #print(all_mates_seg)

        # Identifying the nc sequences that mate with all segments.
        common_mates = all_mates_seg[0]
        for mates in all_mates_seg[1:]:
            common_mates = common(common_mates, mates)
        #print(common_mates)

        # The segment share several k-mers with nc sequences. We have
        # found a parent sequence.
        for mates in common_mates:
            #print(mates)
            #print(non_chimerics[mates])
            parents.append(non_chimerics[mates])
        #print(len(parents))

        if len(parents) >= 2:
            # We use a np.array to hold all the identidy % between the
            # current sequence and each parent sequence.
            # The columns hold the sequences in the same order than the
            # parents list. The rows hold the segments.
            nb_segments = math.ceil(len(current_seq) / chunk_size)
            #print(nb_segments, len(parents))
            identities = np.zeros((nb_segments, len(parents)))

            for i in range(len(parents)):
                parent_segments = get_chunks(parents[i][0], chunk_size)
                aligns = [nw.global_align(current_seg, parent_seg,
                                          gap_open=-1, gap_extend=-1,
                                          matrix=matrix_file)
                          for current_seg, parent_seg
                          in zip(current_segments, parent_segments)]
                identities[:, i] = [get_identity(aligns)]

            # Updating the list of non chimeric sequences.
            is_chimeric = detect_chimera(identities)
            if not is_chimeric:
                non_chimerics.append((current_seq, current_count))
                # And updating the k-mer dictionary.
                kmer_dict = get_unique_kmer(kmer_dict, current_seq,
                                            len(non_chimerics)-1, kmer_size)

        for seq, count in non_chimerics:
            yield (seq, count)



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
    non_chimeric_reader = chimera_removal(amplicon_file, minseqlen, mincount,
                                          chunk_size, kmer_size)
    #reader = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    matrix_file = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                               "MATCH"))
    OTU_list = []

    # We initialise with the sequence with the most abundance.
    OTU_list.append(next(non_chimeric_reader, None))

    for (nc_seq, nc_count) in non_chimeric_reader:
        identities = []
        for otu in OTU_list:

            # First we align the sequences.
            align = nw.global_align(otu[0], nc_seq, gap_open=-1,
                                    gap_extend=-1,
                                    matrix=matrix_file)
            # Then we compute the identity.
            identities.append(get_identity(align))

        # If the sequence is less than 97% identical to all the others,
        # we have found a new OTU.
        # print(identities)
        if all([identity <= 97 for identity in identities]):
            #print('on ajoute')
            OTU_list.append(((nc_seq, nc_count)))

    # print(otu_list)
    return OTU_list


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
        amplicon_file = '/home/laura/Documents/M2BI/Omiques/agc-tp/data/amplicon.fasta.gz'
        #amplicon_file = '/home/laura/Documents/M2BI/Omiques/agc-tp/data/amplicon_head_10000.fasta.gz'
        #amplicon_file = '/home/laura/Documents/M2BI/Omiques/agc-tp/tests/test_sequences.fasta.gz'
        minseqlen = 400
        mincount = 10
        chunk_size = 100
        kmer_size = 8
        output_file = '/home/laura/Documents/M2BI/Omiques/agc-tp/output/OTU_amplicon.fasta'
        #output_file = '/home/laura/Documents/M2BI/Omiques/agc-tp/output/OTU_amplicon_head_10000.fasta'
        #output_file = '/home/laura/Documents/M2BI/Omiques/agc-tp/output/OTU_test_sequences.fasta'

    OTU_list = abundance_greedy_clustering(amplicon_file, minseqlen, mincount,
                                           chunk_size, kmer_size)
    print(f'Nb of OTUs found: {len(OTU_list)}')
    write_OTU(OTU_list, output_file)



if __name__ == '__main__':
    start_time = time.time()
    main()
    end_time = time.time() - start_time
    print(f'OTU search done in {end_time // 60:.0f} min'
          f' {end_time % 60:.2f} s.')
