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
import gzip
import statistics
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
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file',
                        type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen',
                        type=int, default=400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount',
                        type=int, default=10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size',
                        type=int, default=100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size',
                        type=int, default=8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file',
                        type=str, default="OTU.fasta", help="Output file")
    parser.add_argument('-t', '-identity_treshold', dest='identity_treshold',
                        type=int, default=97,
                        help="percentage of identity for clusterization (default 97)")
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
    for seq in fasta_reader:
        if seq in seq_dict:
            seq_dict[seq] += 1
        else:
            seq_dict[seq] = 1
    if not seq_dict:
        sys.exit("No sequence meets minimum size requirement")

    # Getting the number of 'good' sequences for progress bar.
    NB_SEQUENCES = 0
    for key in sorted(seq_dict, key=seq_dict.get, reverse=True):
        if seq_dict[key] >= mincount:
            NB_SEQUENCES += 1

    # Yielding (sequence, count) in decreasing count order.
    for key in sorted(seq_dict, key=seq_dict.get, reverse=True):
        if seq_dict[key] >= mincount:
            yield (key, seq_dict[key])
        else:
            break  # No need to go throught the rest it's ordered.


def get_unique(ids):
    """Create a empty dictionnary with non-redondant ids as keys.
    """
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    """Retrieve a list of all elements included in both lists.
    """
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """Retrieve chunks (fragments) of given size for input sequence.
    """
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError(f"Sequence length ({len_seq}) is too short "
                         f"to be splitted in 4 chunks of size {chunk_size}")
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
    allfound = []
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer not in kmer_dict:
            continue
        allfound += kmer_dict[kmer]
    best_mates = Counter(allfound).most_common(2)
    return [seq_id for seq_id, count in best_mates]


def detect_chimera(perc_identity_matrix):
    """Detecting chimera sequences.
    Si l’écart type moyen des pourcentages d’identité des segments est
    supérieur à 5 et que 2 segments minimum de notre séquence montrent une
    similarité différente à un des deux parents, nous identifierons cette
    séquence comme chimérique.
    """
    # Standard deviations mean is expected above 5.0.
    std_devs = [statistics.stdev(values) for values in perc_identity_matrix]
    ident_std_enough = statistics.mean(std_devs) > 5

    # Either [0,0,0,0] nor [1,1,1,1] are not wanted for below list,
    # not chimeral reads.
    fragment_max_pos = [0 if values[0] == max(values) else 1
                                       for values in perc_identity_matrix]
    diff_segments = not (sum(fragment_max_pos) == 0
                         or sum(fragment_max_pos) == 4)

    # Final answer requires both boolean to be True.
    return ident_std_enough and diff_segments


def align_them(sequence1, sequence2, gap_open=-1, gap_extend=-1):
    """Handle nw.global_align easily.
    Besides, avoid repetitions of matrix and options.
    """
    matrix_file = os.path.abspath(os.path.join(
                                  os.path.dirname(__file__), "MATCH"))
    return nw.global_align(sequence1, sequence2,
                           gap_open=gap_open, gap_extend=gap_extend,
                           matrix=matrix_file)


def get_my_parents(current_segments, kmer_dict, kmer_size):
    """Embedding a part of chimera_removal to fit pylint requirments.
    Aim to identify parent sequences of current one using chunks.

    Args:
        current_segments(list(str)): fragments to study.
        chunk_size (int): Size of the sequence chunks.
        kmer_size (int): Size of the k-mers to be considered.
    Return:
        List: list of 2 parents sequences rank. Empty list if less.
    """
    closest = []
    parents = []

    # Compare to already identified non chimeric sequences.
    for current_seg in current_segments:
        best_mates = search_mates(kmer_dict, current_seg, kmer_size)
        closest.append(best_mates)
    # Perform some cleaning in parent list to keep only 2.
    for i in range(len(closest)-1):
        for j in range(i, len(closest)):
            candidate_parent = common(closest[i], closest[j])
            parents += candidate_parent
    unique_parents = set(parents)
    if len(unique_parents) <= 2:
        # Let's retrieve them.
        return list(unique_parents)
    # Else, sort them by decreasing order of importance.
    parent_counts = {}
    for parent in unique_parents:
        parent_counts[parent] = parents.count(parent)
    parents_counts = sorted(parent_counts,
                            key=parent_counts.get, reverse=True)
    parents = list(parents_counts[:2])
    return parents


def compute_my_identities(current_segments, parents, non_chimerics, chunk_size):
    """Embedding a part of chimera_removal to fit pylint requirments.
    Compute all identities for 2 parent sequences.

    Args:
        current_segments(list(str)): Chunk list for a sequence.
        parents(list(int)): List of 2 identified best mates.
        non_chimerics(list(tuple(str, int))):
              Active list of non chimeral sequences.
        chunk_size(int): Size of the sequence chunks.
    """
    identities = [[0.0,0.0] for current_seg in current_segments]
    for rank, parent in enumerate(parents):
        nc_segments = get_chunks(non_chimerics[parent][0], chunk_size)

        pairs = zip(range(len(current_segments)),
                          current_segments, nc_segments)
        for number, current_seg, nc_seg in pairs:
            aligns = align_them(current_seg, nc_seg)
            identities[number][rank] = get_identity(aligns)
    return identities


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Yield a generator of non chimeric sequences.

    Args:
        amplicon_file (str): Path of the amplicon file.
        minseqlen (int): Minimum sequence length accepted.
        mincount (int): Minimun sequence occurrency accepted.
        chunk_size (int): Size of the sequence chunks.
        kmer_size (int): Size of the k-mers to be considered.

    Yields:
        non_chimera (string, int): Sequence and number of reads.

    """
    reader = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    non_chimerics = []  # Non chimeric sequences (nc), the ones we want.

    # Create a k-mer dictionary to store nc sequences matching k-mers.
    kmer_dict = {}

    # The first 2 sequences are declared non chimeric by default.
    for i in range(2):
        seq_tuple = next(reader, None)
        if seq_tuple is None:
            sys.exit("No sequence meets minimum count requirement")
        yield seq_tuple
        non_chimerics.append(seq_tuple)
        kmer_dict = get_unique_kmer(kmer_dict,
                                    non_chimerics[i][0], i, kmer_size)

    # Operate on following sequence from reader.
    for current_seq, current_count in tqdm(reader, unit='seq',
                                           total=(NB_SEQUENCES-2)):
        if current_seq is None:
            # No more sequence meets minimum count requirement.
            break
        current_segments = get_chunks(current_seq, chunk_size)
        parents = get_my_parents(current_segments, kmer_dict, kmer_size)

        # Less than 2 parents, not a chimera...
        if len(parents) < 2:
            continue

        identities = compute_my_identities(current_segments, parents,
                                           non_chimerics, chunk_size)
        is_chimeric = detect_chimera(identities)

        # Updating the list and kmer dict of non chimeric sequences.
        if not is_chimeric:
            seqtuple = (current_seq, current_count)
            yield seqtuple
            non_chimerics.append(seqtuple)
            kmer_dict = get_unique_kmer(kmer_dict, non_chimerics[-1][0],
                                        len(non_chimerics)-1, kmer_size)


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount,
                                chunk_size, kmer_size, idt_treshold=97):
    """Perform an abundance greedy clustering on a sequences file.

    Args:
        amplicon_file (str): Path to gz.fastq file.
        minseqlen (int): Minimum sequence length accepted.
        mincount (int): Minimun sequence occurrency accepted.
        chunk_size (int): Size of the chunks to be considered.
        kmer_size (int): Size of the k-mers to be considered.
        idt_treshold (int): Identity percentage treshold.

    Returns:
        otu_list (list(str, int)): List of OTUs and their respective number
                                   of occurences.

    """
    reader = chimera_removal(amplicon_file, minseqlen, mincount, chunk_size,
                             kmer_size)
    otu_list = []

    # We initialise with the most abundant sequence.
    otu_list.append(next(reader, None))

    for element in reader:
        identities = []
        for otu in otu_list:
            # Compute alignment identity.
            align = align_them(otu[0], element[0])
            identities.append(get_identity(align))

        # If the sequence is less than 97% identical to all the others,
        # we have found a new OTU.
        if all([identity <= idt_treshold for identity in identities]):
            otu_list.append(element)

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

    # Get arguments
    args = get_arguments()

    # Reads and sizes parameters
    amplicon_file = args.amplicon_file
    minseqlen = args.minseqlen
    mincount = args.mincount
    chunk_size = args.chunk_size
    kmer_size = args.kmer_size
    idt_treshold = args.identity_treshold

    # Output filename for OTU list
    output_file = args.output_file

    # Abundance Greedy Clustering on filtered non chimerics reads
    OTU_list = abundance_greedy_clustering(amplicon_file, minseqlen, mincount,
                                           chunk_size, kmer_size, idt_treshold)
    print(f'Nb of OTUs found: {len(OTU_list)}')
    write_OTU(OTU_list, output_file)


if __name__ == '__main__':

    start_time = time.time()
    main()
    end_time = time.time() - start_time
    print(f'OTU search done in {end_time // 60:.0f} min'
          f' {end_time % 60:.2f} s.')
