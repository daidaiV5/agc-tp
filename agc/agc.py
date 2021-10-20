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

import argparse
import sys
import os
import gzip
import statistics
import operator
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "DAI_EID"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["DAI_EID"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "DAI_EID"
__status__ = "Developpement"


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
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    with gzip.open(amplicon_file, "rt") as file:
        dico= {}
        prot_id = ""
        for ligne in file:
            if ligne.startswith(">"):
                prot_id = ligne[1:].split()[0]
                dico[prot_id] = ""
            else:
                dico[prot_id] += ligne.strip()
        for i in dico:
            sequence = dico[i]
            if len(sequence) >= minseqlen:
                yield sequence


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    dico={}
    list_seq=[]
    for i in read_fasta(amplicon_file, minseqlen) :
        list_seq.append(i)
    list_index=set(list_seq)
    for i in list_index:
        dico[i]=list_seq.count(i)    
    #print (dico)
    #print(new_dico)
    new_dico=dict(sorted(dico.items(), key = lambda x: x[1], reverse = True))
    #print(new_dico)
    for i,j in new_dico.items():
        if j >= mincount:
            yield[i, j]


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """"""
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
    #print(sequence)
    for kmer in cut_kmer(sequence, kmer_size):
        #print (kmer)
        if kmer in kmer_dict:
            kmer_dict[kmer]+=[id_seq]
        else:
            kmer_dict[kmer] = [id_seq]
    #print (kmer_dict)
    return (kmer_dict)


def chimera_removal(amplicon_file , minseqlen ,mincount , chunk_size , kmer_size ):

    list_sequences = []
    list_occurences = []
    for seq_occ in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        list_sequences.append(seq_occ[0])
        list_occurences.append(seq_occ[1])
    chunks = []
    kmer_dict = {}
    for id, seq in enumerate(list_sequences):
        chunks.append(get_chunks(seq, chunk_size))
        kmer_dict = get_unique_kmer(kmer_dict, seq, id, kmer_size)
    parent = []
    for chunk in chunks:
        for seq2 in chunk:
            parent.append(search_mates(kmer_dict, seq2, kmer_size))

    parent_seq = common(parent[0], parent[1])
    chim_id = []
    chunk_list = [get_chunks(list_sequences[parent_seq[0]], chunk_size)]
    chunk_list += [get_chunks(list_sequences[parent_seq[1]], chunk_size)]

    for id, seq in enumerate(list_sequences):
        if seq not in parent_seq:
            chimera = get_chunks(seq, chunk_size)
            matrix = [[] for chunk in range(len(chimera))]
            for id2, chunk2 in enumerate(chimera):
                for chunk in range(len(chunk_list)):
                    matrix[id2].append(get_identity(nw.global_align(chunk2,
                        chunk_list[chunk], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))))
                #print (id2)
            if detect_chimera(matrix):
                chim_id.append(index)

    for index, seq in enumerate(sequences):
        if index not in chim_id:
            yield [seq, occurences[index]]



def detect_chimera(perc_identity_matrix):
    list_sd=[]
    list_proche=[]
    for i in perc_identity_matrix:
        std=round(statistics.stdev(i),1)
        list_sd.append(std)
        if i[0]>i[1]:
            list_proche.append(1)
        else :
            list_proche.append(0)
    moy= round(statistics.mean(list_sd),1)
    if moy >5:
        c1=list_proche.count(1)
        c0=list_proche.count(0)
        if c1 == 0 or c0 ==0 :
            return False
        else :
           return True
    else:
        return False


def search_mates(kmer_dict, sequence, kmer_size):
    list_kmer = []
    for kmer in cut_kmer(sequence, kmer_size):
        list_kmer.append(kmer)
        if kmer in kmer_dict:
            list_kmer += kmer_dict[kmer]
    best_mates = Counter(list_kmer).most_common(2)
    return [seq_id for seq_id, count in best_mates]


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    #normalement il faut changer a ce niveau dereplication_fulllength et metre chimera_removal, mais vue que cette derniere n'a pas reussi le test on prefere garder dereplication_fulllength
    list_otu=[]
    list_occ = [seq for seq in dereplication_fulllength(amplicon_file, minseqlen, mincount)]
    list_otu.append([list_occ[0][0], list_occ[0][1]])
    for i in range(1,len(list_occ)):
        for j in range(0,len(list_occ)):
            #print(j)
            align=nw.global_align(list_occ[i][0],list_occ[j][0], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
            align=list(align)
            #print(align)
            sim=get_identity(align)
            #print(sim)
            if (sim <97):
                list_otu.append([list_occ[i][0], list_occ[i][1]])
    return (list_otu)



def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    with open(output_file, "w") as f:
        count = 1
        for i, (seq, occ) in enumerate(OTU_list):
            f.write(f">OTU_{count} occurrence:{occ}\n")
            f.write(f"{fill(seq)}")
            f.write("\n")
            count += 1

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # notre programme ici
    chimerafree = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)
    print (chimerafree)

if __name__ == '__main__':
    main()
