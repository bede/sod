#!/usr/bin/env python3

# SoD: sample and mutate subsequences from a supplied reference nucleotide sequence
# USAGE ./sod.py -n 10000 -l 250 -s 0.1 -i 0.05 -d 0.05 tests/hxb2.fasta > outfile.fa

__author__ = 'Bede Constantinides'
__license__ = 'GPL v3'
__version__ = '1.0.1'

import sys
import argh
import numpy
import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet


def rand_len():
    return int(numpy.random.geometric(0.5, size=1)[0])


def rand_seq(len):
    return ''.join([random.choice('actg') for base in range(len)])


def rand_base_except(base):
    return random.choice('actg'.replace(base, ''))


def mutate(len_reads, sub_rate, ins_rate, del_rate):
    '''
    Returns masks for each mutation type, as well as the required template length
    '''
    temp_i, read_i = 0, 0
    ins_mask, del_mask, sub_mask = [], [], []
    while read_i < len_reads:
        insertion = rand_len() if random.random() < ins_rate else 0
        deletion = rand_len() if random.random() < del_rate else 0
        substitution = 1 if random.random() < sub_rate else 0
        if insertion:
            if read_i + insertion + 1 > len_reads:
                insertion -= 1 + insertion + read_i - len_reads
            read_i += insertion
        if deletion:
            if read_i + deletion + 1 > len_reads:
                deletion -= 1 + deletion + read_i - len_reads
            temp_i += deletion
        temp_i += 1
        read_i += 1
        ins_mask.append(insertion)
        del_mask.append(deletion)
        sub_mask.append(substitution)
    mut_read_len = temp_i
    return (ins_mask, del_mask, sub_mask), mut_read_len


def simulate(ref_fwd, ref_rev_cmp, ref_len, mut_read_len, ids_masks):
    '''
    Returns SeqRecord of simulated sequence
    '''
    ins_mask, del_mask, sub_mask = ids_masks
    start_pos = random.randint(0, ref_len-mut_read_len)
    direction = 'fwd' if random.getrandbits(1) else 'rev'
    ref = str(ref_fwd) if direction == 'fwd' else str(ref_rev_cmp)
    ref_i = start_pos
    read = ''
    for insertion, deletion, substitution in zip(ins_mask, del_mask, sub_mask):
        base = rand_base_except(ref[ref_i]) if substitution else ref[ref_i]
        if insertion:
            read += base + rand_seq(insertion)
        if deletion:
            ref_i += deletion
        if not insertion:
            read += base
        ref_i += 1
    start_pos = start_pos if direction == 'fwd' else ref_len - start_pos
    end_pos = start_pos + mut_read_len if direction == 'fwd' else start_pos - mut_read_len
    read_header = direction + '_s' + str(start_pos) + '_e' + str(end_pos)
    record = SeqRecord(Seq(read, DNAAlphabet()), id=read_header, description='')
    return record


def main(path_to_ref, n_reads=1, len_reads=250, sub_rate=0.0, ins_rate=0.0, del_rate=0.0, fastq=False):
    ref_fwd = str(SeqIO.read(path_to_ref, 'fasta').seq)
    ref_rev_cmp = Seq(ref_fwd, DNAAlphabet()).reverse_complement()
    ref_len = len(ref_fwd)
    
    records = []
    for i in range(n_reads):
        ids_masks, mut_read_len = mutate(len_reads, sub_rate, ins_rate, del_rate)
        record = simulate(ref_fwd, ref_rev_cmp, ref_len, mut_read_len, ids_masks)
        records.append(record)
        assert len(record.seq) == len_reads

    if fastq:
        for record in records:
            record.letter_annotations['phred_quality'] = [30] * len(record)
        SeqIO.write(records, sys.stdout, 'fastq')
    else:
        SeqIO.write(records, sys.stdout, 'fasta')


argh.dispatch_command(main)
