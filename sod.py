#!/usr/bin/env python3

'''
SoD: sample and mutate subsequences from a supplied reference nucleotide sequence
USAGE ./sod.py -n 10000 -l 250 -s 0.1 -i 0.05 -d 0.05 tests/hxb2.fasta > outfile.fa
HELP sod.py -h

'''
__author__ = 'Bede Constantinides'
__license__ = 'GPL v3'
__version__ = '0.1.0'

import sys
import argh
import numpy
import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet


def rand_len(indel_ext_rate):
    return int(numpy.random.geometric(1-indel_ext_rate, size=1)[0])


def rand_seq(len):
    return ''.join([random.choice('actg') for base in range(len)])


def rand_base_except(base):
    return random.choice('actg'.replace(base, ''))


def mutate(len_reads, sub_rate, ins_rate, del_rate, indel_ext_rate):
    '''
    Returns masks for each mutation type, as well as the required template length
    '''
    temp_i, read_i = 0, 0
    ins_mask, del_mask, sub_mask = [], [], []
    while read_i < len_reads:
        insertion = rand_len(indel_ext_rate) if random.random() < ins_rate else 0
        deletion = rand_len(indel_ext_rate) if random.random() < del_rate else 0
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
    mut_len = temp_i
    return (ins_mask, del_mask, sub_mask), mut_len


def simulate_read(ref_f, ref_rc, ref_len, mut_len, ids_masks, i):
    '''
    Returns simulated SeqRecord
    '''
    ins_mask, del_mask, sub_mask = ids_masks
    ins_count, del_count, sub_count = (sum(map(bool, mask)) for mask in ids_masks)
    direction = 1 if random.getrandbits(1) else 0
    direction_fmt = 'f' if direction else 'r'
    start_pos = random.randint(0, ref_len-mut_len)
    ref = str(ref_f) if direction else str(ref_rc)
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
    start_pos = start_pos if direction else ref_len - start_pos
    end_pos = start_pos + mut_len
    read_header = 'read{}_{}{}_s{}_i{}_d{}'.format(i, direction_fmt, start_pos,
                                                 sub_count, ins_count, del_count)
    record = SeqRecord(Seq(read, DNAAlphabet()), id=read_header, description='')
    return record


def simulate_reads(path_to_ref, n_reads=1, len_reads=250, sub_rate=0.0, ins_rate=0.0,
                   del_rate=0.0, indel_ext_rate=0.5, fastq=False, uppercase=False,):
    '''
    Returns list of simulated SeqRecords
    '''
    assert max(ins_rate, del_rate, sub_rate) <= 0.5

    ref_f = str(SeqIO.read(path_to_ref, 'fasta').seq)
    ref_rc = Seq(ref_f, DNAAlphabet()).reverse_complement()

    records = []
    for i in range(n_reads):
        ids_masks, mut_len = mutate(len_reads, sub_rate, ins_rate, del_rate, indel_ext_rate)
        record = simulate_read(ref_f, ref_rc, len(ref_f), mut_len, ids_masks, i+1)
        if fastq:
            record.letter_annotations['phred_quality'] = [30] * len(record)
        if uppercase:
            record.seq = record.seq.upper()
        records.append(record)
    return records


def main(path_to_ref, n_reads=1, len_reads=250, sub_rate=0.0, ins_rate=0.0, del_rate=0.0,
         indel_ext_rate=0.5, fastq=False, uppercase=False):
    records = simulate_reads(path_to_ref, n_reads, len_reads, sub_rate,
                             ins_rate, del_rate, indel_ext_rate, fastq, uppercase)
    if fastq:
        for record in records:
            record.letter_annotations['phred_quality'] = [30] * len(record)
        SeqIO.write(records, sys.stdout, 'fastq')
    else:
        SeqIO.write(records, sys.stdout, 'fasta')

if __name__ == '__main__':
    argh.dispatch_command(main)
